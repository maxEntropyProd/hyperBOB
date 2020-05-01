   ! This routine is basically an MPI wrapper to call BOBYQA using initial points selected from
   ! a random latin hypercube.  MPI is being used to run BOBYQA in parallel from different starting points
   !
   ! V1.0 1-Feb-2019
   ! V1.1 10-Feb-2019
   !  -  Allow saving of fVal and xMat for all processes (default it no saving)
   ! V1.2 12-Feb-2019
   !  -  Added routine hyperBOB_initialize and hyperBOB_cleanUp, which now must be called before and after calls
   !     to hyperBOB.  This allows xMat and fVal to be accessable, and these names have been changed to
   !     xMat_hyperBOB and fVal_hyperBOB
   ! V1.3 13-Feb-2019
   !  -  This version allows the row dimension of xMat_hyperBOB to be larger than needed, which allows 
   !     post modifcations that may need more space.
   ! V1.4 13-Jan-2020
   !  -  Bug corrected in setting rhoBMin.
   module hyperBOB_module
      ! This module setups up variables needed by hyperBOB
      use kind_module ! sets up double precision and found in BOBYQA
      use bobyqa_module ! this is the F90 version found on github: https://github.com/jacobwilliams/PowellOpt
      use mpi
      implicit none
      real(wp), allocatable:: xMat_hyperBOB(:,:) ! contains the initial search locations from the random hypercube, but only rank 0 has the full matrix on return.
      real(wp), allocatable:: fVal_hyperBOB(:)   ! Stores the function value for all processes      
      ! ****** NOTE, the function must be defined EXACTLY as shown by the interface here.  *******
      ! that is, don't instead declare as "real(wp) x(n)", etc as that will cause problems
      abstract interface
         subroutine func (n, x, f)  ! calfun interface
            import :: wp
            implicit none
            integer,intent(in)               :: n
            real(wp),dimension(:),intent(in) :: x
            real(wp),intent(out)             :: f
         end subroutine func
      end interface
      
   contains
      
      subroutine hyperBOB_initialize(nRows) 
         integer, intent(in) :: nRows  ! 1st dimension of xMat_hyperBOB
         ! local declarations
         integer mpiErr, noProc
         
         ! get the number of processes that can be used
         call MPI_COMM_SIZE(MPI_COMM_WORLD, noProc, mpiErr)
         ! Allocate space
         allocate( xMat_hyperBOB(nRows, noProc) ) ! Used for intial guess and function minimum
         allocate( fVal_hyperBOB(noProc) ) ! Function value at min.
         return
      end subroutine hyperBOB_initialize      
            
      subroutine hyperBOB_cleanUp()
         ! once all calls to hyperBOB are completed, this routine should be called to deallocated space
         ! or if a call to hyperBOB requires redimensioning problem, but then hyperBOB_initialize needs
         ! to be recalled.
         deallocate( xMat_hyperBOB ) 
         deallocate( fVal_hyperBOB ) 
         return
      end subroutine hyperBOB_cleanUp
      
      subroutine hyperBOB(nx, npt, x, xL, xU, rhoBMin, rhoEMax, iprint, maxFun, maxTime, calfun, seed, fMin, err)     
         ! This is the only routine to call
         integer,  intent(in)   :: nx       ! dimension of optimization space
         integer,  intent(in)   :: npt      ! see BOBYQA
         real(wp), intent(out)  :: x(nx)    ! best solution obtained.
         real(wp), intent(in)   :: xL(nx)   ! lower bound on x
         real(wp), intent(in)   :: xU(nx)   ! upper bound on x
         real(wp), intent(inOut):: rhoBMin  ! see BOBYQA; however, that actual value of rhoBMin used depends on the number of processes.  
                                            ! If there are a large number of processes used, the rhoBMin will be reduced to search a local space in the latin hypercube.
         real(wp), intent(inOut):: rhoEMax  ! see BOBYQA.  rhoEMax will be used only if less that 0.1*rhoBMin, where rohBMin is the calculated value
         integer,  intent(in)   :: iprint   ! See BOBYQA.  Typically set to 0 for no output.
         integer,  intent(in)   :: maxFun   ! maximum function evaluations
         real(wp), intent(in)   :: maxTime  ! maximum number of hours to run (currently not implemented)
         procedure (func)       :: calfun   ! Name of objective function passed to BOBYQA
         integer,  intent(inout):: seed     ! Used to start the random num gen. for hypercube.  Should be > 0
         real(wp), intent(out)  :: fMin     ! The minimum function value found.
         integer,  intent(out)  :: err      ! if not 0, then an error and error occured 
                                            ! err = 1 if first dimension of xMat_hyperBOB < nx.

         ! local declarations
         integer nBOBs ! the number of BOBYQA processes run in parallel
         integer nRows ! the first dimension of xMat_hyperBOB, which must be >= nx
         integer i, myRank, mpiErr
         integer, parameter:: master = 0 ! set which process will be the master (typically 0)
         real(wp) temp
         integer bestBOB ! process that finds the best value
         
         err = 0 ! No errors yet
         ! get rank of processes. Note, this assumes hyperBOB_initialize and MPI_INIT( mpierr ) has already been called
         call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, mpiErr) ! get rank of this process in world
         ! Get the number of processes from xMat_hyperBOB that has been previous allocated with call to hyperBOB_initialize 
         nBOBs = size(xMat_hyperBOB, dim=2)
         nRows = size(xMat_hyperBOB, dim=1)
         if (nRows < nx) then
            err = 1
            return
         end if         
         
         ! Generate the starting search locations with nX space sampled as a latin hypercube.  
         ! While latin_random generates a random matrix, if given the same seed, the same matrix will be
         ! generated.  Therefore, each process can do this separatedly and there is no need for just the master followed by a broadcast
         call latin_random (nx, nBOBs, seed, xMat_hyperBOB(1:nx,1:nBOBs))
         
         ! Insure xMat_hyperBOB is within xL and xU by scaling the columns.
         do i=1, nBOBs
            xMat_hyperBOB(1:nx,i) = (xU(1:nx) - xL(1:nx))*xMat_hyperBOB(1:nx,i) + xL(1:nx)
         end do
         
         ! Determine best value to use for rhoBMin and adjust rhoEMax if necessary
         ! BOBYQA issues an error if rhoB > (xU(j) - xL(j))/2, so make it equal to that, but for each subzone
         temp = minval(xU-xL)/2.*(1./real(nBOBs))**(1./real(nx)) ! Each process should search a subvolume that is 1/nx of the full domain (xL, xU)
         if (temp < rhoBMin) rhoBMin = temp ! This allows the user to specify a smaller rhoB, but not larger
         if (rhoEmax > 0.1*rhoBMin) rhoEmax = 0.1*rhoBMin ! Need to insure the final radius is smaller than the starting one. 
         ! For large dimensional systems, like nx > 100 (which is common), even for 10,000 process, rhoB = 1/2(0.91), so still pretty much the whole domain.
         
         ! Now have every process runs it's own BOBYQA, but with a different starting point
         call bobyqa (nx, npt, xMat_hyperBOB(1:nx,myRank+1), xL, xU, rhoBMin, rhoEmax, iprint, maxfun, calfun)
         
         ! Select the best solution and return that and it's stats
         ! call calfun to deterime function value found for each process
         call calfun (nx, xMat_hyperBOB(1:nx,myRank+1), fVal_hyperBOB(myRank+1))
         
         ! Broadcast fVal_hyperBOB to all processes (each just has one element of fVal_hyperBOB)
         do i = 1, nBOBs     
            ! In MPI-BCAST both the sending and receiving processes must specify the same root.
            call MPI_BCAST(fVal_hyperBOB(i), 1, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpierr)
         end do
         ! Determine which process has the best value.  BOBYQA finds minimums
         bestBOB = minloc(fVal_hyperBOB(1:nBOBs),1) ! This will return the first one found if there are more than one min.
         fMin = fVal_hyperBOB(bestBOB) ! Minium function value found for all processes
         x(1:nx) = xMat_hyperBOB(1:nx,myRank+1) ! temporary storing of solution to avoid MPI Buffers must not be aliased in Gather below
         ! Gather solutions from all processes in rank 0 for xMat (all processes have fVal already)
         ! NOTE complete xMat is only available in rank 0, all other processes only have their local solution and bestBOB (copied below)
         call MPI_GATHER(   x, nx, MPI_DOUBLE_PRECISION, &
                         xMat_hyperBOB(1:nx,1:nBOBs), nx, MPI_DOUBLE_PRECISION, &
                          0, MPI_COMM_WORLD, mpierr)
         ! Now broadcast the best solution to all processes, no need to broadcast the others
         call MPI_BCAST(xMat_hyperBOB(1:nx,bestBOB), nx, MPI_DOUBLE_PRECISION, bestBOB-1, MPI_COMM_WORLD, mpierr)
         x(1:nx) = xMat_hyperBOB(1:nx,bestBOB) ! All processes will return this same value
         return
      end subroutine hyperBOB      
      
   end module hyperBOB_module