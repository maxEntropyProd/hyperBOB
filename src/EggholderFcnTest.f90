   ! This is a quick test of hyperBOB using the Eggholder function
   ! see https://en.wikipedia.org/wiki/Rastrigin_function
   ! This uses Ver 1.3 of hyperBOB which allows xMat_hyperBOB to be bigger than needed.  
   
   program EggholderTest
      use hyperBOB_module
      use mpi
      implicit none
      integer mpierr, myRank, nunit, nBOBs, i
      integer nx, npt, iprint, maxFun, seed, err
      real(wp) rhoBMin, rhoEMax, fMin, maxTime
      real(wp), allocatable:: x(:), xL(:), xU(:)
      real(wp) comp_start, comp_end
      character(len=40) fmtString      
      external EggholderFcn
      
      ! Initialize MPI
      call MPI_INIT( mpierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, mpiErr) ! get rank of this process in world
      if (myRank == 0) then
         write(*,'(/,a)') 'Starting Eggholder Test'
         write(*,'(a,$)') 'Enter maxFun evaluations: '
         read(*,*) maxFun
      end if
      ! boradcast maxFun to all processes
      call MPI_BCAST(maxFun, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      nx = 2
      ! allocate space
      allocate( x(nx), xL(nx), xU(nx) )
      npt = 2*nx + 1
      xL = -512.; xU = 512.
      rhoBMin = 512.
      rhoEMax = 0.001
      iprint = 0
      !maxFun = 2000
      maxTime = 0 ! this is not implemented yet
      seed = 5 ! random number start
      comp_start = MPI_Wtime()
      ! intialize hyperBOB
      call hyperBOB_initialize(nx+5)   ! add 5 extra rows    
      call hyperBOB(nx, npt, x, xL, xU, rhoBMin, rhoEMax, iprint, maxFun, maxTime, EggholderFcn, seed, fMin, err)
      ! Save xMat
      if (myRank == 0) then
         open(newunit=nunit, file='allSolns.dat',status='unknown')
         fmtString = '(i4, #######(1x,g13.5))'
         nBOBs = size(xMat_hyperBOB, dim=2)
         write(fmtString(6:12),'(i7)') nBobs
         write(nunit, fmtString) 0, fVal_hyperBOB(1:nBobs)
         do i=1,nx
            write(nunit, fmtString) i, xMat_hyperBOB(i,1:nBobs)
         end do
         close(unit=nunit)
      end if
      ! deallocate space in hyperBOB
      call hyperBOB_cleanUp()
      
      comp_end = MPI_Wtime()
      
      if(myRank == 0) then
         write(*,'(a,g13.5)') 'fMin = ', fMin
         write(*,'(a,2(1x,g13.5))') 'Solution: ', x(1), x(2)
         write(*,'(a,g13.5)') 'Total time (min): ', (comp_end-comp_start)/60.
      end if
   
      call MPI_FINALIZE(mpierr)
      
      deallocate( x, xL, xU )
      stop
      
   end program EggholderTest
   
   subroutine EggholderFcn(n, xy, f)
      use hyperBOB_module
      ! This function's solution: f(512, 404.2319) = -959.6407
      integer, intent(in)              :: n  ! deminsion of problem
      real(wp),dimension(:),intent(in) :: xy ! Value of x were f is to be evaluated at
      real(wp), intent(out)            :: f  ! function value
      ! local declarations
      integer i
      real(wp) x, y
      x = xy(1)
      y = xy(2)
         
      f = -(y + 47.)*sin(sqrt(abs(x/2. + (y + 47.)))) - x*sin(sqrt(abs(x - (y + 47.))))

      return
   end subroutine EggholderFcn
   