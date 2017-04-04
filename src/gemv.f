        program gemv1
              use mpi
              use example
              implicit none
        
              integer :: n, nb    ! problem size and block size
              integer :: myArows, myAcols   ! size of local subset of global matrix
              integer :: myXrows, myXcols   ! size of local subset of global vector
              integer :: i,j, myi, myj
              real, dimension(:,:), allocatable :: myA,myX,myY
              integer :: ierr
        
              integer, external :: numroc   ! blacs routine
              integer :: me, procs, icontxt, prow, pcol, myrow, mycol  ! blacs data
              integer :: info    ! scalapack return value
              integer, dimension(9)   :: ides_a, ides_x, ides_y ! matrix descriptors
              integer, dimension(2) :: dims
              real :: error, globerror
        
        ! try to call subroutine from my module
              call printerup
        ! Initialize blacs processor grid
        
              call blacs_pinfo   (me,procs)
        
        ! create as square as possible a grid of processors
        
              dims = 0
              call MPI_Dims_create(procs, 2, dims, ierr)
              prow = dims(1)
              pcol = dims(2)
        
        ! create the BLACS context
        
              call blacs_get     (0, 0, icontxt)
              call blacs_gridinit(icontxt, 'R', prow, pcol)
              call blacs_gridinfo(icontxt, prow, pcol, myrow, mycol)
        
        ! Construct local arrays
        ! Global structure:  matrix A of n rows and n columns
        
              n = int(2500.*sqrt(dble(procs)))
        
        ! blocksize - a free parameter.
        
              nb = 100
        
        ! how big is "my" chunk of matrix A?
        
              myArows = numroc(n, nb, myrow, 0, prow)
              myAcols = numroc(n, nb, mycol, 0, pcol)
        
        ! how big is "my" chunk of vector x?
        
              myXrows = numroc(n, nb, myrow, 0, prow)
              myXcols = 1
        
        ! Initialize local arrays    
        
              allocate(myA(myArows,myAcols)) 
              allocate(myX(myXrows,myXcols)) 
              allocate(myY(myXrows,myXcols)) 
        
              myA = 0.
              do myj=1,myAcols
                  ! get global index from local index
                  call l2g(myj,mycol,n,pcol,nb,j)
                  do myi=1,myArows
                      ! get global index from local index
                      call l2g(myi,myrow,n,prow,nb,i)
                      if (i == j) myA(myi,myj) = i
                  enddo
              enddo
        
              myX = 0.
              call l2g(1,mycol,n,pcol,nb,j)
              if (j == 1) then
                  do myi=1,myXrows
                      call l2g(myi,myrow,n,prow,nb,i)
                      myX(myi,1) = 1.
                  enddo
              endif
        
              myY = 0.
        
        ! Prepare array descriptors for ScaLAPACK 
        
        call descinit(ides_a, n, n, nb, nb, 0, 0, icontxt, myArows,info)
        call descinit(ides_x, n, 1, nb, nb, 0, 0, icontxt, myXrows,info)
        call descinit(ides_y, n, 1, nb, nb, 0, 0, icontxt, myXrows,info)
        
        ! Call ScaLAPACK library routine
        
        call psgemv('N',n,n,1.,mya,1,1,ides_a,myx,1,1,
     & ides_x,1,0.,myy,1,1,ides_y,1)
              if (me == 0) then
                if (info /= 0) then
                     print *, 'Error -- info = ', info
                endif
              endif
        
        ! Deallocate the local arrays
        
              deallocate(myA, myX)
        
        ! Print results - Y should be 1,2,3...
        
              error = 0.
              call l2g(1,mycol,n,pcol,nb,j)
              if (j == 1) then
                  do myi=1,myXrows
                      call l2g(myi,myrow,n,prow,nb,i)
                      error = error + (myY(myi,1)-1.*i)**2.
                      !print *, myi, myY(myi,1)
                  enddo
              endif
        
              !    MPI_REDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, ROOT, COMM, IERROR)
              call MPI_Reduce(error, globerror, 1, MPI_REAL, MPI_SUM,
     & 0, MPI_COMM_WORLD, ierr)
        
              if (me == 0) then
                print *,'Y l2 error = ', sqrt(globerror/n)
              endif
        
              deallocate(myY)
        
        ! End blacs for processors that are used
        
              call blacs_gridexit(icontxt)
              call blacs_exit(0)
        
        contains
        
        ! convert global index to local index in block-cyclic distribution
        
           subroutine g2l(i,n,np,nb,p,il)
        
           implicit none
           integer, intent(in) :: i    ! global array index, input
           integer, intent(in) :: n    ! global array dimension, input
           integer, intent(in) :: np   ! processor array dimension, input
           integer, intent(in) :: nb   ! block size, input
           integer, intent(out):: p    ! processor array index, output
           integer, intent(out):: il   ! local array index, output
           integer :: im1   
        
           im1 = i-1
           p   = mod((im1/nb),np)
           il  = (im1/(np*nb))*nb + mod(im1,nb) + 1
        
           return
           end subroutine g2l
        
        ! convert local index to global index in block-cyclic distribution
        
           subroutine l2g(il,p,n,np,nb,i)
        
           implicit none
           integer :: il   ! local array index, input
           integer :: p    ! processor array index, input
           integer :: n    ! global array dimension, input
           integer :: np   ! processor array dimension, input
           integer :: nb   ! block size, input
           integer :: i    ! global array index, output
           integer :: ilm1   
        
           ilm1 = il-1
           i    = (((ilm1/nb) * np) + p)*nb + mod(ilm1,nb) + 1
        
           return
           end subroutine l2g
        
        end program gemv1
