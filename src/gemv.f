        program gemv1
              use mpi
              !use example
              implicit none
        
              integer :: n, nb    ! problem size and block size
              integer :: myArows, myAcols   ! size of local subset of global matrix
              integer :: myXrows, myXcols   ! size of local subset of global vector
              integer :: i,j, myi, myj
              double precision, dimension(:,:), allocatable :: myA
              integer, dimension(:), allocatable :: myX
              integer :: ierr
              character (len=10) :: arg
              character (len=5), parameter :: fmtstr = "e11.5"
        
              integer, external :: numroc   ! blacs routine
              integer :: me, procs, icontxt, prow, pcol, myrow, mycol  ! blacs data
              integer :: info    ! scalapack return value
              integer, dimension(9)  :: ides_a, ides_x ! matrix descriptors
              integer, dimension(2) :: dims
              double precision :: det, globdet ! determinant and global
              double precision :: starttime, laptime, stoptime ! for Wtime
              
              integer, allocatable :: seed(:) ! random seed array
              integer :: randsize, r

        ! get problem size from input 
              call get_command_argument(1, arg)
              READ (arg(:),'(i10)') n

        ! start first timer
              starttime = MPI_Wtime()
              laptime = starttime

        ! Initialize blacs processor grid
        
              call blacs_pinfo  (me,procs)

        ! blocksize - a free parameter.
              nb = 10

        ! create as square as possible a grid of processors
        
              dims = 0
              call MPI_Dims_create(procs, 2, dims, ierr)
              prow = dims(1)
              pcol = dims(2)
              !print *, "prow",prow,"pcol",pcol
        
        ! create the BLACS context
        
              call blacs_get     (0, 0, icontxt)
              call blacs_gridinit(icontxt, 'R', prow, pcol)
              call blacs_gridinfo(icontxt, prow, pcol, myrow, mycol)
        
        ! Construct local arrays
        
        ! how big is "my" chunk of matrix A?
        
              myArows = numroc(n, nb, myrow, 0, prow)
              myAcols = numroc(n, nb, mycol, 0, pcol)
              !print *, "myrows",myarows,"mycols",myacols
        
        ! how big is "my" chunk of vector x?
        
              myXrows = numroc(n, nb, myrow, 0, prow) + nb
              myXcols = 1
        
        if (me == 0) write (*,'(I3A1I4A1)',advance="no"),
     & procs,',',n,','

        ! Initialize local arrays    
        
              allocate(myA(myArows,myAcols)) 
              allocate(myX(myXrows)) 
        
            ! get random seed
              call random_seed(size=randsize)
              allocate(seed(randsize))
              do r=1,randsize
                  call system_clock(seed(r))
              enddo
              call random_seed(put=seed)

              call random_number(myA)
              myA = myA - 0.5d+0
              !myA = 0.d+0
              !do myj=1,myAcols
                  ! get global index from local index
              !    call l2g(myj,mycol,n,pcol,nb,j)
              !    do myi=1,myArows
                      ! get global index from local index
              !        call l2g(myi,myrow,n,prow,nb,i)
              !        if (i <= j) myA(myi,myj) = 1.d+0
              !    enddo
              !enddo

              !print matrix out
            if (me == 0 .AND. n < 10) then
              do myi=1,myArows
                do myj=1,myAcols 
                      write (*,"(f8.3)",advance="no"),myA(myi,myj)
                enddo
                print *,""
              enddo
            endif


              myX = 0

        
        ! Prepare array descriptors for ScaLAPACK 
            call descinit( ides_a, n, n, nb, nb, 0, 0, icontxt, 
     & myArows, info )
            call descinit( ides_x, n, 1, nb, nb, 0, 0, icontxt,
     & myXrows, info )
        
            if (me == 0) then
                stoptime = MPI_Wtime()
                write (*,"("//fmtstr//"a1)",advance="no"),
     & stoptime-laptime,','
                !print *, 'Setup time: ',stoptime-laptime
                laptime = stoptime
            endif

        ! Call ScaLAPACK library routine
        !  SUBROUTINE PDGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
        call pdgetrf(n,n,myA,1,1,ides_a,myX,info)

              if (me == 0) then
                if (info /= 0) then
                     print *, 'Error -- info = ', info
                endif
              endif
        
              !print matrix out
            if (me == 0 .AND. n < 10) then
              print *,'Decomposition'
              do myi=1,myArows
                do myj=1,myAcols 
                      write (*,"(f8.3)",advance="no"),myA(myi,myj)
                enddo
                print *,""
              enddo
            endif
            if (me == 0) then
                stoptime = MPI_Wtime()
                write (*,"("//fmtstr//"a1)",advance="no"),
     & stoptime-laptime,','
                !print *, 'Matrix decomp time: ',stoptime-laptime
                laptime = stoptime
            endif

        
        ! Calculate determinant based on diagonal

            det = 1.d+0
            do myj=1,myAcols
                call l2g(myj,mycol,n,pcol,nb,j)
                do myi=1,myArows
                    call l2g(myi,myrow,n,prow,nb,i)
                    if (j == i) det = det * myA(myi,myj)
                enddo
            enddo

            !print *, 'Local determinant ', det, 'proc', me
            if (me == 0) then
                stoptime = MPI_Wtime()
                write (*,"("//fmtstr//"a1)",advance="no"),
     & stoptime-laptime,','
                !print *, 'Determinant calc time: ',stoptime-laptime
                laptime = stoptime
            endif
        
        
              !    MPI_REDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, ROOT, COMM, IERROR)
              call MPI_Reduce(det, globdet, 1, MPI_DOUBLE, MPI_PROD, 0,
     &  MPI_COMM_WORLD, ierr)
        
            if (me == 0) then
                stoptime = MPI_Wtime()
                write (*,"("//fmtstr//"a1)",advance="no"),
     & stoptime-laptime,','
                !print *, 'Reduction time: ',stoptime-laptime
                laptime = stoptime
            endif
        ! Deallocate X
        
              deallocate(myA,myX)

        ! End blacs for processors that are used
        
              call blacs_gridexit(icontxt)
              call blacs_exit(0)
        
            if (me == 0) then
                stoptime = MPI_Wtime()
                write (*,"("//fmtstr//")"),stoptime - starttime
                !print *, 'Total time: ',stoptime - starttime
            endif

            if (me == 0) then
              print *,'Determinant = ', globdet
            endif
        

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
