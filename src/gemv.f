        program gemv1
              use mpi
              !use example
              implicit none
        
              integer :: n, nb    ! problem size and block size
              integer :: myArows, myAcols   ! size of local subset of global matrix
              integer :: myXrows, myXcols   ! size of local subset of global vector
              integer :: i,j, myi, myj
              double precision, dimension(:,:), allocatable :: myA, myX
              integer :: ierr
              character (len=10) :: arg
              character (len=5), parameter :: fmtstr = "e11.5"
              character (len=10) :: fname
        
              integer, external :: numroc   ! blacs routine
              integer :: me, procs, icontxt, prow, pcol, myrow, mycol  ! blacs data
              integer :: info    ! scalapack return value
              integer, dimension(9)  :: ides_a, ides_x ! matrix descriptors
              integer, dimension(2) :: dims
              double precision :: det, globdet ! determinant and global
              double precision :: starttime, laptime, stoptime ! for Wtime
              
              integer, allocatable :: seed(:) ! random seed array
              integer :: randsize

        ! get problem size from input
              call get_command_argument(1, arg)
              READ (arg(:),'(i10)') n

        ! start first timer
              starttime = MPI_Wtime()
              laptime = starttime

        ! Initialize blacs processor grid
        
              call blacs_pinfo  (me,procs)

        ! use file to write extra output
        if (n < 10) then
              write (fname, "(A4I1)") "proc", me
              OPEN (unit=7,file= trim(fname))
        endif

        ! blocksize - a free parameter.
              nb = 4

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
        
              myXrows = numroc(n, nb, myrow, 0, prow)
              myXcols = numroc(n, nb, mycol, 0, pcol)
        
        if (me == 0) write (*,'(I3A1I4A1)',advance="no"),
     & procs,',',n,','

        ! Initialize local arrays
        
              allocate(myA(myArows,myAcols))
              allocate(myX(myXrows,myXcols))
        
            ! get random seed
              call random_seed(size=randsize)
              allocate(seed(randsize))
         OPEN(89,FILE='/dev/urandom',ACCESS='stream',FORM='UNFORMATTED')
              READ(89) seed
              CLOSE(89)
              call random_seed(put=seed)
              call random_number(myA)
              myA = myA - 0.5d+0
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
            if (n < 10) then
              do myi=1,myArows
                do myj=1,myAcols
                      write (7,"(f8.3)",advance="no"),myA(myi,myj)
                enddo
                write (7,"(A1)") ' '
              enddo
            endif


             ! myX = 0.d0

        
        ! Prepare array descriptors for ScaLAPACK
            call descinit( ides_a, n, n, nb, nb, 0, 0, icontxt,
     & myArows, info )
            call descinit( ides_x, n, n, nb, nb, 0, 0, icontxt,
     & myXrows, info )
        
            if (me == 0) then
                stoptime = MPI_Wtime()
                write (*,"("//fmtstr//"a1)",advance="no"),
     & stoptime-laptime,','
                !print *, 'Setup time: ',stoptime-laptime
                laptime = stoptime
            endif
        ! Call ScaLAPACK library routine
        !void   pdgemm_ (F_CHAR_T TRANSA, F_CHAR_T TRANSB, int *M, int
        !*N, int *K, double *ALPHA, double *A, int *IA, int *JA, int
        !*DESCA, double *B, int *IB, int *JB, int *DESCB, double *BETA,
        !double *C, int *IC, int *JC, int *DESCC)

        !multiplies A by A^T and stores in X

        call pdgemm('N','T',n,n,n,1.d+0,myA,1,1,ides_a,myA,
     & 1,1,ides_a,0.d+0,myX,1,1,ides_x)

              if (me == 0) then
                if (info /= 0) then
                     print *, 'Error -- info = ', info
                endif
              endif
        
              if (n < 10) then
              write (7,"(A6)"),'Result'
              !print matrix out
              do myi=1,myXrows
                do myj=1,myXcols
                      write (7,"(f8.3)",advance="no"),myX(myi,myj)
                enddo
                write (7,"(A1)") ' '
              enddo
              endif

        ! Deallocate A
        
              deallocate(myA)

            if (me == 0) then
                stoptime = MPI_Wtime()
                write (*,"("//fmtstr//"a1)",advance="no"),
     & stoptime-laptime,','
                !print *, 'Matrix mult time: ',stoptime-laptime
                laptime = stoptime
            endif

        ! Now do Cholesky decomposition
        ! PDPOTRF (UPLO, N, A, IA, JA, DESCA, INFO)
        call pdpotrf( 'L', n, myX, 1, 1, ides_x, info )

        !if info is greater than 0 this means that determinant is zero?

              if (me == 0) then
                if (info /= 0) then
                     print *, 'Error -- info = ', info
                endif
              endif
        
              if (n < 10) then
              write (7,"(A13)"),'Decomposition'
              !print matrix out
              do myi=1,myXrows
                do myj=1,myXcols
                      write (7,"(f8.3)",advance="no"),myX(myi,myj)
                enddo
                write (7,"(A1)") ' '
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
            do myj=1,myXcols
                call l2g(myj,mycol,n,pcol,nb,j)
                do myi=1,myXrows
                    call l2g(myi,myrow,n,prow,nb,i)
                    if (j == i) det = det * myX(myi,myj)
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
                print *,'Determinant = ', globdet
              endif
        
            if (me == 0) then
                stoptime = MPI_Wtime()
                write (*,"("//fmtstr//"a1)",advance="no"),
     & stoptime-laptime,','
                !print *, 'Reduction time: ',stoptime-laptime
                laptime = stoptime
            endif
        ! Deallocate X
        
              deallocate(myX)

        ! End blacs for processors that are used
        
              call blacs_gridexit(icontxt)
              call blacs_exit(0)
        
            if (me == 0) then
                stoptime = MPI_Wtime()
                write (*,"("//fmtstr//")"),stoptime - starttime
                !print *, 'Total time: ',stoptime - starttime
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
