        program gemv1
              use mpi
              use readMatrix
              implicit none
        
              integer :: n, nb    ! problem size and block size
              integer :: myArows, myAcols   ! size of local subset of global matrix
              integer :: myXrows, myXcols   ! size of local subset of global vector
              integer :: i,j, myi, myj
              double precision, dimension(:,:), allocatable :: myA,myX
              double precision, dimension(:,:), allocatable :: myM
              integer :: ierr
              character (len=10) :: arg
              character (*), parameter :: fmtstr = "e11.5"
              character (len=10) :: fout
              character (len=50) :: fname
        
              integer, external :: numroc   ! blacs routine
              integer :: me, procs, icontxt, prow, pcol, myrow, mycol  ! blacs data
              integer :: info    ! scalapack return value
              integer, dimension(9)  :: ides_a, ides_x ! matrix descriptors
              integer, dimension(2) :: dims
              double precision :: det, globdet ! determinant and global
              double precision :: starttime, laptime, stoptime ! for Wtime

              integer :: r,nr,c,nc
              integer :: sendr, sendc, recvr, recvc !for scatter calcs
              
        ! get problem size from input
        if (command_argument_count() > 1) then
              call get_command_argument(1, arg)
              READ (arg(:),'(i10)') n
        else
              print *, 'Need matrix size as parameter'
        endif

        ! start first timer
              starttime = MPI_Wtime()
              laptime = starttime

        ! Initialize blacs processor grid
        
              call blacs_pinfo  (me,procs)

        ! use file to write extra output
        if (n < 17) then
              write (fout, "(A4I1)") "proc", me
              OPEN (unit=7,file= trim(fout))
        endif

        ! create as square as possible a grid of processors
        
              dims = 0
              call MPI_Dims_create(procs, 2, dims, ierr)
              prow = dims(1)
              pcol = dims(2)
              !print *, "prow",prow,"pcol",pcol
        
        ! blocksize - a free parameter.
              nb = 2 !n/prow

        ! create the BLACS context
        
              call blacs_get     (0, 0, icontxt)
              call blacs_gridinit(icontxt, 'R', prow, pcol)
              call blacs_gridinfo(icontxt, prow, pcol, myrow, mycol)
              !print *, 'myrow',myrow,'mycol',mycol
        
        ! Construct local arrays
        
        ! how big is "my" chunk of matrix A?
        
              myArows = numroc(n, nb, myrow, 0, prow)
              myAcols = numroc(n, nb, mycol, 0, pcol)
              !print *, "myrows",myarows,"mycols",myacols
        
        ! how big is "my" chunk of vector x?
        
              myXrows = numroc(n, nb, myrow, 0, prow)
              myXcols = numroc(n, nb, mycol, 0, pcol)
        
        if (me == 0) write (*,'(I3A1I4A1I4A1)',advance="no"),
     & procs,',',n,',',nb,','

        ! Initialize local arrays
        
              allocate(myA(myArows,myAcols))
              allocate(myX(myXrows,myXcols))

        ! Populate the matrix at root
        if (me == 0) then
            allocate(myM(n,n))
            if (command_argument_count() > 2) then
              ! use the filename to read matrix in
              call get_command_argument(2, arg)
              READ (arg(:),'(a50)') fname
              call readM(fname,n,myM)
            else
             ! Fill matrix with rand values between -0.5 and 0.5
              call getRandM(n,n,myM)
            endif

            if (n < 17) then
                write (7,"(A)"),'Matrix on root proc:'
              do myi=1,n
                do myj=1,n
                      write (7,"(f8.3)",advance="no"),myM(myi,myj)
                enddo
                write (7,"(A1)") ' '
              enddo
            endif

        endif

        ! Now scatter to processes
        
        ! Scatter matrix 
        sendr = 0
        recvr = 1
        recvc = 1
        do r=1, n, nb
          sendc = 0
          ! Number of rows to be sent
          ! Is this the last row block?
          nr = nb
          if (n-r+1 < nb) nr = n-r+1
   
          do c=1, n, nb
              ! Number of cols to be sent
              ! Is this the last col block?
              nc = nb
              if (n-c+1 < nb) nc = n-c+1
   
              ! Send a nr-by-nc submatrix to process (sendr, sendc)
              if (me == 0) call dgesd2d(icontxt, nr, nc, myM(r,c),
     &  n, sendr, sendc)
   
              if (myrow == sendr .AND. mycol == sendc) then
                  ! Receive the same data
                  ! The leading dimension of the local matrix is nrows!
                  call dgerv2d(icontxt, nr,nc, myA(recvr,recvc),
     &  myArows, 0, 0)
                  recvc = mod(recvc+nc-1,myAcols) + 1
              endif
   
              sendc = mod(sendc+1,pcol)
          enddo
   
          if (myrow == sendr) recvr = mod(recvr+nr-1,myArows) + 1
          sendr = mod(sendr+1,prow)

        enddo

              !print matrix out
            if (n < 17) then
                write (7,"(A)"),'Local Matrix:'
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
        
              if (n < 17) then
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
        
              if (n < 17) then
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
                call l2g(myj,mycol,pcol,nb,j)
                do myi=1,myXrows
                    call l2g(myi,myrow,prow,nb,i)
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

              if (me == 0) then
                print *,'Determinant = ', globdet
              endif

        contains
        
        ! convert global index to local index in block-cyclic distribution
        
           subroutine g2l(i,np,nb,p,il)
        
           implicit none
           integer, intent(in) :: i    ! global array index, input
           !integer, intent(in) :: n    ! global array dimension, input
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
        
           subroutine l2g(il,p,np,nb,i)
        
           implicit none
           integer :: il   ! local array index, input
           integer :: p    ! processor array index, input
           !integer :: n    ! global array dimension, input
           integer :: np   ! processor array dimension, input
           integer :: nb   ! block size, input
           integer :: i    ! global array index, output
           integer :: ilm1
        
           ilm1 = il-1
           i    = (((ilm1/nb) * np) + p)*nb + mod(ilm1,nb) + 1
        
           return
           end subroutine l2g
        
        end program gemv1
