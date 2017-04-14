        program gemv1
              use mpi
              use readMatrix
              implicit none
        
              integer :: n, nb    ! problem size and block size
              integer :: myArows, myAcols   ! size of local subset of global matrix
              integer :: myXrows, myXcols   ! size of local subset of global vector
              integer :: i,j, myi, myj
              double precision, dimension(:,:), allocatable :: myA,myM
              integer, dimension(:), allocatable :: myX
              integer :: ierr
              character (len=10) :: arg
              character (len=10) :: fout
              character (len=50) :: fname
        
              integer, external :: numroc   ! blacs routine
              integer :: me, procs, icontxt, prow, pcol, myrow, mycol  ! blacs data
              integer :: info    ! scalapack return value
              integer, dimension(9)  :: ides_a, ides_x ! matrix descriptors
              integer, dimension(2) :: dims
              double precision :: det, globdet ! determinant and global
              integer :: starttime, laptime, cr ! for timing

              integer :: r,nr,c,nc
              integer :: sendr, sendc, recvr, recvc !for scatter calcs
              
        ! get problem size from input
        if (command_argument_count() > 0) then
              call get_command_argument(1, arg)
              READ (arg(:),'(i10)') n
        else
              print *, 'Need matrix size as parameter'
              call EXIT(1)
        endif

        ! start first timer
        call system_clock(count_rate=cr)
        call system_clock(starttime)
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
              nb = 8 !n/prow

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
        
              myXrows = numroc(n, nb, myrow, 0, prow) + nb
              myXcols = 1
        
        if (me == 0) write (*,'(I3A1I4A1I4A1)',advance="no"),
     & procs,',',n,',',nb,','

        ! Initialize local arrays
        
              allocate(myA(myArows,myAcols))
              allocate(myX(myXrows))

        ! Populate the matrix at root
        if (me == 0) then
            allocate(myM(n,n))
            if (command_argument_count() > 1) then
              ! use the filename to read matrix in
              call get_command_argument(2, fname)
              !READ (arg(:),'(a50)') fname
              !print *, 'Reading file ', fname
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

        if(me == 0) deallocate(myM)

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


              myX = 0

        
        ! Prepare array descriptors for ScaLAPACK 
            call descinit( ides_a, n, n, nb, nb, 0, 0, icontxt, 
     & myArows, info )
            call descinit( ides_x, n, 1, nb, nb, 0, 0, icontxt,
     & myXrows, info )
        
            if (me == 0) call logtime(laptime,cr,.false.)
        ! Call ScaLAPACK library routine
        !  SUBROUTINE PDGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
        call pdgetrf(n,n,myA,1,1,ides_a,myX,info)

              if (me == 0) then
                if (info /= 0) then
                     print *, 'Error -- info = ', info
                endif
              endif
        
              !print matrix out
            if (n < 17) then
              write (7,"(A13)") 'Decomposition'
              do myi=1,myArows
                do myj=1,myAcols 
                      write (7,"(f8.3)",advance="no"),myA(myi,myj)
                enddo
                write (7,"(A1)") ' '
              enddo

              write (7, "(A6)") 'Pivots'
              do myi=1,myXrows
                write (7, "(i3)",advance="no"),myX(myi)
              enddo
              write (7,"(A1)") ' '
            endif

            if (me == 0) call logtime(laptime,cr,.false.)
        
        ! Calculate determinant based on diagonal

            det = 0.d+0
            do myj=1,myAcols
                call l2g(myj,mycol,pcol,nb,j)
                do myi=1,myArows
                    call l2g(myi,myrow,prow,nb,i)
                    if (j == i) det = det + log10(abs(myA(myi,myj)))
                enddo
            enddo

            !print *, 'Local determinant ', det, 'proc', me
            if (me == 0) call logtime(laptime,cr,.false.)
        
        
              !    MPI_REDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, ROOT, COMM, IERROR)
              call MPI_Reduce(det, globdet, 1, MPI_DOUBLE, MPI_SUM, 0,
     &  MPI_COMM_WORLD, ierr)
        
            if (me == 0) call logtime(laptime,cr,.false.)

        ! Deallocate X
        
              deallocate(myA,myX)

        ! End blacs for processors that are used
        
              call blacs_gridexit(icontxt)
              call blacs_exit(0)
        
            if (me == 0) call logtime(starttime,cr,.true.)

              if (me == 0) then
                print *,'Log(|Det|) = ', globdet
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

        subroutine logtime(laptime,cr,eol)

                implicit none
                integer, intent(inout) :: laptime
                integer, intent(in)    :: cr
                logical, optional      :: eol
                integer :: stoptime

                call system_clock(stoptime)

                if(eol) then
                        write (*,"(e11.5)"),
     &                      real(stoptime-laptime)/cr
                else
                        write(*,"(e11.5a1)",advance="no"),
     &                      real(stoptime-laptime)/cr,','
                        laptime = stoptime
                endif

                return
                end subroutine logtime

        
        end program gemv1
