      module readMatrix

      implicit none

      ! Define interface
      interface readMatrixIface
           module procedure readM, getRandM, getZeroM
      end interface readMatrixIface

      contains

      subroutine readM(fname, n, Mat)
            implicit none
            character(*), intent(in)    :: fname
            integer, intent(in)           :: n
            double precision, intent(out) :: Mat(n,n)
      open(99,file=fname,form='unformatted',access='stream')
      read(99) Mat
      end subroutine readM

      subroutine getRandM(n,m,Mat)
              implicit none
              integer, intent(in) :: n,m
              double precision, intent(out) :: Mat(n,m)
              integer, allocatable :: seed(:)
              integer :: randsize

              call random_seed(size=randsize)
              allocate(seed(randsize))
         OPEN(89,FILE='/dev/urandom',ACCESS='stream',FORM='UNFORMATTED')
              READ(89) seed
              CLOSE(89)
              call random_seed(put=seed)
              call random_number(Mat)
              Mat = Mat - 0.5d+0
      end subroutine getRandM

      subroutine getZeroM(m,Mat)
        implicit none
        integer, intent(in) :: n,m
        integer :: k,i,j,l,p
        double precision, intent(out) :: Mat(n,m)
        integer, allocatable :: seed(:)
        integer :: randsize
        real :: r

        call random_seed(size=randsize)
        allocate(seed(randsize))
        OPEN(89,FILE='/dev/urandom',ACCESS='stream',FORM='UNFORMATTED')
        READ(89) seed
        CLOSE(89)
        call random_seed(put=seed)
        call random_number(Mat)
        Mat = Mat - 0.5d+0

        ! Get a random number
        call random_number(r)
        i = floor( r*n )

        ! Get another
        call random_number(r)
        j = floor( r*n )

        ! Copy a row
        do k = 1,5
          Mat(i,k) = Mat(j,k)
        end do

        ! Output matrix to debug
        !do l = 1,5
        !    Print *, (Mat(l,p), p=1,5)
        !end do
      end subroutine getZeroM

      end module readMatrix
