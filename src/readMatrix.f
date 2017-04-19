      module readMatrix

      implicit none

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

      subroutine getZeroM(n,m,Mat)
        implicit none
        integer, intent(in) :: n,m
        integer :: i,j
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

        Print *, Mat
        ! Copy a row
        Mat(i,1:m) = Mat(j,1:m)

        print*,"row",j,"copied to row",i

        ! Output matrix to debug
        Print *, Mat
      end subroutine getZeroM


      end module readMatrix
