      module readMatrix

      ! Used to make reading the bin file easy
      !USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE
 
      implicit none
       
      ! Define interface
      interface readMatrixIface
           module procedure readM
      end interface readMatrixIface

      contains
      subroutine readM(fname, n, M)
            implicit none
            character(*), intent(in)    :: fname
            integer, intent(in)           :: n
            double precision, intent(out) :: M(n,n)
      open(99,file=fname,form='unformatted',access='stream')
      read(99) M

      end subroutine readM

      end module readMatrix
