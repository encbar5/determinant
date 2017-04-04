!comments here
        module example
            implicit none
            real (kind=8), parameter :: z = 3.14d+0
            interface dostuff
                module procedure printerup
            end interface

        contains
            subroutine printerup
                print *, z
            end subroutine printerup

        end module example
