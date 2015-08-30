        module crandom_mod

        interface crandom_number
        module procedure crandom_number_d, crandom_number_z
        end interface
        
        contains

        subroutine crandom_number_z(n, A)
        implicit none
        integer n
        complex*16 A(n)

        integer i,iend,isize
        integer nb
        parameter(nb=1024)
        complex*16 zi
        real*8 x(nb),y(nb)

        zi = dcmplx(0.0,1.0)

        do i=1,n,nb
                iend = min(n,i+nb-1)
                isize = iend - i + 1
                call random_number(x)
                call random_number(y)
                A(i:iend) = x(1:isize) + y(1:isize)*zi
        enddo

        return
        end subroutine crandom_number_z


        subroutine crandom_number_d(n, A)
        implicit none
        integer n
        real*8 A(n)

        call random_number(A)

        return
        end subroutine crandom_number_d

        end module crandom_mod
