




module DRP_D0s
    use boundaries_types
    implicit none

contains

    subroutine D0s_DRP5_3Dx(f, ff, n1,n2,n2e,n3,n3e, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: ff

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.234102595369109d0     /2.d0
        A(2) =  -0.3184044667712d0     /2.d0
        A(3) =  0.11029870162898d0      /2.d0
        A(4) =  -0.030619166038291d0    /2.d0
        A(5) =  0.0046223358114003d0    /2.d0

        do k=1,n3e
            do j=1,n2e
                do i=5,n1-6
                    ff(i, j, k)= A(5)*(f(i+5, j, k) + f(i-4, j, k))     &
                    + A(4)*(f(i+4, j, k) + f(i-3, j, k))     &
                    + A(3)*(f(i+3, j, k) + f(i-2, j, k))     &
                    + A(2)*(f(i+2, j, k) + f(i-1, j, k))     &
                    + A(1)*(f(i+1, j, k) + f(i, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do k=1,n3e
                do j=1,n2e

                    do i=1, 4
                        ff(i, j, k)=   A(5)*(f(i+5, j, k) + f(mod(n1+i-6,n1-1)+1, j, k))  &
                        + A(4)*(f(i+4, j, k) + f(mod(n1+i-5,n1-1)+1, j, k))  &
                        + A(3)*(f(i+3, j, k) + f(mod(n1+i-4,n1-1)+1, j, k))  &
                        + A(2)*(f(i+2, j, k) + f(mod(n1+i-3,n1-1)+1, j, k))  &
                        + A(1)*(f(i+1, j, k) + f(i, j, k))
                    enddo

                    do i=n1-5,n1-1
                        ff(i, j, k)= A(5)*(f(mod(i+4,n1-1)+1, j, k) + f(i-4, j, k))             &
                        + A(4)*(f(mod(i+3,n1-1)+1, j, k) + f(i-3, j, k))             &
                        + A(3)*(f(mod(i+2,n1-1)+1, j, k) + f(i-2, j, k))             &
                        + A(2)*(f(mod(i+1,n1-1)+1, j, k) + f(i-1, j, k))             &
                        + A(1)*(f(mod(i,n1-1)+1, j, k) + f(i, j, k))
                    enddo
                enddo

            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1,n3e
                do j=1,n2e
                    ff(1,j,k)=0.5d0*(f(2,j,k)+f(1,j,k))
                    ff(2,j,k)  =   9.0d0/16.d0  * (f(3,j,k)+f(2,j,k))        &
                    -   0.0625d0  * (f(4,j,k)+f(1,j,k))
                    ff(3,j,k)  =   0.59395104312381d0  *   (f(4,j,k)+f(3,j,k))       &
                    -   0.10967656468571d0  *   (f(5,j,k)+f(2,j,k))             &
                    +   0.015725521561903d0 *   (f(6,j,k) + f(1,j,k))
                    ff(4,j,k)  =   0.59395104312381d0  *   (f(5,j,k)+f(4,j,k))       &
                    -   0.10967656468571d0  *   (f(6,j,k)+f(3,j,k))             &
                    +   0.015725521561903d0 *   (f(7,j,k) + f(2,j,k))
                    ff(5,j,k)  =   A(1)  *   (f(6,j,k)+f(5,j,k))       &
                    +   A(2)  *   (f(7,j,k)+f(4,j,k))     &
                    +   A(3) *   (f(8,j,k) + f(3,j,k))   &
                    +   A(4) *   (f(9,j,k) + f(2,j,k))   &
                    +   A(5)*   (f(10,j,k) + f(1,j,k))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    ff(n1-1,j,k)=0.5d0*(f(n1,j,k)+f(n1-1,j,k))
                    ff(n1-2,j,k)  =   9.0d0/16.d0 * (f(n1-1,j,k)+f(n1-2,j,k))     &
                    -   0.0625d0    * (f(n1,j,k)+f(n1-3,j,k))
                    ff(n1-3,j,k)  =   0.59395104312381d0  *   (f(n1-2,j,k)+f(n1-3,j,k))       &
                    -   0.10967656468571d0  *   (f(n1-1,j,k)+f(n1-4,j,k))             &
                    +   0.015725521561903d0 *   (f(n1,j,k) + f(n1-5,j,k))
                    ff(n1-4,j,k)  =   0.59395104312381d0  *   (f(n1-3,j,k)+f(n1-4,j,k))       &
                    -   0.10967656468571d0  *   (f(n1-2,j,k)+f(n1-5,j,k))             &
                    +   0.015725521561903d0 *   (f(n1-1,j,k) + f(n1-6,j,k))
                    ff(n1-5,j,k)  =   A(1)  *   (f(n1-4,j,k)+f(n1-5,j,k))       &
                    +   A(2)  *   (f(n1-3,j,k)+f(n1-6,j,k))             &
                    +   A(3) *   (f(n1-2,j,k) + f(n1-7,j,k))           &
                    +   A(4) *   (f(n1-1,j,k) + f(n1-8,j,k))           &
                    +   A(5) *   (f(n1,j,k) + f(n1-9,j,k))
                enddo
            enddo


        endif

        return

    end subroutine D0s_DRP5_3Dx



    subroutine D0s_DRP5_3Dy(f, ff, n1,n1e,n2,n3,n3e, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: ff

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.234102595369109d0     /2.d0
        A(2) =  -0.3184044667712d0     /2.d0
        A(3) =  0.11029870162898d0      /2.d0
        A(4) =  -0.030619166038291d0    /2.d0
        A(5) =  0.0046223358114003d0    /2.d0

        do k=1, n3e
            do i=1, n1e
                do j=5,n2-6
                    ff(i, j, k)= A(5)*(f(i, j+5, k) + f(i, j-4, k))     &
                    + A(4)*(f(i, j+4, k) + f(i, j-3, k))     &
                    + A(3)*(f(i, j+3, k) + f(i, j-2, k))     &
                    + A(2)*(f(i, j+2, k) + f(i, j-1, k))     &
                    + A(1)*(f(i, j+1, k) + f(i, j, k))
                enddo
            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=1, 4
                        ff(i, j, k)=   A(5)*(f(i, j+5, k) + f(i, mod(n2+j-6,n2-1)+1, k))  &
                        + A(4)*(f(i, j+4, k) + f(i, mod(n2+j-5,n2-1)+1, k))  &
                        + A(3)*(f(i, j+3, k) + f(i, mod(n2+j-4,n2-1)+1, k))  &
                        + A(2)*(f(i, j+2, k) + f(i, mod(n2+j-3,n2-1)+1, k))  &
                        + A(1)*(f(i, j+1, k) + f(i, j, k))
                    enddo

                    do j=n2-5,n2-1
                        ff(i, j, k)= A(5)*(f(i, mod(j+4,n2-1)+1, k) + f(i, j-4, k))             &
                        + A(4)*(f(i, mod(j+3,n2-1)+1, k) + f(i, j-3, k))             &
                        + A(3)*(f(i, mod(j+2,n2-1)+1, k) + f(i, j-2, k))             &
                        + A(2)*(f(i, mod(j+1,n2-1)+1, k) + f(i, j-1, k))             &
                        + A(1)*(f(i, mod(j,n2-1)+1, k) + f(i, j, k))
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e
                    ff(i,1,k)=0.5d0*(f(i,2,k)+f(i,1,k))
                    ff(i,2,k)  =   9.0d0/16.d0  * (f(i,3,k)+f(i,2,k))        &
                    -   0.0625d0  * (f(i,4,k)+f(i,1,k))
                    ff(i,3,k)  =   0.59395104312381d0  *   (f(i,4,k)+f(i,3,k))       &
                    -   0.10967656468571d0  *   (f(i,5,k)+f(i,2,k))             &
                    +   0.015725521561903d0 *   (f(i,6,k) + f(i,1,k))
                    ff(i,4,k)  =   0.59395104312381d0  *   (f(i,5,k)+f(i,4,k))       &
                    -   0.10967656468571d0  *   (f(i,6,k)+f(i,3,k))             &
                    +   0.015725521561903d0 *   (f(i,7,k) + f(i,2,k))
                    ff(i,5,k)  =   A(1)  *   (f(i,6,k)+f(i,5,k))       &
                    +   A(2)  *   (f(i,7,k)+f(i,4,k))     &
                    +   A(3) *   (f(i,8,k) + f(i,3,k))   &
                    +   A(4) *   (f(i,9,k) + f(i,2,k))   &
                    +   A(5) *   (f(i,10,k) + f(i,1,k))

                enddo
            enddo

            do k=1, n3e
                do i=1, n1e
                    ff(i,n2-1,k)=0.5d0*(f(i,n2,k)+f(i,n2-1,k))
                    ff(i,n2-2,k)  =   9.0d0/16.d0 * (f(i,n2-1,k)+f(i,n2-2,k))     &
                    -   0.0625d0    * (f(i,n2,k)+f(i,n2-3,k))
                    ff(i,n2-3,k)  =   0.59395104312381d0  *   (f(i,n2-2,k)+f(i,n2-3,k))       &
                    -   0.10967656468571d0  *   (f(i,n2-1,k)+f(i,n2-4,k))             &
                    +   0.015725521561903d0 *   (f(i,n2,k) + f(i,n2-5,k))
                    ff(i,n2-4,k)  =   0.59395104312381d0  *   (f(i,n2-3,k)+f(i,n2-4,k))       &
                    -   0.10967656468571d0  *   (f(i,n2-2,k)+f(i,n2-5,k))             &
                    +   0.015725521561903d0 *   (f(i,n2-1,k) + f(i,n2-6,k))
                    ff(i,n2-5,k)  =   A(1)  *   (f(i,n2-4,k)+f(i,n2-5,k))       &
                    +   A(2)  *   (f(i,n2-3,k)+f(i,n2-6,k))             &
                    +   A(3) *   (f(i,n2-2,k) + f(i,n2-7,k))           &
                    +   A(4) *   (f(i,n2-1,k) + f(i,n2-8,k))           &
                    +   A(5)*   (f(i,n2,k) + f(i,n2-9,k))

                enddo
            enddo

        endif

        return

    end subroutine D0s_DRP5_3Dy

    subroutine D0s_DRP5_3Dz(f, ff, n1,n1e,n2,n2e,n3, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: ff

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.234102595369109d0     /2.d0
        A(2) =  -0.3184044667712d0     /2.d0
        A(3) =  0.11029870162898d0      /2.d0
        A(4) =  -0.030619166038291d0    /2.d0
        A(5) =  0.0046223358114003d0    /2.d0

        do j=1, n2e
            do i=1, n1e
                do k=5,n3-6
                    ff(i, j, k)= A(5)*(f(i, j, k+5) + f(i, j, k-4))     &
                    + A(4)*(f(i, j, k+4) + f(i, j, k-3))     &
                    + A(3)*(f(i, j, k+3) + f(i, j, k-2))     &
                    + A(2)*(f(i, j, k+2) + f(i, j, k-1))     &
                    + A(1)*(f(i, j, k+1) + f(i, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e

                    do k=1, 4
                        ff(i, j, k)=   A(5)*(f(i, j, k+5) + f(i, j, mod(n3+k-6,n3-1)+1))  &
                        + A(4)*(f(i, j, k+4) + f(i, j, mod(n3+k-5,n3-1)+1))  &
                        + A(3)*(f(i, j, k+3) + f(i, j, mod(n3+k-4,n3-1)+1))  &
                        + A(2)*(f(i, j, k+2) + f(i, j, mod(n3+k-3,n3-1)+1))  &
                        + A(1)*(f(i, j, k+1) + f(i, j, k))
                    enddo

                    do k=n3-5,n3-1
                        ff(i, j, k)= A(5)*(f(i, j, mod(k+4,n3-1)+1) + f(i, j, k-4))             &
                        + A(4)*(f(i, j, mod(k+3,n3-1)+1) + f(i, j, k-3))             &
                        + A(3)*(f(i, j, mod(k+2,n3-1)+1) + f(i, j, k-2))             &
                        + A(2)*(f(i, j, mod(k+1,n3-1)+1) + f(i, j, k-1))             &
                        + A(1)*(f(i, j, mod(k,n3-1)+1) + f(i, j, k))
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e
                    ff(i,j,1)=0.5d0*(f(i,j,2)+f(i,j,1))
                    ff(i,j,2)  =   9.0d0/16.d0  * (f(i,j,3)+f(i,j,2))        &
                    -   0.0625d0  * (f(i,j,4)+f(i,j,1))
                    ff(i,j,3)  =   0.59395104312381d0  *   (f(i,j,4)+f(i,j,3))       &
                    -   0.10967656468571d0  *   (f(i,j,5)+f(i,j,2))             &
                    +   0.015725521561903d0 *   (f(i,j,6) + f(i,j,1))
                    ff(i,j,4)  =   0.59395104312381d0  *   (f(i,j,5)+f(i,j,4))       &
                    -   0.10967656468571d0  *   (f(i,j,6)+f(i,j,3))             &
                    +   0.015725521561903d0 *   (f(i,j,7) + f(i,j,2))
                    ff(i,j,5)  =   A(1)  *   (f(i,j,6)+f(i,j,5))       &
                    +   A(2)  *   (f(i,j,7)+f(i,j,4))     &
                    +   A(3) *   (f(i,j,8) + f(i,j,3))   &
                    +   A(4) *   (f(i,j,9) + f(i,j,2))   &
                    +   A(5) *   (f(i,j,10) + f(i,j,1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    ff(i,j,n3-1)=0.5d0*(f(i,j,n3)+f(i,j,n3-1))
                    ff(i,j,n3-2)  =   9.0d0/16.d0 * (f(i,j,n3-1)+f(i,j,n3-2))     &
                    -   0.0625d0    * (f(i,j,n3)+f(i,j,n3-3))
                    ff(i,j,n3-3)  =   0.59395104312381d0  *   (f(i,j,n3-2)+f(i,j,n3-3))       &
                    -   0.10967656468571d0  *   (f(i,j,n3-1)+f(i,j,n3-4))             &
                    +   0.015725521561903d0 *   (f(i,j,n3) + f(i,j,n3-5))
                    ff(i,j,n3-4)  =   0.59395104312381d0  *   (f(i,j,n3-3)+f(i,j,n3-4))       &
                    -   0.10967656468571d0  *   (f(i,j,n3-2)+f(i,j,n3-5))             &
                    +   0.015725521561903d0 *   (f(i,j,n3-1) + f(i,j,n3-6))
                    ff(i,j,n3-5)  =   A(1)  *   (f(i,j,n3-4)+f(i,j,n3-5))       &
                    +   A(2)  *   (f(i,j,n3-3)+f(i,j,n3-6))             &
                    +   A(3) *   (f(i,j,n3-2) + f(i,j,n3-7))           &
                    +   A(4) *   (f(i,j,n3-1) + f(i,j,n3-8))           &
                    +   A(5) *   (f(i,j,n3) + f(i,j,n3-9))
                enddo
            enddo

        endif

        return

    end subroutine D0s_DRP5_3Dz

end module DRP_D0s

module DRP_D0s_sh
    use boundaries_types
    implicit none

contains

    subroutine D0ssh_DRP5_3Dx(f, ff, n1,n2,n2e,n3,n3e, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: ff

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.234102595369109d0     /2.d0
        A(2) =  -0.3184044667712d0     /2.d0
        A(3) =  0.11029870162898d0      /2.d0
        A(4) =  -0.030619166038291d0    /2.d0
        A(5) =  0.0046223358114003d0    /2.d0

        do k=1,n3e
            do j=1,n2e
                do i=6,n1-5
                    ff(i, j, k)= A(5)*(f(i+4, j, k) + f(i-5, j, k))     &
                    + A(4)*(f(i+3, j, k) + f(i-4, j, k))     &
                    + A(3)*(f(i+2, j, k) + f(i-3, j, k))     &
                    + A(2)*(f(i+1, j, k) + f(i-2, j, k))     &
                    + A(1)*(f(i, j, k) + f(i-1, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do k=1,n3e
                do j=1,n2e

                    do i=1, 5
                        ff(i, j, k)=   A(5)*(f(i+4, j, k) + f(mod(n1+i-7,n1-1)+1, j, k))  &
                        + A(4)*(f(i+3, j, k) + f(mod(n1+i-6,n1-1)+1, j, k))  &
                        + A(3)*(f(i+2, j, k) + f(mod(n1+i-5,n1-1)+1, j, k))  &
                        + A(2)*(f(i+1, j, k) + f(mod(n1+i-4,n1-1)+1, j, k))  &
                        + A(1)*(f(i, j, k) + f(mod(n1+i-3,n1-1)+1, j, k))
                    enddo

                    do i=n1-4,n1-1
                        ff(i, j, k)= A(5)*(f(mod(i+3,n1-1)+1, j, k) + f(i-5, j, k))             &
                        + A(4)*(f(mod(i+2,n1-1)+1, j, k) + f(i-4, j, k))             &
                        + A(3)*(f(mod(i+1,n1-1)+1, j, k) + f(i-3, j, k))             &
                        + A(2)*(f(mod(i,n1-1)+1, j, k) + f(i-2, j, k))             &
                        + A(1)*(f(i, j, k) + f(i-1, j, k))
                    enddo
                enddo

            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1,n3e
                do j=1,n2e
                    ff(2,j,k)=0.5d0*(f(2,j,k)+f(1,j,k))
                    ff(3,j,k)  =   9.0d0/16.d0  * (f(3,j,k)+f(2,j,k))        &
                    -   0.0625d0  * (f(4,j,k)+f(1,j,k))
                    ff(4,j,k)  =   0.59395104312381d0  *   (f(4,j,k)+f(3,j,k))       &
                    -   0.10967656468571d0  *   (f(5,j,k)+f(2,j,k))             &
                    +   0.015725521561903d0 *   (f(6,j,k) + f(1,j,k))
                    ff(5,j,k)  =   0.59395104312381d0  *   (f(5,j,k)+f(4,j,k))       &
                    -   0.10967656468571d0  *   (f(6,j,k)+f(3,j,k))             &
                    +   0.015725521561903d0 *   (f(7,j,k) + f(2,j,k))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    ff(n1-1,j,k)=0.5d0*(f(n1-1,j,k)+f(n1-2,j,k))
                    ff(n1-2,j,k)  =   9.0d0/16.d0 * (f(n1-2,j,k)+f(n1-3,j,k))     &
                    -   0.0625d0    * (f(n1-1,j,k)+f(n1-4,j,k))
                    ff(n1-3,j,k)  =   0.59395104312381d0  *   (f(n1-3,j,k)+f(n1-4,j,k))       &
                    -   0.10967656468571d0  *   (f(n1-2,j,k)+f(n1-5,j,k))             &
                    +   0.015725521561903d0 *   (f(n1-1,j,k) + f(n1-6,j,k))
                    ff(n1-4,j,k)  =   0.59395104312381d0  *   (f(n1-4,j,k)+f(n1-5,j,k))       &
                    -   0.10967656468571d0  *   (f(n1-3,j,k)+f(n1-6,j,k))             &
                    +   0.015725521561903d0 *   (f(n1-2,j,k) + f(n1-7,j,k))
                enddo
            enddo


        endif

        return

    end subroutine D0ssh_DRP5_3Dx



    subroutine D0ssh_DRP5_3Dy(f, ff, n1,n1e,n2,n3,n3e, bc_type)

        implicit none

        integer, intent(in)     ::n1,n1e,n2,n3,n3e, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: ff

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.234102595369109d0     /2.d0
        A(2) =  -0.3184044667712d0     /2.d0
        A(3) =  0.11029870162898d0      /2.d0
        A(4) =  -0.030619166038291d0    /2.d0
        A(5) =  0.0046223358114003d0    /2.d0

        do k=1, n3e
            do i=1, n1e
                do j=5,n2-6
                    ff(i, j+1, k)= A(5)*(f(i, j+5, k) + f(i, j-4, k))     &
                    + A(4)*(f(i, j+4, k) + f(i, j-3, k))     &
                    + A(3)*(f(i, j+3, k) + f(i, j-2, k))     &
                    + A(2)*(f(i, j+2, k) + f(i, j-1, k))     &
                    + A(1)*(f(i, j+1, k) + f(i, j, k))
                enddo
            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=1, 5
                        ff(i, j, k)=   A(5)*(f(i, j+4, k) + f(i, mod(n2+j-7,n2-1)+1, k))  &
                        + A(4)*(f(i, j+3, k) + f(i, mod(n2+j-6,n2-1)+1, k))  &
                        + A(3)*(f(i, j+2, k) + f(i, mod(n2+j-5,n2-1)+1, k))  &
                        + A(2)*(f(i, j+1, k) + f(i, mod(n2+j-4,n2-1)+1, k))  &
                        + A(1)*(f(i, j, k) + f(i, mod(n2+j-3,n2-1)+1, k))
                    enddo

                    do j=n2-4,n2-1
                        ff(i, j, k)= A(5)*(f(i, mod(j+3,n2-1)+1, k) + f(i, j-5, k))             &
                        + A(4)*(f(i, mod(j+2,n2-1)+1, k) + f(i, j-4, k))             &
                        + A(3)*(f(i, mod(j+1,n2-1)+1, k) + f(i, j-3, k))             &
                        + A(2)*(f(i, mod(j,n2-1)+1, k) + f(i, j-2, k))             &
                        + A(1)*(f(i, j, k) + f(i, j-1, k))
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e
                    ff(i,2,k)=0.5d0*(f(i,2,k)+f(i,1,k))
                    ff(i,3,k)  =   9.0d0/16.d0  * (f(i,3,k)+f(i,2,k))        &
                    -   0.0625d0  * (f(i,4,k)+f(i,1,k))
                    ff(i,4,k)  =   0.59395104312381d0  *   (f(i,4,k)+f(i,3,k))       &
                    -   0.10967656468571d0  *   (f(i,5,k)+f(i,2,k))             &
                    +   0.015725521561903d0 *   (f(i,6,k) + f(i,1,k))
                    ff(i,5,k)  =   0.59395104312381d0  *   (f(i,5,k)+f(i,4,k))       &
                    -   0.10967656468571d0  *   (f(i,6,k)+f(i,3,k))             &
                    +   0.015725521561903d0 *   (f(i,7,k) + f(i,2,k))

                enddo
            enddo

            do k=1, n3e
                do i=1, n1e
                    ff(i,n2-1,k)=0.5d0*(f(i,n2-1,k)+f(i,n2-2,k))
                    ff(i,n2-2,k)  =   9.0d0/16.d0 * (f(i,n2-2,k)+f(i,n2-3,k))     &
                    -   0.0625d0    * (f(i,n2-1,k)+f(i,n2-4,k))
                    ff(i,n2-3,k)  =   0.59395104312381d0  *   (f(i,n2-3,k)+f(i,n2-4,k))       &
                    -   0.10967656468571d0  *   (f(i,n2-2,k)+f(i,n2-5,k))             &
                    +   0.015725521561903d0 *   (f(i,n2-1,k) + f(i,n2-6,k))
                    ff(i,n2-4,k)  =   0.59395104312381d0  *   (f(i,n2-4,k)+f(i,n2-5,k))       &
                    -   0.10967656468571d0  *   (f(i,n2-3,k)+f(i,n2-6,k))             &
                    +   0.015725521561903d0 *   (f(i,n2-2,k) + f(i,n2-7,k))

                enddo
            enddo

        endif

        return

    end subroutine D0ssh_DRP5_3Dy



    subroutine D0ssh_DRP5_MULT_3Dy(f, ff, n1,n1e,n2,n3,n3e, bc_type)

        implicit none

        integer, intent(in)     ::n1,n1e,n2,n3,n3e, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: ff

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.234102595369109d0     /2.d0
        A(2) =  -0.3184044667712d0     /2.d0
        A(3) =  0.11029870162898d0      /2.d0
        A(4) =  -0.030619166038291d0    /2.d0
        A(5) =  0.0046223358114003d0    /2.d0

        do k=1, n3e
            do i=1, n1e
                do j=5,n2-6
                    ff(i, j+1, k)= (A(5)*(f(i, j+5, k) + f(i, j-4, k))     &
                    + A(4)*(f(i, j+4, k) + f(i, j-3, k))     &
                    + A(3)*(f(i, j+3, k) + f(i, j-2, k))     &
                    + A(2)*(f(i, j+2, k) + f(i, j-1, k))     &
                    + A(1)*(f(i, j+1, k) + f(i, j, k)))*ff(i, j+1, k)
                enddo
            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=1, 5
                        ff(i, j, k)=   (A(5)*(f(i, j+4, k) + f(i, mod(n2+j-7,n2-1)+1, k))  &
                        + A(4)*(f(i, j+3, k) + f(i, mod(n2+j-6,n2-1)+1, k))  &
                        + A(3)*(f(i, j+2, k) + f(i, mod(n2+j-5,n2-1)+1, k))  &
                        + A(2)*(f(i, j+1, k) + f(i, mod(n2+j-4,n2-1)+1, k))  &
                        + A(1)*(f(i, j, k) + f(i, mod(n2+j-3,n2-1)+1, k)))*ff(i, j, k)
                    enddo

                    do j=n2-4,n2-1
                        ff(i, j, k)= (A(5)*(f(i, mod(j+3,n2-1)+1, k) + f(i, j-5, k))             &
                        + A(4)*(f(i, mod(j+2,n2-1)+1, k) + f(i, j-4, k))             &
                        + A(3)*(f(i, mod(j+1,n2-1)+1, k) + f(i, j-3, k))             &
                        + A(2)*(f(i, mod(j,n2-1)+1, k) + f(i, j-2, k))             &
                        + A(1)*(f(i, j, k) + f(i, j-1, k)))*ff(i, j, k)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e
                    ff(i,2,k)=(0.5d0*(f(i,2,k)+f(i,1,k)))*ff(i,2,k)
                    ff(i,3,k)  =   (9.0d0/16.d0  * (f(i,3,k)+f(i,2,k))        &
                    -   0.0625d0  * (f(i,4,k)+f(i,1,k)))*ff(i,3,k)
                    ff(i,4,k)  =   (0.59395104312381d0  *   (f(i,4,k)+f(i,3,k))       &
                    -   0.10967656468571d0  *   (f(i,5,k)+f(i,2,k))             &
                    +   0.015725521561903d0 *   (f(i,6,k) + f(i,1,k)))*ff(i,4,k)
                    ff(i,5,k)  =   (0.59395104312381d0  *   (f(i,5,k)+f(i,4,k))       &
                    -   0.10967656468571d0  *   (f(i,6,k)+f(i,3,k))             &
                    +   0.015725521561903d0 *   (f(i,7,k) + f(i,2,k)))*ff(i,5,k)

                enddo
            enddo

            do k=1, n3e
                do i=1, n1e
                    ff(i,n2-1,k)=(0.5d0*(f(i,n2-1,k)+f(i,n2-2,k)))*ff(i,n2-1,k)
                    ff(i,n2-2,k)  =   (9.0d0/16.d0 * (f(i,n2-2,k)+f(i,n2-3,k))     &
                    -   0.0625d0    * (f(i,n2-1,k)+f(i,n2-4,k)))*ff(i,n2-2,k)
                    ff(i,n2-3,k)  =   (0.59395104312381d0  *   (f(i,n2-3,k)+f(i,n2-4,k))       &
                    -   0.10967656468571d0  *   (f(i,n2-2,k)+f(i,n2-5,k))            &
                    +   0.015725521561903d0 *   (f(i,n2-1,k) + f(i,n2-6,k)))*ff(i,n2-3,k)
                    ff(i,n2-4,k)  =   (0.59395104312381d0  *   (f(i,n2-4,k)+f(i,n2-5,k))       &
                    -   0.10967656468571d0  *   (f(i,n2-3,k)+f(i,n2-6,k))            &
                    +   0.015725521561903d0 *   (f(i,n2-2,k) + f(i,n2-7,k)))*ff(i,n2-4,k)

                enddo
            enddo

        endif

        return

    end subroutine D0ssh_DRP5_MULT_3Dy

    subroutine D0ssh_DRP5_3Dz(f, ff, n1,n1e,n2,n2e,n3, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: ff

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.234102595369109d0     /2.d0
        A(2) =  -0.3184044667712d0     /2.d0
        A(3) =  0.11029870162898d0      /2.d0
        A(4) =  -0.030619166038291d0    /2.d0
        A(5) =  0.0046223358114003d0    /2.d0

        do j=1, n2e
            do i=1, n1e
                do k=5,n3-6
                    ff(i, j, k+1)= A(5)*(f(i, j, k+5) + f(i, j, k-4))     &
                    + A(4)*(f(i, j, k+4) + f(i, j, k-3))     &
                    + A(3)*(f(i, j, k+3) + f(i, j, k-2))     &
                    + A(2)*(f(i, j, k+2) + f(i, j, k-1))     &
                    + A(1)*(f(i, j, k+1) + f(i, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e

                    do k=1, 5
                        ff(i, j, k)=   A(5)*(f(i, j, k+4) + f(i, j, mod(n3+k-7,n3-1)+1))  &
                        + A(4)*(f(i, j, k+3) + f(i, j, mod(n3+k-6,n3-1)+1))  &
                        + A(3)*(f(i, j, k+2) + f(i, j, mod(n3+k-5,n3-1)+1))  &
                        + A(2)*(f(i, j, k+1) + f(i, j, mod(n3+k-4,n3-1)+1))  &
                        + A(1)*(f(i, j, k) + f(i, j, mod(n3+k-3,n3-1)+1))
                    enddo

                    do k=n3-4,n3-1
                        ff(i, j, k)= A(5)*(f(i, j, mod(k+3,n3-1)+1) + f(i, j, k-5))             &
                        + A(4)*(f(i, j, mod(k+2,n3-1)+1) + f(i, j, k-4))             &
                        + A(3)*(f(i, j, mod(k+1,n3-1)+1) + f(i, j, k-3))             &
                        + A(2)*(f(i, j, mod(k,n3-1)+1) + f(i, j, k-2))             &
                        + A(1)*(f(i, j, k) + f(i, j, k-1))
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e
                    ff(i,j,2)=0.5d0*(f(i,j,2)+f(i,j,1))
                    ff(i,j,3)  =   9.0d0/16.d0  * (f(i,j,3)+f(i,j,2))        &
                    -   0.0625d0  * (f(i,j,4)+f(i,j,1))
                    ff(i,j,4)  =   0.59395104312381d0  *   (f(i,j,4)+f(i,j,3))       &
                    -   0.10967656468571d0  *   (f(i,j,5)+f(i,j,2))             &
                    +   0.015725521561903d0 *   (f(i,j,6) + f(i,j,1))
                    ff(i,j,5)  =   0.59395104312381d0  *   (f(i,j,5)+f(i,j,4))       &
                    -   0.10967656468571d0  *   (f(i,j,6)+f(i,j,3))             &
                    +   0.015725521561903d0 *   (f(i,j,7) + f(i,j,2))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    ff(i,j,n3-1)=0.5d0*(f(i,j,n3-1)+f(i,j,n3-2))
                    ff(i,j,n3-2)  =   9.0d0/16.d0 * (f(i,j,n3-2)+f(i,j,n3-3))     &
                    -   0.0625d0    * (f(i,j,n3-1)+f(i,j,n3-4))
                    ff(i,j,n3-3)  =   0.59395104312381d0  *   (f(i,j,n3-3)+f(i,j,n3-4))       &
                    -   0.10967656468571d0  *   (f(i,j,n3-2)+f(i,j,n3-5))             &
                    +   0.015725521561903d0 *   (f(i,j,n3-1) + f(i,j,n3-6))
                    ff(i,j,n3-4)  =   0.59395104312381d0  *   (f(i,j,n3-4)+f(i,j,n3-5))       &
                    -   0.10967656468571d0  *   (f(i,j,n3-3)+f(i,j,n3-6))             &
                    +   0.015725521561903d0 *   (f(i,j,n3-2) + f(i,j,n3-7))
                enddo
            enddo

        endif

        return

    end subroutine D0ssh_DRP5_3Dz

    subroutine D0ssh_DRP5_MULT_3Dz(f, ff, n1,n1e,n2,n2e,n3, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: ff

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.234102595369109d0     /2.d0
        A(2) =  -0.3184044667712d0     /2.d0
        A(3) =  0.11029870162898d0      /2.d0
        A(4) =  -0.030619166038291d0    /2.d0
        A(5) =  0.0046223358114003d0    /2.d0

        do j=1, n2e
            do i=1, n1e
                do k=5,n3-6
                    ff(i, j, k+1)= (A(5)*(f(i, j, k+5) + f(i, j, k-4))     &
                    + A(4)*(f(i, j, k+4) + f(i, j, k-3))     &
                    + A(3)*(f(i, j, k+3) + f(i, j, k-2))     &
                    + A(2)*(f(i, j, k+2) + f(i, j, k-1))     &
                    + A(1)*(f(i, j, k+1) + f(i, j, k)))*ff(i, j, k+1)
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e

                    do k=1, 5
                        ff(i, j, k)=   (A(5)*(f(i, j, k+4) + f(i, j, mod(n3+k-7,n3-1)+1))  &
                        + A(4)*(f(i, j, k+3) + f(i, j, mod(n3+k-6,n3-1)+1))  &
                        + A(3)*(f(i, j, k+2) + f(i, j, mod(n3+k-5,n3-1)+1))  &
                        + A(2)*(f(i, j, k+1) + f(i, j, mod(n3+k-4,n3-1)+1))  &
                        + A(1)*(f(i, j, k) + f(i, j, mod(n3+k-3,n3-1)+1)))*ff(i, j, k)
                    enddo

                    do k=n3-4,n3-1
                        ff(i, j, k)= (A(5)*(f(i, j, mod(k+3,n3-1)+1) + f(i, j, k-5))             &
                        + A(4)*(f(i, j, mod(k+2,n3-1)+1) + f(i, j, k-4))             &
                        + A(3)*(f(i, j, mod(k+1,n3-1)+1) + f(i, j, k-3))             &
                        + A(2)*(f(i, j, mod(k,n3-1)+1) + f(i, j, k-2))             &
                        + A(1)*(f(i, j, k) + f(i, j, k-1)))*ff(i, j, k)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e
                    ff(i,j,2)=(0.5d0*(f(i,j,2)+f(i,j,1)))*ff(i,j,2)
                    ff(i,j,3)  =   (9.0d0/16.d0  * (f(i,j,3)+f(i,j,2))        &
                    -   0.0625d0  * (f(i,j,4)+f(i,j,1)))*ff(i,j,3)
                    ff(i,j,4)  =   (0.59395104312381d0  *   (f(i,j,4)+f(i,j,3))       &
                    -   0.10967656468571d0  *   (f(i,j,5)+f(i,j,2))             &
                    +   0.015725521561903d0 *   (f(i,j,6) + f(i,j,1)))*ff(i,j,4)
                    ff(i,j,5)  =   (0.59395104312381d0  *   (f(i,j,5)+f(i,j,4))       &
                    -   0.10967656468571d0  *   (f(i,j,6)+f(i,j,3))             &
                    +   0.015725521561903d0 *   (f(i,j,7) + f(i,j,2)))*ff(i,j,5)
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    ff(i,j,n3-1)=(0.5d0*(f(i,j,n3-1)+f(i,j,n3-2)))*ff(i,j,n3-1)
                    ff(i,j,n3-2)  =   (9.0d0/16.d0 * (f(i,j,n3-2)+f(i,j,n3-3))     &
                    -   0.0625d0    * (f(i,j,n3-1)+f(i,j,n3-4)))*ff(i,j,n3-2)
                    ff(i,j,n3-3)  =   (0.59395104312381d0  *   (f(i,j,n3-3)+f(i,j,n3-4))       &
                    -   0.10967656468571d0  *   (f(i,j,n3-2)+f(i,j,n3-5))             &
                    +   0.015725521561903d0 *   (f(i,j,n3-1) + f(i,j,n3-6)))*ff(i,j,n3-3)
                    ff(i,j,n3-4)  =   (0.59395104312381d0  *   (f(i,j,n3-4)+f(i,j,n3-5))       &
                    -   0.10967656468571d0  *   (f(i,j,n3-3)+f(i,j,n3-6))             &
                    +   0.015725521561903d0 *   (f(i,j,n3-2) + f(i,j,n3-7)))*ff(i,j,n3-4)
                enddo
            enddo

        endif

        return

    end subroutine D0ssh_DRP5_MULT_3Dz

end module DRP_D0s_sh
