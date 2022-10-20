module DRP_D2c
    use boundaries_types
    implicit none

contains

    subroutine D2c_DRP6_3Dx(f, d2f, n1,n2,n2e,n3,n3e, h, shifted, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: d2f

        integer :: i,j,k
        real(kind=8) A(0:6)

        A(1) =  1.812906831737684d0         / h**2
        A(2) =  -0.33473276510983d0         / h**2
        A(3) =  0.087046253270279d0         / h**2
        A(4) =  -0.021406336650382d0        / h**2
        A(5) =  0.0039703382061228d0        / h**2
        A(6) = -3.9303109660650135d-4       / h**2
        A(0) = -2.d0*(A(1) + A(2) + A(3) + A(4) + A(5) + A(6))

        do k=1,n3e
            do j=1,n2e
                do i=7,n1-7
                    d2f(i, j, k)= A(6)*(f(i+6, j, k) + f(i-6, j, k))     &
                    + A(5)*(f(i+5, j, k) + f(i-5, j, k))     &
                    + A(4)*(f(i+4, j, k) + f(i-4, j, k))     &
                    + A(3)*(f(i+3, j, k) + f(i-3, j, k))     &
                    + A(2)*(f(i+2, j, k) + f(i-2, j, k))     &
                    + A(1)*(f(i+1, j, k) + f(i-1, j, k))     &
                    + A(0)*f(i, j, k)
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do k=1,n3e
                do j=1,n2e

                    do i=1, 6
                        d2f(i, j, k)=   A(6)*(f(i+6, j, k) + f(mod(n1+i-8,n1-1)+1, j, k))  &
                        + A(5)*(f(i+5, j, k) + f(mod(n1+i-7,n1-1)+1, j, k))  &
                        + A(4)*(f(i+4, j, k) + f(mod(n1+i-6,n1-1)+1, j, k))  &
                        + A(3)*(f(i+3, j, k) + f(mod(n1+i-5,n1-1)+1, j, k))  &
                        + A(2)*(f(i+2, j, k) + f(mod(n1+i-4,n1-1)+1, j, k))  &
                        + A(1)*(f(i+1, j, k) + f(mod(n1+i-3,n1-1)+1, j, k))  &
                        + A(0)*f(i, j, k)
                    enddo

                    do i=n1-6, n1-1
                        d2f(i, j, k)= A(6)*(f(mod(i+5,n1-1)+1, j, k) + f(i-6, j, k))             &
                        + A(5)*(f(mod(i+4,n1-1)+1, j, k) + f(i-5, j, k))             &
                        + A(4)*(f(mod(i+3,n1-1)+1, j, k) + f(i-4, j, k))             &
                        + A(3)*(f(mod(i+2,n1-1)+1, j, k) + f(i-3, j, k))             &
                        + A(2)*(f(mod(i+1,n1-1)+1, j, k) + f(i-2, j, k))             &
                        + A(1)*(f(mod(i,n1-1)+1, j, k) + f(i-1, j, k))               &
                        + A(0)*f(i, j, k)
                    enddo
                enddo

            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1,n3e
                do j=1,n2e
                    d2f(2,j,k)  =   ( -2.0d0*f(2,j,k) + (f(3,j,k)+f(1,j,k)) ) / h**2
                    d2f(3,j,k)  =   -5.0d0  * f(3,j,k)  / (2.d0*h**2)   &
                    +   4.d0  * (f(4,j,k)+f(2,j,k)) / (3.d0*h**2)   &
                    -   1.d0 * (f(5,j,k) + f(1,j,k))    / (12.d0*h**2)
                    d2f(4,j,k)  =   - 2.814728882213931d0 * f(4,j,k)    / h**2      &
                    +   1.569379994993781d0*(f(5,j,k)+f(3,j,k))     / h**2      &
                    -   0.177751998d0*(f(6,j,k)+f(2,j,k))           / h**2      &
                    +   0.0157364441d0* (f(7,j,k) + f(1,j,k))           / h**2
                    d2f(5,j,k)  =   - 2.814728882213931d0 * f(5,j,k)    / h**2      &
                    +   1.569379994993781d0*(f(6,j,k)+f(4,j,k))     / h**2      &
                    -   0.177751998d0*(f(7,j,k)+f(3,j,k))           / h**2      &
                    +   0.0157364441d0* (f(8,j,k) + f(2,j,k))           / h**2
                    d2f(6,j,k)  =   (                                           &
                    -   3.01824892312857d0      * f(6,j,k)                      &
                    +   1.744086462579369d0     * (f(7,j,k)+f(5,j,k))           &
                    -   0.2851293540782d0      * (f(8,j,k)+f(4,j,k))           &
                    +   0.059207636829005d0      * (f(9,j,k) + f(3,j,k))         &
                    -   0.0099521462688359d0     * (f(10,j,k) + f(2,j,k))        &
                    +   9.1186250295123312d-4    * (f(11,j,k) + f(1,j,k))        &
                    ) /h**2

                enddo
            enddo



            if (shifted) then

                do k=1,n3e
                    do j=1,n2e
                        d2f(n1-2,j,k)  =   ( -2.0d0*f(n1-2,j,k) + (f(n1-1,j,k)+f(n1-3,j,k)) )   / h**2
                        d2f(n1-3,j,k)  =   -5.0d0  * f(n1-3,j,k)                                / (2.d0*h**2)   &
                        +   4.d0  * (f(n1-2,j,k)+f(n1-4,j,k))                                   / (3.d0*h**2)   &
                        -   1.d0 * (f(n1-1,j,k) + f(n1-5,j,k))                                  / (12.d0*h**2)
                        d2f(n1-4,j,k)  =   - 2.814728882213931d0 * f(n1-4,j,k)                  / h**2      &
                        +   1.569379994993781d0*(f(n1-3,j,k)+f(n1-5,j,k))                       / h**2      &
                        -   0.177751998d0*(f(n1-2,j,k)+f(n1-6,j,k))                             / h**2      &
                        +   0.0157364441d0* (f(n1-1,j,k) + f(n1-7,j,k))                         / h**2
                        d2f(n1-5,j,k)  =   - 2.814728882213931d0 * f(n1-5,j,k)                  / h**2      &
                        +   1.569379994993781d0*(f(n1-4,j,k)+f(n1-6,j,k))                       / h**2      &
                        -   0.177751998d0*(f(n1-3,j,k)+f(n1-7,j,k))                             / h**2      &
                        +   0.0157364441d0* (f(n1-2,j,k) + f(n1-8,j,k))                         / h**2
                        d2f(n1-6,j,k)  =   (                                                                &
                        -   3.01824892312857d0      * f(n1-6,j,k)                                           &
                        +   1.744086462579369d0     * (f(n1-5,j,k)+f(n1-7,j,k))                             &
                        -   0.2851293540782d0      * (f(n1-4,j,k)+f(n1-8,j,k))                             &
                        +   0.059207636829005d0      * (f(n1-3,j,k) + f(n1-9,j,k))                           &
                        -   0.0099521462688359d0     * (f(n1-2,j,k) + f(n1-10,j,k))                          &
                        +   9.1186250295123312d-4    * (f(n1-1,j,k) + f(n1-11,j,k))                          &
                        )                                                                       /h**2

                    enddo
                enddo

            else

                do k=1,n3e
                    do j=1,n2e
                        d2f(n1-1,j,k)  =   ( -2.0d0*f(n1-1,j,k) + (f(n1,j,k)+f(n1-2,j,k)) )     / h**2
                        d2f(n1-2,j,k)  =   -5.0d0  * f(n1-2,j,k)                                / (2.d0*h**2)   &
                        +   4.d0  * (f(n1-1,j,k)+f(n1-3,j,k))                                   / (3.d0*h**2)   &
                        -   1.d0 * (f(n1,j,k) + f(n1-4,j,k))                                    / (12.d0*h**2)
                        d2f(n1-3,j,k)  =   - 2.814728882213931d0 * f(n1-3,j,k)                  / h**2      &
                        +   1.569379994993781d0*(f(n1-2,j,k)+f(n1-4,j,k))                       / h**2      &
                        -   0.177751998d0*(f(n1-1,j,k)+f(n1-5,j,k))                             / h**2      &
                        +   0.0157364441d0* (f(n1,j,k) + f(n1-6,j,k))                           / h**2
                        d2f(n1-4,j,k)  =   - 2.814728882213931d0 * f(n1-4,j,k)                  / h**2      &
                        +   1.569379994993781d0*(f(n1-3,j,k)+f(n1-5,j,k))                       / h**2      &
                        -   0.177751998d0*(f(n1-2,j,k)+f(n1-6,j,k))                             / h**2      &
                        +   0.0157364441d0* (f(n1-1,j,k) + f(n1-7,j,k))                         / h**2
                        d2f(n1-5,j,k)  =   (                                                                &
                        -   3.01824892312857d0      * f(n1-5,j,k)                                           &
                        +   1.744086462579369d0     * (f(n1-4,j,k)+f(n1-6,j,k))                             &
                        -   0.2851293540782d0      * (f(n1-3,j,k)+f(n1-7,j,k))                             &
                        +   0.059207636829005d0      * (f(n1-2,j,k) + f(n1-8,j,k))                           &
                        -   0.0099521462688359d0     * (f(n1-1,j,k) + f(n1-9,j,k))                           &
                        +   9.1186250295123312d-4    * (f(n1,j,k) + f(n1-10,j,k))                            &
                        )                                                                       /h**2
                        d2f(n1-6,j,k)  =   (                                                                &
                        +   A(0)    * f(n1-6,j,k)                                           &
                        +   A(1)     * (f(n1-5,j,k)+f(n1-7,j,k))                             &
                        +   A(2)      * (f(n1-4,j,k)+f(n1-8,j,k))                             &
                        +   A(3)      * (f(n1-3,j,k) + f(n1-9,j,k))                           &
                        +   A(4)      * (f(n1-2,j,k) + f(n1-10,j,k))                          &
                        +   A(5)     * (f(n1-1,j,k) + f(n1-11,j,k))                          &
                        +   A(6)    * (f(n1,j,k) + f(n1-12,j,k))                            &
                        )

                    enddo
                enddo

            end if

        endif

        return

    end subroutine D2c_DRP6_3Dx

    subroutine D2c_DRP6_ACC_3Dx(f, d2f, n1,n2,n2e,n3,n3e, h, shifted, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: d2f

        integer :: i,j,k
        real(kind=8) A(0:6)

        A(1) =  1.812906831737684d0         / h**2
        A(2) =  -0.33473276510983d0         / h**2
        A(3) =  0.087046253270279d0         / h**2
        A(4) =  -0.021406336650382d0        / h**2
        A(5) =  0.0039703382061228d0        / h**2
        A(6) = -3.9303109660650135d-4       / h**2
        A(0) = -2.d0*(A(1) + A(2) + A(3) + A(4) + A(5) + A(6))

        do k=1,n3e
            do j=1,n2e
                do i=7,n1-7
                    d2f(i, j, k)= d2f(i, j, k)+(A(6)*(f(i+6, j, k) + f(i-6, j, k))     &
                    + A(5)*(f(i+5, j, k) + f(i-5, j, k))     &
                    + A(4)*(f(i+4, j, k) + f(i-4, j, k))     &
                    + A(3)*(f(i+3, j, k) + f(i-3, j, k))     &
                    + A(2)*(f(i+2, j, k) + f(i-2, j, k))     &
                    + A(1)*(f(i+1, j, k) + f(i-1, j, k))     &
                    + A(0)*f(i, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do k=1,n3e
                do j=1,n2e

                    do i=1, 6
                        d2f(i, j, k)=   d2f(i, j, k)+(A(6)*(f(i+6, j, k) + f(mod(n1+i-8,n1-1)+1, j, k))  &
                        + A(5)*(f(i+5, j, k) + f(mod(n1+i-7,n1-1)+1, j, k))  &
                        + A(4)*(f(i+4, j, k) + f(mod(n1+i-6,n1-1)+1, j, k))  &
                        + A(3)*(f(i+3, j, k) + f(mod(n1+i-5,n1-1)+1, j, k))  &
                        + A(2)*(f(i+2, j, k) + f(mod(n1+i-4,n1-1)+1, j, k))  &
                        + A(1)*(f(i+1, j, k) + f(mod(n1+i-3,n1-1)+1, j, k))  &
                        + A(0)*f(i, j, k))
                    enddo

                    do i=n1-6, n1-1
                        d2f(i, j, k)= d2f(i, j, k)+(A(6)*(f(mod(i+5,n1-1)+1, j, k) + f(i-6, j, k))             &
                        + A(5)*(f(mod(i+4,n1-1)+1, j, k) + f(i-5, j, k))             &
                        + A(4)*(f(mod(i+3,n1-1)+1, j, k) + f(i-4, j, k))             &
                        + A(3)*(f(mod(i+2,n1-1)+1, j, k) + f(i-3, j, k))             &
                        + A(2)*(f(mod(i+1,n1-1)+1, j, k) + f(i-2, j, k))             &
                        + A(1)*(f(mod(i,n1-1)+1, j, k) + f(i-1, j, k))               &
                        + A(0)*f(i, j, k))
                    enddo
                enddo

            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1,n3e
                do j=1,n2e
                    d2f(2,j,k)  =   d2f(2, j, k)+(( -2.0d0*f(2,j,k) + (f(3,j,k)+f(1,j,k)) ) / h**2)
                    d2f(3,j,k)  =   d2f(3, j, k)+(-5.0d0  * f(3,j,k)  / (2.d0*h**2)   &
                    +   4.d0  * (f(4,j,k)+f(2,j,k)) / (3.d0*h**2)   &
                    -   1.d0 * (f(5,j,k) + f(1,j,k))    / (12.d0*h**2))
                    d2f(4,j,k)  =   d2f(4, j, k)+(- 2.814728882213931d0 * f(4,j,k)    / h**2      &
                    +   1.569379994993781d0*(f(5,j,k)+f(3,j,k))     / h**2      &
                    -   0.177751998d0*(f(6,j,k)+f(2,j,k))           / h**2      &
                    +   0.0157364441d0* (f(7,j,k) + f(1,j,k))           / h**2)
                    d2f(5,j,k)  =   d2f(5, j, k)+(- 2.814728882213931d0 * f(5,j,k)    / h**2      &
                    +   1.569379994993781d0*(f(6,j,k)+f(4,j,k))     / h**2      &
                    -   0.177751998d0*(f(7,j,k)+f(3,j,k))           / h**2      &
                    +   0.0157364441d0* (f(8,j,k) + f(2,j,k))           / h**2)
                    d2f(6,j,k)  =   d2f(6, j, k)+((                                           &
                    -   3.01824892312857d0      * f(6,j,k)                      &
                    +   1.744086462579369d0     * (f(7,j,k)+f(5,j,k))           &
                    -   0.2851293540782d0      * (f(8,j,k)+f(4,j,k))           &
                    +   0.059207636829005d0      * (f(9,j,k) + f(3,j,k))         &
                    -   0.0099521462688359d0     * (f(10,j,k) + f(2,j,k))        &
                    +   9.1186250295123312d-4    * (f(11,j,k) + f(1,j,k))        &
                    ) /h**2)

                enddo
            enddo



            if (shifted) then

                do k=1,n3e
                    do j=1,n2e
                        d2f(n1-2,j,k)  =   d2f(n1-2,j,k)+( -2.0d0*f(n1-2,j,k) + (f(n1-1,j,k)+f(n1-3,j,k)) )   / h**2
                        d2f(n1-3,j,k)  =   d2f(n1-3,j,k)-5.0d0  * f(n1-3,j,k)                                / (2.d0*h**2)   &
                        +   4.d0  * (f(n1-2,j,k)+f(n1-4,j,k))                                   / (3.d0*h**2)   &
                        -   1.d0 * (f(n1-1,j,k) + f(n1-5,j,k))                                  / (12.d0*h**2)
                        d2f(n1-4,j,k)  =   d2f(n1-4,j,k)- 2.814728882213931d0 * f(n1-4,j,k)                  / h**2      &
                        +   1.569379994993781d0*(f(n1-3,j,k)+f(n1-5,j,k))                       / h**2      &
                        -   0.177751998d0*(f(n1-2,j,k)+f(n1-6,j,k))                             / h**2      &
                        +   0.0157364441d0* (f(n1-1,j,k) + f(n1-7,j,k))                         / h**2
                        d2f(n1-5,j,k)  =   d2f(n1-5,j,k)- 2.814728882213931d0 * f(n1-5,j,k)                  / h**2      &
                        +   1.569379994993781d0*(f(n1-4,j,k)+f(n1-6,j,k))                       / h**2      &
                        -   0.177751998d0*(f(n1-3,j,k)+f(n1-7,j,k))                             / h**2      &
                        +   0.0157364441d0* (f(n1-2,j,k) + f(n1-8,j,k))                         / h**2
                        d2f(n1-6,j,k)  =   d2f(n1-6,j,k)+(                                                                &
                        -   3.01824892312857d0      * f(n1-6,j,k)                                           &
                        +   1.744086462579369d0     * (f(n1-5,j,k)+f(n1-7,j,k))                             &
                        -   0.2851293540782d0      * (f(n1-4,j,k)+f(n1-8,j,k))                             &
                        +   0.059207636829005d0      * (f(n1-3,j,k) + f(n1-9,j,k))                           &
                        -   0.0099521462688359d0     * (f(n1-2,j,k) + f(n1-10,j,k))                          &
                        +   9.1186250295123312d-4    * (f(n1-1,j,k) + f(n1-11,j,k))                          &
                        )                                                                       /h**2

                    enddo
                enddo

            else

                do k=1,n3e
                    do j=1,n2e
                        d2f(n1-1,j,k)  =   d2f(n1-1,j,k)+( -2.0d0*f(n1-1,j,k) + (f(n1,j,k)+f(n1-2,j,k)) )     / h**2
                        d2f(n1-2,j,k)  =   d2f(n1-2,j,k)-5.0d0  * f(n1-2,j,k)                                / (2.d0*h**2)   &
                        +   4.d0  * (f(n1-1,j,k)+f(n1-3,j,k))                                   / (3.d0*h**2)   &
                        -   1.d0 * (f(n1,j,k) + f(n1-4,j,k))                                    / (12.d0*h**2)
                        d2f(n1-3,j,k)  =   d2f(n1-3,j,k)- 2.814728882213931d0 * f(n1-3,j,k)                  / h**2      &
                        +   1.569379994993781d0*(f(n1-2,j,k)+f(n1-4,j,k))                       / h**2      &
                        -   0.177751998d0*(f(n1-1,j,k)+f(n1-5,j,k))                             / h**2      &
                        +   0.0157364441d0* (f(n1,j,k) + f(n1-6,j,k))                           / h**2
                        d2f(n1-4,j,k)  =   d2f(n1-4,j,k)- 2.814728882213931d0 * f(n1-4,j,k)                  / h**2      &
                        +   1.569379994993781d0*(f(n1-3,j,k)+f(n1-5,j,k))                       / h**2      &
                        -   0.177751998d0*(f(n1-2,j,k)+f(n1-6,j,k))                             / h**2      &
                        +   0.0157364441d0* (f(n1-1,j,k) + f(n1-7,j,k))                         / h**2
                        d2f(n1-5,j,k)  =   d2f(n1-5,j,k)+(                                                                &
                        -   3.01824892312857d0      * f(n1-5,j,k)                                           &
                        +   1.744086462579369d0     * (f(n1-4,j,k)+f(n1-6,j,k))                             &
                        -   0.2851293540782d0      * (f(n1-3,j,k)+f(n1-7,j,k))                             &
                        +   0.059207636829005d0      * (f(n1-2,j,k) + f(n1-8,j,k))                           &
                        -   0.0099521462688359d0     * (f(n1-1,j,k) + f(n1-9,j,k))                           &
                        +   9.1186250295123312d-4    * (f(n1,j,k) + f(n1-10,j,k))                            &
                        )                                                                       /h**2
                        d2f(n1-6,j,k)  =   d2f(n1-6,j,k)+(                                                                &
                        +   A(0)    * f(n1-6,j,k)                                           &
                        +   A(1)     * (f(n1-5,j,k)+f(n1-7,j,k))                             &
                        +   A(2)      * (f(n1-4,j,k)+f(n1-8,j,k))                             &
                        +   A(3)      * (f(n1-3,j,k) + f(n1-9,j,k))                           &
                        +   A(4)      * (f(n1-2,j,k) + f(n1-10,j,k))                          &
                        +   A(5)     * (f(n1-1,j,k) + f(n1-11,j,k))                          &
                        +   A(6)    * (f(n1,j,k) + f(n1-12,j,k))                            &
                        )

                    enddo
                enddo

            end if

        endif

        return

    end subroutine D2c_DRP6_ACC_3Dx



    subroutine D2c_DRP6_3Dy(f, d2f, n1,n1e,n2,n3,n3e, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: d2f

        integer :: i,j,k
        real(kind=8) A(0:6)

        A(1) =  1.812906831737684d0         / h**2
        A(2) =  -0.33473276510983d0         / h**2
        A(3) =  0.087046253270279d0         / h**2
        A(4) =  -0.021406336650382d0        / h**2
        A(5) =  0.0039703382061228d0        / h**2
        A(6) = -3.9303109660650135d-4       / h**2
        A(0) = -2.d0*(A(1) + A(2) + A(3) + A(4) + A(5) + A(6))

        do k=1, n3e
            do i=1, n1e

                do j=7,n2-7
                    d2f(i,j,k)= A(6)*(f(i,j+6,k) + f(i,j-6,k))     &
                    + A(5)*(f(i,j+5,k) + f(i,j-5,k))     &
                    + A(4)*(f(i,j+4,k) + f(i,j-4,k))     &
                    + A(3)*(f(i,j+3,k) + f(i,j-3,k))     &
                    + A(2)*(f(i,j+2,k) + f(i,j-2,k))     &
                    + A(1)*(f(i,j+1,k) + f(i,j-1,k))     &
                    + A(0)*f(i,j,k)
                enddo

            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=1, 6
                        d2f(i,j,k)=   A(6)*(f(i,j+6,k) + f(i,mod(n2+j-8,n2-1)+1,k))  &
                        + A(5)*(f(i,j+5,k) + f(i,mod(n2+j-7,n2-1)+1,k))  &
                        + A(4)*(f(i,j+4,k) + f(i,mod(n2+j-6,n2-1)+1,k))  &
                        + A(3)*(f(i,j+3,k) + f(i,mod(n2+j-5,n2-1)+1,k))  &
                        + A(2)*(f(i,j+2,k) + f(i,mod(n2+j-4,n2-1)+1,k))  &
                        + A(1)*(f(i,j+1,k) + f(i,mod(n2+j-3,n2-1)+1,k)) &
                        + A(0)*f(i,j,k)
                    enddo

                    do j=n2-6, n2-1
                        d2f(i,j,k)= A(6)*(f(i,mod(j+5,n2-1)+1,k) + f(i,j-6,k))             &
                        + A(5)*(f(i,mod(j+4,n2-1)+1,k) + f(i,j-5,k))             &
                        + A(4)*(f(i,mod(j+3,n2-1)+1,k) + f(i,j-4,k))             &
                        + A(3)*(f(i,mod(j+2,n2-1)+1,k) + f(i,j-3,k))             &
                        + A(2)*(f(i,mod(j+1,n2-1)+1,k) + f(i,j-2,k))             &
                        + A(1)*(f(i,mod(j,n2-1)+1,k) + f(i,j-1,k))               &
                        + A(0)*f(i,j,k)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e
                    d2f(i,2,k)  =   ( -2.0d0*f(i,2,k) + (f(i,3,k)+f(i,1,k)) ) / h**2
                    d2f(i,3,k)  =   -5.0d0  * f(i,3,k)  / (2.d0*h**2)   &
                    +   4.d0  * (f(i,4,k)+f(i,2,k)) / (3.d0*h**2)   &
                    -   1.d0 * (f(i,5,k) + f(i,1,k))    / (12.d0*h**2)
                    d2f(i,4,k)  =   - 2.814728882213931d0 * f(i,4,k)    / h**2      &
                    +   1.569379994993781d0*(f(i,5,k)+f(i,3,k))     / h**2      &
                    -   0.177751998d0*(f(i,6,k)+f(i,2,k))           / h**2      &
                    +   0.0157364441d0* (f(i,7,k) + f(i,1,k))           / h**2
                    d2f(i,5,k)  =   - 2.814728882213931d0 * f(i,5,k)    / h**2      &
                    +   1.569379994993781d0*(f(i,6,k)+f(i,4,k))     / h**2      &
                    -   0.177751998d0*(f(i,7,k)+f(i,3,k))           / h**2      &
                    +   0.0157364441d0* (f(i,8,k) + f(i,2,k))           / h**2
                    d2f(i,6,k)  =   (                                           &
                    -   3.01824892312857d0      * f(i,6,k)                      &
                    +   1.744086462579369d0     * (f(i,7,k)+f(i,5,k))           &
                    -   0.2851293540782d0      * (f(i,8,k)+f(i,4,k))           &
                    +   0.059207636829005d0      * (f(i,9,k) + f(i,3,k))         &
                    -   0.0099521462688359d0     * (f(i,10,k) + f(i,2,k))        &
                    +   9.1186250295123312d-4    * (f(i,11,k) + f(i,1,k))        &
                    ) /h**2
                enddo
            enddo



            if (shifted) then

                do k=1, n3e
                    do i=1, n1e
                        d2f(i,n2-2,k)  =   ( -2.0d0*f(i,n2-2,k) + (f(i,n2-1,k)+f(i,n2-3,k)) ) / h**2
                        d2f(i,n2-3,k)  =   -5.0d0  * f(i,n2-3,k)  / (2.d0*h**2)   &
                        +   4.d0  * (f(i,n2-2,k)+f(i,n2-4,k)) / (3.d0*h**2)   &
                        -   1.d0 * (f(i,n2-1,k) + f(i,n2-5,k))    / (12.d0*h**2)
                        d2f(i,n2-4,k)  =   - 2.814728882213931d0 * f(i,n2-4,k)    / h**2      &
                        +   1.569379994993781d0*(f(i,n2-3,k)+f(i,n2-5,k))     / h**2      &
                        -   0.177751998d0*(f(i,n2-2,k)+f(i,n2-6,k))           / h**2      &
                        +   0.0157364441d0* (f(i,n2-1,k) + f(i,n2-7,k))           / h**2
                        d2f(i,n2-5,k)  =   - 2.814728882213931d0 * f(i,n2-5,k)    / h**2      &
                        +   1.569379994993781d0*(f(i,n2-4,k)+f(i,n2-6,k))     / h**2      &
                        -   0.177751998d0*(f(i,n2-3,k)+f(i,n2-7,k))           / h**2      &
                        +   0.0157364441d0* (f(i,n2-2,k) + f(i,n2-8,k))           / h**2
                        d2f(i,n2-6,k)  =   (                                               &
                        -   3.01824892312857d0      * f(i,n2-6,k)                       &
                        +   1.744086462579369d0     * (f(i,n2-5,k)+f(i,n2-7,k))         &
                        -   0.2851293540782d0      * (f(i,n2-4,k)+f(i,n2-8,k))         &
                        +   0.059207636829005d0      * (f(i,n2-3,k)+f(i,n2-9,k))         &
                        -   0.0099521462688359d0     * (f(i,n2-2,k)+f(i,n2-10,k))         &
                        +   9.1186250295123312d-4    * (f(i,n2-1,k) + f(i,n2-11,k))        &
                        ) /h**2
                    enddo
                enddo

            else

                do k=1, n3e
                    do i=1, n1e
                        d2f(i,n2-1,k)  =   ( -2.0d0*f(i,n2-1,k) + (f(i,n2,k)+f(i,n2-2,k)) ) / h**2
                        d2f(i,n2-2,k)  =   -5.0d0  * f(i,n2-2,k)  / (2.d0*h**2)   &
                        +   4.d0  * (f(i,n2-1,k)+f(i,n2-3,k)) / (3.d0*h**2)   &
                        -   1.d0 * (f(i,n2,k) + f(i,n2-4,k))    / (12.d0*h**2)
                        d2f(i,n2-3,k)  =   - 2.814728882213931d0 * f(i,n2-3,k)    / h**2      &
                        +   1.569379994993781d0*(f(i,n2-2,k)+f(i,n2-4,k))     / h**2      &
                        -   0.177751998d0*(f(i,n2-1,k)+f(i,n2-5,k))           / h**2      &
                        +   0.0157364441d0* (f(i,n2,k) + f(i,n2-6,k))           / h**2
                        d2f(i,n2-4,k)  =   - 2.814728882213931d0 * f(i,n2-4,k)    / h**2      &
                        +   1.569379994993781d0*(f(i,n2-3,k)+f(i,n2-5,k))     / h**2      &
                        -   0.177751998d0*(f(i,n2-2,k)+f(i,n2-6,k))           / h**2      &
                        +   0.0157364441d0* (f(i,n2-1,k) + f(i,n2-7,k))           / h**2
                        d2f(i,n2-5,k)  =   (                                               &
                        -   3.01824892312857d0      * f(i,n2-5,k)                       &
                        +   1.744086462579369d0     * (f(i,n2-4,k)+f(i,n2-6,k))         &
                        -   0.2851293540782d0      * (f(i,n2-3,k)+f(i,n2-7,k))         &
                        +   0.059207636829005d0      * (f(i,n2-2,k)+f(i,n2-8,k))         &
                        -   0.0099521462688359d0     * (f(i,n2-1,k)+f(i,n2-9,k))         &
                        +   9.1186250295123312d-4    * (f(i,n2,k) + f(i,n2-10,k))        &
                        ) /h**2
                        d2f(i,n2-6,k)  =   (                                           &
                        +   A(0)    * f(i,n2-6,k)                     &
                        +   A(1)     * (f(i,n2-5,k)+f(i,n2-7,k))           &
                        +   A(2)      * (f(i,n2-4,k)+f(i,n2-8,k))           &
                        +   A(3)      * (f(i,n2-3,k)+f(i,n2-9,k))         &
                        +   A(4)      * (f(i,n2-2,k)+f(i,n2-10,k))         &
                        +   A(5)     * (f(i,n2-1,k)+f(i,n2-11,k))       &
                        +   A(6)    * (f(i,n2,k) + f(i,n2-12,k))         &
                        )
                    enddo
                enddo
            end if

        endif

        return

    end subroutine D2c_DRP6_3Dy



    subroutine D2c_DRP6_MULT_3Dy(f, d2f, n1,n1e,n2,n3,n3e, h, shifted, bc_type, g)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: d2f
        real(kind=8), dimension(:), intent(in)                       :: g

        integer :: i,j,k
        real(kind=8) A(0:6)

        A(1) =  1.812906831737684d0         / h**2
        A(2) =  -0.33473276510983d0         / h**2
        A(3) =  0.087046253270279d0         / h**2
        A(4) =  -0.021406336650382d0        / h**2
        A(5) =  0.0039703382061228d0        / h**2
        A(6) = -3.9303109660650135d-4       / h**2
        A(0) = -2.d0*(A(1) + A(2) + A(3) + A(4) + A(5) + A(6))

        do k=1, n3e
            do i=1, n1e

                do j=7,n2-7
                    d2f(i,j,k)= (A(6)*(f(i,j+6,k) + f(i,j-6,k))     &
                    + A(5)*(f(i,j+5,k) + f(i,j-5,k))     &
                    + A(4)*(f(i,j+4,k) + f(i,j-4,k))     &
                    + A(3)*(f(i,j+3,k) + f(i,j-3,k))     &
                    + A(2)*(f(i,j+2,k) + f(i,j-2,k))     &
                    + A(1)*(f(i,j+1,k) + f(i,j-1,k))     &
                    + A(0)*f(i,j,k))*g(j)
                enddo

            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=1, 6
                        d2f(i,j,k)=   (A(6)*(f(i,j+6,k) + f(i,mod(n2+j-8,n2-1)+1,k))  &
                        + A(5)*(f(i,j+5,k) + f(i,mod(n2+j-7,n2-1)+1,k))  &
                        + A(4)*(f(i,j+4,k) + f(i,mod(n2+j-6,n2-1)+1,k))  &
                        + A(3)*(f(i,j+3,k) + f(i,mod(n2+j-5,n2-1)+1,k))  &
                        + A(2)*(f(i,j+2,k) + f(i,mod(n2+j-4,n2-1)+1,k))  &
                        + A(1)*(f(i,j+1,k) + f(i,mod(n2+j-3,n2-1)+1,k)) &
                        + A(0)*f(i,j,k))*g(j)
                    enddo

                    do j=n2-6, n2-1
                        d2f(i,j,k)= (A(6)*(f(i,mod(j+5,n2-1)+1,k) + f(i,j-6,k))             &
                        + A(5)*(f(i,mod(j+4,n2-1)+1,k) + f(i,j-5,k))             &
                        + A(4)*(f(i,mod(j+3,n2-1)+1,k) + f(i,j-4,k))             &
                        + A(3)*(f(i,mod(j+2,n2-1)+1,k) + f(i,j-3,k))             &
                        + A(2)*(f(i,mod(j+1,n2-1)+1,k) + f(i,j-2,k))             &
                        + A(1)*(f(i,mod(j,n2-1)+1,k) + f(i,j-1,k))               &
                        + A(0)*f(i,j,k))*g(j)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e
                    d2f(i,2,k)  =   (( -2.0d0*f(i,2,k) + (f(i,3,k)+f(i,1,k)) ) / h**2)*g(2)
                    d2f(i,3,k)  =   (-5.0d0  * f(i,3,k)  / (2.d0*h**2)   &
                    +   4.d0  * (f(i,4,k)+f(i,2,k)) / (3.d0*h**2)   &
                    -   1.d0 * (f(i,5,k) + f(i,1,k))    / (12.d0*h**2))*g(3)
                    d2f(i,4,k)  =   (- 2.814728882213931d0 * f(i,4,k)    / h**2      &
                    +   1.569379994993781d0*(f(i,5,k)+f(i,3,k))     / h**2      &
                    -   0.177751998d0*(f(i,6,k)+f(i,2,k))           / h**2      &
                    +   0.0157364441d0* (f(i,7,k) + f(i,1,k))           / h**2)*g(4)
                    d2f(i,5,k)  =   (- 2.814728882213931d0 * f(i,5,k)    / h**2      &
                    +   1.569379994993781d0*(f(i,6,k)+f(i,4,k))     / h**2      &
                    -   0.177751998d0*(f(i,7,k)+f(i,3,k))           / h**2      &
                    +   0.0157364441d0* (f(i,8,k) + f(i,2,k))           / h**2)*g(5)
                    d2f(i,6,k)  =   ((                                           &
                    -   3.01824892312857d0      * f(i,6,k)                      &
                    +   1.744086462579369d0     * (f(i,7,k)+f(i,5,k))           &
                    -   0.2851293540782d0      * (f(i,8,k)+f(i,4,k))           &
                    +   0.059207636829005d0      * (f(i,9,k) + f(i,3,k))         &
                    -   0.0099521462688359d0     * (f(i,10,k) + f(i,2,k))        &
                    +   9.1186250295123312d-4    * (f(i,11,k) + f(i,1,k))        &
                    ) /h**2)*g(6)
                enddo
            enddo



            if (shifted) then

                do k=1, n3e
                    do i=1, n1e
                        d2f(i,n2-2,k)  =   (( -2.0d0*f(i,n2-2,k) + (f(i,n2-1,k)+f(i,n2-3,k)) )  / h**2)*g(n2-2)
                        d2f(i,n2-3,k)  =   (-5.0d0  * f(i,n2-3,k)                               / (2.d0*h**2)   &
                        +   4.d0  * (f(i,n2-2,k)+f(i,n2-4,k))                                   / (3.d0*h**2)   &
                        -   1.d0 * (f(i,n2-1,k) + f(i,n2-5,k))                                  / (12.d0*h**2))*g(n2-3)
                        d2f(i,n2-4,k)  =   (- 2.814728882213931d0 * f(i,n2-4,k)                 / h**2      &
                        +   1.569379994993781d0*(f(i,n2-3,k)+f(i,n2-5,k))                       / h**2      &
                        -   0.177751998d0*(f(i,n2-2,k)+f(i,n2-6,k))                             / h**2      &
                        +   0.0157364441d0* (f(i,n2-1,k) + f(i,n2-7,k))                         / h**2)*g(n2-4)
                        d2f(i,n2-5,k)  =   (- 2.814728882213931d0 * f(i,n2-5,k)                 / h**2      &
                        +   1.569379994993781d0*(f(i,n2-4,k)+f(i,n2-6,k))                       / h**2      &
                        -   0.177751998d0*(f(i,n2-3,k)+f(i,n2-7,k))                             / h**2      &
                        +   0.0157364441d0* (f(i,n2-2,k) + f(i,n2-8,k))                         / h**2)*g(n2-5)
                        d2f(i,n2-6,k)  =   ((                                                               &
                        -   3.01824892312857d0      * f(i,n2-6,k)                                           &
                        +   1.744086462579369d0     * (f(i,n2-5,k)+f(i,n2-7,k))                             &
                        -   0.2851293540782d0      * (f(i,n2-4,k)+f(i,n2-8,k))                             &
                        +   0.059207636829005d0      * (f(i,n2-3,k)+f(i,n2-9,k))                             &
                        -   0.0099521462688359d0     * (f(i,n2-2,k)+f(i,n2-10,k))                            &
                        +   9.1186250295123312d-4    * (f(i,n2-1,k) + f(i,n2-11,k))                          &
                        )                                                                       /h**2)*g(n2-6)
                    enddo
                enddo

            else

                do k=1, n3e
                    do i=1, n1e
                        d2f(i,n2-1,k)  =   (( -2.0d0*f(i,n2-1,k) + (f(i,n2,k)+f(i,n2-2,k)) )    / h**2)*g(n2-1)
                        d2f(i,n2-2,k)  =   (-5.0d0  * f(i,n2-2,k)                               / (2.d0*h**2)   &
                        +   4.d0  * (f(i,n2-1,k)+f(i,n2-3,k))                                   / (3.d0*h**2)   &
                        -   1.d0 * (f(i,n2,k) + f(i,n2-4,k))                                    / (12.d0*h**2))*g(n2-2)
                        d2f(i,n2-3,k)  =   (- 2.814728882213931d0 * f(i,n2-3,k)                 / h**2      &
                        +   1.569379994993781d0*(f(i,n2-2,k)+f(i,n2-4,k))                       / h**2      &
                        -   0.177751998d0*(f(i,n2-1,k)+f(i,n2-5,k))                             / h**2      &
                        +   0.0157364441d0* (f(i,n2,k) + f(i,n2-6,k))                           / h**2)*g(n2-3)
                        d2f(i,n2-4,k)  =   (- 2.814728882213931d0 * f(i,n2-4,k)                 / h**2      &
                        +   1.569379994993781d0*(f(i,n2-3,k)+f(i,n2-5,k))                       / h**2      &
                        -   0.177751998d0*(f(i,n2-2,k)+f(i,n2-6,k))                             / h**2      &
                        +   0.0157364441d0* (f(i,n2-1,k) + f(i,n2-7,k))                         / h**2)*g(n2-4)
                        d2f(i,n2-5,k)  =   ((                                                               &
                        -   3.01824892312857d0      * f(i,n2-5,k)                                           &
                        +   1.744086462579369d0     * (f(i,n2-4,k)+f(i,n2-6,k))                             &
                        -   0.2851293540782d0      * (f(i,n2-3,k)+f(i,n2-7,k))                             &
                        +   0.059207636829005d0      * (f(i,n2-2,k)+f(i,n2-8,k))                             &
                        -   0.0099521462688359d0     * (f(i,n2-1,k)+f(i,n2-9,k))                             &
                        +   9.1186250295123312d-4    * (f(i,n2,k) + f(i,n2-10,k))                            &
                        )                                                                           /h**2)*g(n2-5)
                        d2f(i,n2-6,k)  =   ((                                                               &
                        +   A(0)    * f(i,n2-6,k)                                           &
                        +   A(1)     * (f(i,n2-5,k)+f(i,n2-7,k))                             &
                        +   A(2)      * (f(i,n2-4,k)+f(i,n2-8,k))                             &
                        +   A(3)      * (f(i,n2-3,k)+f(i,n2-9,k))                             &
                        +   A(4)      * (f(i,n2-2,k)+f(i,n2-10,k))                            &
                        +   A(5)     * (f(i,n2-1,k)+f(i,n2-11,k))                            &
                        +   A(6)    * (f(i,n2,k) + f(i,n2-12,k))                            &
                        )                         )*g(n2-6)
                    enddo
                enddo
            end if

        endif

        return

    end subroutine D2c_DRP6_MULT_3Dy



    subroutine D2c_DRP6_MULTACC_3Dy(f, d2f, n1,n1e,n2,n3,n3e, h, shifted, bc_type, g)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: d2f
        real(kind=8), dimension(:), intent(in)                       :: g

        integer :: i,j,k
        real(kind=8) A(0:6)

        A(1) =  1.812906831737684d0         / h**2
        A(2) =  -0.33473276510983d0         / h**2
        A(3) =  0.087046253270279d0         / h**2
        A(4) =  -0.021406336650382d0        / h**2
        A(5) =  0.0039703382061228d0        / h**2
        A(6) = -3.9303109660650135d-4       / h**2
        A(0) = -2.d0*(A(1) + A(2) + A(3) + A(4) + A(5) + A(6))

        do k=1, n3e
            do i=1, n1e

                do j=7,n2-7
                    d2f(i,j,k)= d2f(i,j,k)+(A(6)*(f(i,j+6,k) + f(i,j-6,k))     &
                    + A(5)*(f(i,j+5,k) + f(i,j-5,k))     &
                    + A(4)*(f(i,j+4,k) + f(i,j-4,k))     &
                    + A(3)*(f(i,j+3,k) + f(i,j-3,k))     &
                    + A(2)*(f(i,j+2,k) + f(i,j-2,k))     &
                    + A(1)*(f(i,j+1,k) + f(i,j-1,k))     &
                    + A(0)*f(i,j,k))*g(j)
                enddo

            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=1, 6
                        d2f(i,j,k)=   d2f(i,j,k)+(A(6)*(f(i,j+6,k) + f(i,mod(n2+j-8,n2-1)+1,k))  &
                        + A(5)*(f(i,j+5,k) + f(i,mod(n2+j-7,n2-1)+1,k))  &
                        + A(4)*(f(i,j+4,k) + f(i,mod(n2+j-6,n2-1)+1,k))  &
                        + A(3)*(f(i,j+3,k) + f(i,mod(n2+j-5,n2-1)+1,k))  &
                        + A(2)*(f(i,j+2,k) + f(i,mod(n2+j-4,n2-1)+1,k))  &
                        + A(1)*(f(i,j+1,k) + f(i,mod(n2+j-3,n2-1)+1,k)) &
                        + A(0)*f(i,j,k))*g(j)
                    enddo

                    do j=n2-6, n2-1
                        d2f(i,j,k)= d2f(i,j,k)+(A(6)*(f(i,mod(j+5,n2-1)+1,k) + f(i,j-6,k))             &
                        + A(5)*(f(i,mod(j+4,n2-1)+1,k) + f(i,j-5,k))             &
                        + A(4)*(f(i,mod(j+3,n2-1)+1,k) + f(i,j-4,k))             &
                        + A(3)*(f(i,mod(j+2,n2-1)+1,k) + f(i,j-3,k))             &
                        + A(2)*(f(i,mod(j+1,n2-1)+1,k) + f(i,j-2,k))             &
                        + A(1)*(f(i,mod(j,n2-1)+1,k) + f(i,j-1,k))               &
                        + A(0)*f(i,j,k))*g(j)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e
                    d2f(i,2,k)  =   d2f(i,2,k)+(( -2.0d0*f(i,2,k) + (f(i,3,k)+f(i,1,k)) ) / h**2)*g(2)
                    d2f(i,3,k)  =   d2f(i,3,k)+(-5.0d0  * f(i,3,k)  / (2.d0*h**2)   &
                    +   4.d0  * (f(i,4,k)+f(i,2,k)) / (3.d0*h**2)   &
                    -   1.d0 * (f(i,5,k) + f(i,1,k))    / (12.d0*h**2))*g(3)
                    d2f(i,4,k)  =   d2f(i,4,k)+(- 2.814728882213931d0 * f(i,4,k)    / h**2      &
                    +   1.569379994993781d0*(f(i,5,k)+f(i,3,k))     / h**2      &
                    -   0.177751998d0*(f(i,6,k)+f(i,2,k))           / h**2      &
                    +   0.0157364441d0* (f(i,7,k) + f(i,1,k))           / h**2)*g(4)
                    d2f(i,5,k)  =   d2f(i,5,k)+(- 2.814728882213931d0 * f(i,5,k)    / h**2      &
                    +   1.569379994993781d0*(f(i,6,k)+f(i,4,k))     / h**2      &
                    -   0.177751998d0*(f(i,7,k)+f(i,3,k))           / h**2      &
                    +   0.0157364441d0* (f(i,8,k) + f(i,2,k))           / h**2)*g(5)
                    d2f(i,6,k)  =   d2f(i,6,k)+((                                           &
                    -   3.01824892312857d0      * f(i,6,k)                      &
                    +   1.744086462579369d0     * (f(i,7,k)+f(i,5,k))           &
                    -   0.2851293540782d0      * (f(i,8,k)+f(i,4,k))           &
                    +   0.059207636829005d0      * (f(i,9,k) + f(i,3,k))         &
                    -   0.0099521462688359d0     * (f(i,10,k) + f(i,2,k))        &
                    +   9.1186250295123312d-4    * (f(i,11,k) + f(i,1,k))        &
                    ) /h**2)*g(6)
                enddo
            enddo

            if (shifted) then

                do k=1, n3e
                    do i=1, n1e
                        d2f(i,n2-2,k)  =   d2f(i,n2-2,k)+(( -2.0d0*f(i,n2-2,k) + (f(i,n2-1,k)+f(i,n2-3,k)) )  / h**2)*g(n2-2)
                        d2f(i,n2-3,k)  =   d2f(i,n2-3,k)+(-5.0d0  * f(i,n2-3,k)                               / (2.d0*h**2)   &
                        +   4.d0  * (f(i,n2-2,k)+f(i,n2-4,k))                                   / (3.d0*h**2)   &
                        -   1.d0 * (f(i,n2-1,k) + f(i,n2-5,k))                                  / (12.d0*h**2))*g(n2-3)
                        d2f(i,n2-4,k)  =   d2f(i,n2-4,k)+(- 2.814728882213931d0 * f(i,n2-4,k)                 / h**2      &
                        +   1.569379994993781d0*(f(i,n2-3,k)+f(i,n2-5,k))                       / h**2      &
                        -   0.177751998d0*(f(i,n2-2,k)+f(i,n2-6,k))                             / h**2      &
                        +   0.0157364441d0* (f(i,n2-1,k) + f(i,n2-7,k))                         / h**2)*g(n2-4)
                        d2f(i,n2-5,k)  =   d2f(i,n2-5,k)+(- 2.814728882213931d0 * f(i,n2-5,k)                 / h**2      &
                        +   1.569379994993781d0*(f(i,n2-4,k)+f(i,n2-6,k))                       / h**2      &
                        -   0.177751998d0*(f(i,n2-3,k)+f(i,n2-7,k))                             / h**2      &
                        +   0.0157364441d0* (f(i,n2-2,k) + f(i,n2-8,k))                         / h**2)*g(n2-5)
                        d2f(i,n2-6,k)  =   d2f(i,n2-6,k)+((                                                               &
                        -   3.01824892312857d0      * f(i,n2-6,k)                                           &
                        +   1.744086462579369d0     * (f(i,n2-5,k)+f(i,n2-7,k))                             &
                        -   0.2851293540782d0      * (f(i,n2-4,k)+f(i,n2-8,k))                             &
                        +   0.059207636829005d0      * (f(i,n2-3,k)+f(i,n2-9,k))                             &
                        -   0.0099521462688359d0     * (f(i,n2-2,k)+f(i,n2-10,k))                            &
                        +   9.1186250295123312d-4    * (f(i,n2-1,k) + f(i,n2-11,k))                          &
                        )                                                                       /h**2)*g(n2-6)
                    enddo
                enddo

            else

                do k=1, n3e
                    do i=1, n1e
                        d2f(i,n2-1,k)  =   d2f(i,n2-1,k)+(( -2.0d0*f(i,n2-1,k) + (f(i,n2,k)+f(i,n2-2,k)) )    / h**2)*g(n2-1)
                        d2f(i,n2-2,k)  =   d2f(i,n2-2,k)+(-5.0d0  * f(i,n2-2,k)                               / (2.d0*h**2)   &
                        +   4.d0  * (f(i,n2-1,k)+f(i,n2-3,k))                                   / (3.d0*h**2)   &
                        -   1.d0 * (f(i,n2,k) + f(i,n2-4,k))                                    / (12.d0*h**2))*g(n2-2)
                        d2f(i,n2-3,k)  =   d2f(i,n2-3,k)+(- 2.814728882213931d0 * f(i,n2-3,k)                 / h**2      &
                        +   1.569379994993781d0*(f(i,n2-2,k)+f(i,n2-4,k))                       / h**2      &
                        -   0.177751998d0*(f(i,n2-1,k)+f(i,n2-5,k))                             / h**2      &
                        +   0.0157364441d0* (f(i,n2,k) + f(i,n2-6,k))                           / h**2)*g(n2-3)
                        d2f(i,n2-4,k)  =   d2f(i,n2-4,k)+(- 2.814728882213931d0 * f(i,n2-4,k)                 / h**2      &
                        +   1.569379994993781d0*(f(i,n2-3,k)+f(i,n2-5,k))                       / h**2      &
                        -   0.177751998d0*(f(i,n2-2,k)+f(i,n2-6,k))                             / h**2      &
                        +   0.0157364441d0* (f(i,n2-1,k) + f(i,n2-7,k))                         / h**2)*g(n2-4)
                        d2f(i,n2-5,k)  =   d2f(i,n2-5,k)+((                                                               &
                        -   3.01824892312857d0      * f(i,n2-5,k)                                           &
                        +   1.744086462579369d0     * (f(i,n2-4,k)+f(i,n2-6,k))                             &
                        -   0.2851293540782d0      * (f(i,n2-3,k)+f(i,n2-7,k))                             &
                        +   0.059207636829005d0      * (f(i,n2-2,k)+f(i,n2-8,k))                             &
                        -   0.0099521462688359d0     * (f(i,n2-1,k)+f(i,n2-9,k))                             &
                        +   9.1186250295123312d-4    * (f(i,n2,k) + f(i,n2-10,k))                            &
                        )                                                                           /h**2)*g(n2-5)
                        d2f(i,n2-6,k)  =   d2f(i,n2-6,k)+((                                                               &
                        +   A(0)    * f(i,n2-6,k)                                           &
                        +   A(1)     * (f(i,n2-5,k)+f(i,n2-7,k))                             &
                        +   A(2)      * (f(i,n2-4,k)+f(i,n2-8,k))                             &
                        +   A(3)      * (f(i,n2-3,k)+f(i,n2-9,k))                             &
                        +   A(4)      * (f(i,n2-2,k)+f(i,n2-10,k))                            &
                        +   A(5)     * (f(i,n2-1,k)+f(i,n2-11,k))                            &
                        +   A(6)    * (f(i,n2,k) + f(i,n2-12,k))                            &
                        )                                   )*g(n2-6)
                    enddo
                enddo

            end if

        endif

        return

    end subroutine D2c_DRP6_MULTACC_3Dy



    subroutine D2c_DRP6_3Dz(f, d2f, n1,n1e,n2,n2e,n3, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: d2f

        integer :: i,j,k
        real(kind=8) A(0:6)

        A(1) =  1.812906831737684d0         / h**2
        A(2) =  -0.33473276510983d0         / h**2
        A(3) =  0.087046253270279d0         / h**2
        A(4) =  -0.021406336650382d0        / h**2
        A(5) =  0.0039703382061228d0        / h**2
        A(6) = -3.9303109660650135d-4       / h**2
        A(0) = -2.d0*(A(1) + A(2) + A(3) + A(4) + A(5) + A(6))

        do j=1, n2e
            do i=1, n1e
                do k=7,n3-7
                    d2f(i,j,k)= A(6)*(f(i,j,k+6) + f(i,j,k-6))     &
                    + A(5)*(f(i,j,k+5) + f(i,j,k-5))        &
                    + A(4)*(f(i,j,k+4) + f(i,j,k-4))        &
                    + A(3)*(f(i,j,k+3) + f(i,j,k-3))        &
                    + A(2)*(f(i,j,k+2) + f(i,j,k-2))        &
                    + A(1)*(f(i,j,k+1) + f(i,j,k-1))        &
                    + A(0)*f(i,j,k)
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e

                    do k=1, 6
                        d2f(i,j,k)=   A(6)*(f(i,j,k+6) + f(i,j,mod(n3+k-8,n3-1)+1))  &
                        + A(5)*(f(i,j,k+5) + f(i,j,mod(n3+k-7,n3-1)+1))      &
                        + A(4)*(f(i,j,k+4) + f(i,j,mod(n3+k-6,n3-1)+1))     &
                        + A(3)*(f(i,j,k+3) + f(i,j,mod(n3+k-5,n3-1)+1))     &
                        + A(2)*(f(i,j,k+2) + f(i,j,mod(n3+k-4,n3-1)+1))     &
                        + A(1)*(f(i,j,k+1) + f(i,j,mod(n3+k-3,n3-1)+1))     &
                        + A(0)*f(i,j,k)
                    enddo

                    do k=n3-6, n3-1
                        d2f(i,j,k)= A(6)*(f(i,j,mod(k+5,n3-1)+1) + f(i,j,k-6))             &
                        + A(5)*(f(i,j,mod(k+4,n3-1)+1) + f(i,j,k-5))            &
                        + A(4)*(f(i,j,mod(k+3,n3-1)+1) + f(i,j,k-4))            &
                        + A(3)*(f(i,j,mod(k+2,n3-1)+1) + f(i,j,k-3))            &
                        + A(2)*(f(i,j,mod(k+1,n3-1)+1) + f(i,j,k-2))            &
                        + A(1)*(f(i,j,mod(k,n3-1)+1) + f(i,j,k-1))              &
                        + A(0)*f(i,j,k)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e
                    d2f(i,j,2)  =   ( -2.0d0*f(i,j,2) + (f(i,j,3)+f(i,j,1)) ) / h**2
                    d2f(i,j,3)  =   -5.0d0  * f(i,j,3)  / (2.d0*h**2)   &
                    +   4.d0  * (f(i,j,4)+f(i,j,2)) / (3.d0*h**2)   &
                    -   1.d0 * (f(i,j,5) + f(i,j,1))    / (12.d0*h**2)
                    d2f(i,j,4)  =   - 2.814728882213931d0 * f(i,j,4)    / h**2      &
                    +   1.569379994993781d0*(f(i,j,5)+f(i,j,3))     / h**2      &
                    -   0.177751998d0*(f(i,j,6)+f(i,j,2))           / h**2      &
                    +   0.0157364441d0* (f(i,j,7) + f(i,j,1))           / h**2
                    d2f(i,j,5)  =   - 2.814728882213931d0 * f(i,j,5)    / h**2      &
                    +   1.569379994993781d0*(f(i,j,6)+f(i,j,4))     / h**2      &
                    -   0.177751998d0*(f(i,j,7)+f(i,j,3))           / h**2      &
                    +   0.0157364441d0* (f(i,j,8) + f(i,j,2))           / h**2
                    d2f(i,j,6)  =   (                                           &
                    -   3.01824892312857d0      * f(i,j,6)                      &
                    +   1.744086462579369d0     * (f(i,j,7)+f(i,j,5))           &
                    -   0.2851293540782d0      * (f(i,j,8)+f(i,j,4))           &
                    +   0.059207636829005d0      * (f(i,j,9) + f(i,j,3))         &
                    -   0.0099521462688359d0     * (f(i,j,10) + f(i,j,2))        &
                    +   9.1186250295123312d-4    * (f(i,j,11) + f(i,j,1))        &
                    ) /h**2
                enddo
            enddo

            if (shifted) then

                do j=1, n2e
                    do i=1, n1e
                        d2f(i,j,n3-2)  =   ( -2.0d0*f(i,j,n3-2) + (f(i,j,n3-1)+f(i,j,n3-3)) ) / h**2
                        d2f(i,j,n3-3)  =   -5.0d0  * f(i,j,n3-3)  / (2.d0*h**2)   &
                        +   4.d0  * (f(i,j,n3-2)+f(i,j,n3-4)) / (3.d0*h**2)   &
                        -   1.d0 * (f(i,j,n3-1) + f(i,j,n3-5))    / (12.d0*h**2)
                        d2f(i,j,n3-4)  =   - 2.814728882213931d0 * f(i,j,n3-4)    / h**2      &
                        +   1.569379994993781d0*(f(i,j,n3-3)+f(i,j,n3-5))     / h**2      &
                        -   0.177751998d0*(f(i,j,n3-2)+f(i,j,n3-6))           / h**2      &
                        +   0.0157364441d0* (f(i,j,n3-1) + f(i,j,n3-7))           / h**2
                        d2f(i,j,n3-5)  =   - 2.814728882213931d0 * f(i,j,n3-5)    / h**2      &
                        +   1.569379994993781d0*(f(i,j,n3-4)+f(i,j,n3-6))     / h**2      &
                        -   0.177751998d0*(f(i,j,n3-3)+f(i,j,n3-7))           / h**2      &
                        +   0.0157364441d0* (f(i,j,n3-2) + f(i,j,n3-8))           / h**2
                        d2f(i,j,n3-6)  =   (                                               &
                        -   3.01824892312857d0      * f(i,j,n3-6)                  &
                        +   1.744086462579369d0     * (f(i,j,n3-5)+f(i,j,n3-7))       &
                        -   0.2851293540782d0      * (f(i,j,n3-4)+f(i,j,n3-8))       &
                        +   0.059207636829005d0      * (f(i,j,n3-3)+f(i,j,n3-9))     &
                        -   0.0099521462688359d0     * (f(i,j,n3-2)+f(i,j,n3-10))     &
                        +   9.1186250295123312d-4    * (f(i,j,n3-1) + f(i,j,n3-11))     &
                        ) /h**2
                    enddo
                enddo

            else

                do j=1, n2e
                    do i=1, n1e
                        d2f(i,j,n3-1)  =   ( -2.0d0*f(i,j,n3-1) + (f(i,j,n3)+f(i,j,n3-2)) ) / h**2
                        d2f(i,j,n3-2)  =   -5.0d0  * f(i,j,n3-2)  / (2.d0*h**2)   &
                        +   4.d0  * (f(i,j,n3-1)+f(i,j,n3-3)) / (3.d0*h**2)   &
                        -   1.d0 * (f(i,j,n3) + f(i,j,n3-4))    / (12.d0*h**2)
                        d2f(i,j,n3-3)  =   - 2.814728882213931d0 * f(i,j,n3-3)    / h**2      &
                        +   1.569379994993781d0*(f(i,j,n3-2)+f(i,j,n3-4))     / h**2      &
                        -   0.177751998d0*(f(i,j,n3-1)+f(i,j,n3-5))           / h**2      &
                        +   0.0157364441d0* (f(i,j,n3) + f(i,j,n3-6))           / h**2
                        d2f(i,j,n3-4)  =   - 2.814728882213931d0 * f(i,j,n3-4)    / h**2      &
                        +   1.569379994993781d0*(f(i,j,n3-3)+f(i,j,n3-5))     / h**2      &
                        -   0.177751998d0*(f(i,j,n3-2)+f(i,j,n3-6))           / h**2      &
                        +   0.0157364441d0* (f(i,j,n3-1) + f(i,j,n3-7))           / h**2
                        d2f(i,j,n3-5)  =   (                                               &
                        -   3.01824892312857d0      * f(i,j,n3-5)                  &
                        +   1.744086462579369d0     * (f(i,j,n3-4)+f(i,j,n3-6))       &
                        -   0.2851293540782d0      * (f(i,j,n3-3)+f(i,j,n3-7))       &
                        +   0.059207636829005d0      * (f(i,j,n3-2)+f(i,j,n3-8))     &
                        -   0.0099521462688359d0     * (f(i,j,n3-1)+f(i,j,n3-9))     &
                        +   9.1186250295123312d-4    * (f(i,j,n3) + f(i,j,n3-10))     &
                        ) /h**2
                        d2f(i,j,n3-6)  =   (                                           &
                        +   A(0)    * f(i,j,n3-6)                     &
                        +   A(1)     * (f(i,j,n3-5)+f(i,j,n3-7))           &
                        +   A(2)      * (f(i,j,n3-4)+f(i,j,n3-8))           &
                        +   A(3)      * (f(i,j,n3-3)+f(i,j,n3-9))         &
                        +   A(4)      * (f(i,j,n3-2)+f(i,j,n3-10))         &
                        +   A(5)     * (f(i,j,n3-1)+f(i,j,n3-11))       &
                        +   A(6)    * (f(i,j,n3) + f(i,j,n3-12))         &
                        )
                    enddo
                enddo

            end if

        endif

        return

    end subroutine D2c_DRP6_3Dz

end module DRP_D2c
