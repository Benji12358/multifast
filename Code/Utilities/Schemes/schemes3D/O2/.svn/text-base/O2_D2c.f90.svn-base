module O2_D2c
    use boundaries_types
    implicit none

contains

    subroutine D2c_O2_3Dx(f, d2f, n1,n2,n2e,n3,n3e, h, shifted, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: d2f

        integer :: i,j,k
        real(kind=8) A(0:1)

        A(1) =  1.d0         / h**2
        A(0) = -2.d0        /h**2
        do k=1,n3e
            do j=1,n2e
                do i=2,n1-2
                    d2f(i, j, k)= A(1)*(f(i+1, j, k) + f(i-1, j, k))     &
                    + A(0)*f(i, j, k)
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do k=1,n3e
                do j=1,n2e
                    d2f(1, j, k)=   A(1)*(f(2, j, k) + f(n1-1, j, k))  &
                    + A(0)*f(1, j, k)
                    d2f(n1-1, j, k)=   A(1)*(f(1, j, k) + f(n1-2, j, k))  &
                    + A(0)*f(n1-1, j, k)
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            if (.not. shifted) then

                do k=1,n3e
                    do j=1,n2e
                        d2f(n1-1,j,k)  =   ( -2.0d0*f(n1-1,j,k) + (f(n1,j,k)+f(n1-2,j,k)) )     / h**2
                    enddo
                enddo

            end if

        endif

        return

    end subroutine D2c_O2_3Dx

    subroutine D2c_O2_ACC_3Dx(f, d2f, n1,n2,n2e,n3,n3e, h, shifted, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: d2f

        integer :: i,j,k
        real(kind=8) A(0:1)

        A(1) =  1.d0    / h**2
        A(0) = -2.d0*(A(1))

        do k=1,n3e
            do j=1,n2e
                do i=2,n1-2
                    d2f(i, j, k)= d2f(i, j, k)+(A(1)*(f(i+1, j, k) + f(i-1, j, k))     &
                    + A(0)*f(i, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do k=1,n3e
                do j=1,n2e

                    d2f(1, j, k)= d2f(1, j, k)+(A(1)*(f(2, j, k) + f(n1-1, j, k))               &
                        + A(0)*f(1, j, k))
                    d2f(n1-1, j, k)= d2f(n1-1, j, k)+(A(1)*(f(1, j, k) + f(n1-2, j, k))               &
                        + A(0)*f(n1-1, j, k))
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            if (.not. shifted) then

                do k=1,n3e
                    do j=1,n2e
                        d2f(n1-1,j,k)  =   d2f(n1-1,j,k)+( -2.0d0*f(n1-1,j,k) + (f(n1,j,k)+f(n1-2,j,k)) )     / h**2
                    enddo
                enddo

            end if

        endif

        return

    end subroutine D2c_O2_ACC_3Dx



    subroutine D2c_O2_3Dy(f, d2f, n1,n1e,n2,n3,n3e, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: d2f

        integer :: i,j,k
        real(kind=8) A(0:1)

        A(1) =  1.d0        / h**2
        A(0) = -2.d0        / h**2

        do k=1, n3e
            do i=1, n1e

                do j=2,n2-2
                    d2f(i,j,k)= A(1)*(f(i,j+1,k) + f(i,j-1,k))     &
                    + A(0)*f(i,j,k)
                enddo

            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e
                    d2f(i,1,k)= A(1)*(f(i,2,k) + f(i,n2-1,k))               &
                        + A(0)*f(i,1,k)
                    d2f(i,n2-1,k)= A(1)*(f(i,1,k) + f(i,n2-2,k))               &
                        + A(0)*f(i,n2-1,k)
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            if (.not. shifted) then

                do k=1, n3e
                    do i=1, n1e
                        d2f(i,n2-1,k)  =   ( -2.0d0*f(i,n2-1,k) + (f(i,n2,k)+f(i,n2-2,k)) ) / h**2
                    enddo
                enddo

            end if

        endif

        return

    end subroutine D2c_O2_3Dy



    subroutine D2c_O2_MULT_3Dy(f, d2f, n1,n1e,n2,n3,n3e, h, shifted, bc_type, g)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: d2f
        real(kind=8), dimension(:), intent(in)                       :: g

        integer :: i,j,k
        real(kind=8) A(0:1)

        A(1) =  1.0d0         / h**2
        A(0) = -2.d0*(A(1))

        do k=1, n3e
            do i=1, n1e

                do j=2,n2-2
                    d2f(i,j,k)= (A(1)*(f(i,j+1,k) + f(i,j-1,k))     &
                    + A(0)*f(i,j,k))*g(j)
                enddo

            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    d2f(i,1,k)= (A(1)*(f(i,2,k) + f(i,n2-1,k))               &
                        + A(0)*f(i,1,k))*g(1)
                    d2f(i,n2-1,k)= (A(1)*(f(i,1,k) + f(i,n2-2,k))               &
                        + A(0)*f(i,n2-1,k))*g(n2-1)

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            if (.not. shifted) then

                do k=1, n3e
                    do i=1, n1e
                        d2f(i,n2-1,k)  =   (( -2.0d0*f(i,n2-1,k) + (f(i,n2,k)+f(i,n2-2,k)) )    / h**2)*g(n2-1)
                    enddo
                enddo

            end if

        endif

        return

    end subroutine D2c_O2_MULT_3Dy



    subroutine D2c_O2_MULTACC_3Dy(f, d2f, n1,n1e,n2,n3,n3e, h, shifted, bc_type, g)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: d2f
        real(kind=8), dimension(:), intent(in)                       :: g

        integer :: i,j,k
        real(kind=8) A(0:1)

        A(1) =  1.d0        / h**2
        A(0) = -2.d0*(A(1))

        do k=1, n3e
            do i=1, n1e

                do j=2,n2-2
                    d2f(i,j,k)= d2f(i,j,k)+(A(1)*(f(i,j+1,k) + f(i,j-1,k))     &
                    + A(0)*f(i,j,k))*g(j)
                enddo

            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    d2f(i,1,k)= d2f(i,1,k)+(A(1)*(f(i,2,k) + f(i,n2-1,k))               &
                        + A(0)*f(i,1,k))*g(1)
                    d2f(i,n2-1,k)= d2f(i,n2-1,k)+(A(1)*(f(i,1,k) + f(i,n2-2,k))               &
                        + A(0)*f(i,n2-1,k))*g(n2-1)

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then


            if (.not. shifted) then

                do k=1, n3e
                    do i=1, n1e
                        d2f(i,n2-1,k)  =   d2f(i,n2-1,k)+(( -2.0d0*f(i,n2-1,k) + (f(i,n2,k)+f(i,n2-2,k)) )    / h**2)*g(n2-1)
                    enddo
                enddo

            end if

        endif

        return

    end subroutine D2c_O2_MULTACC_3Dy



    subroutine D2c_O2_3Dz(f, d2f, n1,n1e,n2,n2e,n3, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: d2f

        integer :: i,j,k
        real(kind=8) A(0:1)

        A(1) =  1.d0         / h**2
        A(0) = -2.d0*(A(1))

        do j=1, n2e
            do i=1, n1e
                do k=2,n3-2
                    d2f(i,j,k)= A(1)*(f(i,j,k+1) + f(i,j,k-1))        &
                    + A(0)*f(i,j,k)
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e
                    d2f(i,j,1)=   A(1)*(f(i,j,2) + f(i,j,n3-1))     &
                        + A(0)*f(i,j,1)
                    d2f(i,j,n3-1)=   A(1)*(f(i,j,1) + f(i,j,n3-2))     &
                        + A(0)*f(i,j,n3-1)
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            if (.not. shifted) then

                do j=1, n2e
                    do i=1, n1e
                        d2f(i,j,n3-1)  =   ( -2.0d0*f(i,j,n3-1) + (f(i,j,n3)+f(i,j,n3-2)) ) / h**2
                    enddo
                enddo

            end if

        endif

        return

    end subroutine D2c_O2_3Dz

end module O2_D2c
