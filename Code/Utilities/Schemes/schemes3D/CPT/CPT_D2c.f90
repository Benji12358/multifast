module CPT_D2c
    use boundaries_types
    implicit none

    private

    real(kind=8), allocatable, dimension(:) :: q,v
    real(kind=8), allocatable, dimension(:) :: m, bp, cp
    real(kind=8), dimension(:), allocatable :: A1, A2, A3
    real(kind=8)                                 :: vy, vq, coef_q
    real(kind=8)    :: B(0:2), Bl(-1:3)

    public::D2c_CPT_3Dx,D2c_CPT_ACC_3Dx, D2c_CPT_3Dy,D2c_CPT_MULT_3Dy,D2c_CPT_MULTACC_3Dy, D2c_CPT_3Dz

contains

    subroutine fill_matrix(n_max, shifted, bc_type, h)
        use CPT_INIT
        implicit none
        integer     :: n_max
        logical     :: shifted
        integer     :: bc_type
        real*8      :: h
        integer     :: i

        if (allocated(A1)) then
            deallocate(A1)
            deallocate(A2)
            deallocate(A3)
        end if

        allocate(A1(n_max))
        allocate(A2(n_max))
        allocate(A3(n_max))

        B(0)=-102.d0    / (44.d0*h**2)
        B(1)=12.d0      / (11.d0*h**2)
        B(2)= 3.d0      / (44.d0*h**2)

        do i=1, n_max-1
            A1(i) = 2.d0/11.d0
            A2(i) = 1.d0
            A3(i) = 2.d0/11.d0
        enddo

        select case (bc_type)

            case (periodic)

                if (allocated(q)) then
                    deallocate(q)
                    deallocate(v)
                    deallocate(m, bp, cp)
                end if

                allocate(q(n_max-1))
                allocate(v(n_max-1))
                allocate(m(n_max-1), bp(n_max-1), cp(n_max-1))

                call TS_initPR(A1(1:n_max-1), A2(1:n_max-1), A3(1:n_max-1), bp, cp, m, q,v,vq, n_max-1)

            case (Dirichlet)

                A1(2)=0.d0
                A2(2)=1.d0
                A3(2)=-1.d0

                Bl(-1)   =   1.d0    / h**2
                Bl(0)    =   -3.d0   / h**2
                Bl(1)    =   3.d0    / h**2
                Bl(2)    =   -1.d0   / h**2

                if (shifted) then

                    A1(n_max-2)=-1.d0
                    A2(n_max-2)=1.d0
                    A3(n_max-2)=0.d0


                    if (allocated(m)) then
                        deallocate(m, bp)
                    end if

                    allocate(m(2:n_max-2), bp(2:n_max-2))

                    call TS_initDIR(A1(2:n_max-2), A2(2:n_max-2), A3(2:n_max-2), bp(2:n_max-2), m(2:n_max-2), n_max-3)

                else

                    A1(n_max-2) = 2.d0/11.d0
                    A2(n_max-2) = 1.d0
                    A3(n_max-2) = 2.d0/11.d0

                    A1(n_max-1)=-1.d0
                    A2(n_max-1)=1.d0
                    A3(n_max-1)=0.d0

                    if (allocated(m)) then
                        deallocate(m, bp)
                    end if

                    allocate(m(2:n_max-1), bp(2:n_max-1))

                    call TS_initDIR(A1(2:n_max-1), A2(2:n_max-1), A3(2:n_max-1), bp(2:n_max-1), m(2:n_max-1), n_max-2)
                end if

            case default

        end select

    end subroutine fill_matrix

    subroutine D2c_CPT_3Dx(f, d2f, n1, n2,n2e, n3,n3e, h, shifted, bc_type)

        implicit none


        integer, intent(in) :: n1, n2,n2e, n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: d2f

        real(kind=8), dimension(n2,n3)                      :: sx

        integer :: i, j, k
        real(kind=8) rhs(n1)


        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n1, shifted, bc_type, h)

        if (bc_type.eq.periodic) then

            do k=1,n3e
                do j=1,n2e
                    d2f(1, j,k)=     B(0)*f(1, j,k) + B(1)*( f(2, j,k)+f(n1-1, j,k) )  &
                    +   B(2)*( f(3, j,k)+f(n1-2, j,k) )
                    d2f(2, j,k)=     B(0)*f(2, j,k) + B(1)*( f(3, j,k)+f(1, j,k) )  &
                    +   B(2)*( f(4, j,k)+f(n1-1, j,k) )
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i=3, n1-3
                        d2f(i, j,k)=     B(0)*f(i,j,k) + B(2)*( f(i+2, j,k)+f(i-2, j,k))   &
                        +   B(1)*( f(i+1, j,k)+f(i-1, j,k) )
                    enddo
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    d2f(n1-1, j,k)=   B(0)*f(n1-1, j,k) + B(1)*( f(1, j,k)+f(n1-2, j,k) )  &
                    +   B(2)*( f(2, j,k)+f(n1-3, j,k) )
                    d2f(n1-2, j,k)=   B(0)*f(n1-2, j,k) + B(1)*( f(n1-1, j,k)+f(n1-3, j,k) )  &
                    +   B(2)*( f(1, j,k)+f(n1-4, j,k) )
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i = 2,n1-1
                        d2f(i, j,k)=d2f(i, j,k) - m(i)*d2f(i-1, j,k)
                    enddo
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    d2f(n1-1, j,k) = d2f(n1-1, j,k)*bp(n1-1)
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i = n1-2, 1, -1
                        d2f(i, j,k)=(d2f(i, j,k) - A3(i)*d2f(i+1, j,k))*bp(i)
                    enddo
                enddo
            enddo

            do k=1, n3e
                do j=1, n2e
                    sx(j,k) = (v(1)*d2f(1,j,k)+v(n1-1)*d2f(n1-1,j,k))/(1+v(1)*q(1)+v(n1-1)*q(n1-1))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i=1, n1-1
                        d2f(i, j,k)=d2f(i, j,k)-sx(j,k)*q(i)
                    enddo
                enddo
            enddo


        end if


        if ((bc_type.eq.Dirichlet)) then

            if (shifted) then

                do k=1,n3e
                    do j=1,n2e

                        d2f(2, j,k)= Bl(-1)*f(1, j,k) + Bl(0)*f(2, j,k) + Bl(1)*f(3, j,k) + Bl(2)*f(4, j,k)

                        do i=3, n1-3
                            d2f(i, j,k)=     B(0)*f(i, j,k) + B(2)*( f(i+2, j,k)+f(i-2, j,k))   &
                            +   B(1)*( f(i+1, j,k)+f(i-1, j,k) )
                        enddo

                        d2f(n1-2, j,k)= Bl(-1)*f(n1-1, j,k) + Bl(0)*f(n1-2, j,k)     &
                        + Bl(1)*f(n1-3, j,k)+ Bl(2)*f(n1-4, j,k)


                        do i=3, n1-3
                            d2f(i, j,k)= d2f(i, j,k) - m(i)*d2f(i-1, j,k)
                        enddo

                        d2f(n1-2, j,k)= (d2f(n1-2, j,k) - m(n1-2)*d2f(n1-3, j,k))*bp(n1-2)


                        do i = n1-3, 2, -1
                            d2f(i, j,k) = (d2f(i, j,k) - A3(i)*d2f(i+1, j,k))*bp(i)
                        enddo

                    enddo
                enddo

            else

                do k=1,n3e
                    do j=1,n2e

                        d2f(2, j,k)= Bl(-1)*f(1, j,k) + Bl(0)*f(2, j,k) + Bl(1)*f(3, j,k) + Bl(2)*f(4, j,k)

                        do i=3, n1-2
                            d2f(i, j,k)=     B(0)*f(i, j,k) + B(2)*( f(i+2, j,k)+f(i-2, j,k))   &
                            +   B(1)*( f(i+1, j,k)+f(i-1, j,k) )
                        enddo

                        d2f(n1-1, j,k)= Bl(-1)*f(n1, j,k) + Bl(0)*f(n1-1, j,k)     &
                        + Bl(1)*f(n1-2, j,k)+ Bl(2)*f(n1-3, j,k)


                        do i=3, n1-2
                            d2f(i, j,k)= d2f(i, j,k) - m(i)*d2f(i-1, j,k)
                        enddo

                        d2f(n1-1, j,k)= (d2f(n1-1, j,k) - m(n1-1)*d2f(n1-2, j,k))*bp(n1-1)


                        do i = n1-2, 2, -1
                            d2f(i, j,k) = (d2f(i, j,k) - A3(i)*d2f(i+1, j,k))*bp(i)
                        enddo

                    enddo
                enddo

            end if
        end if

        return

    end subroutine D2c_CPT_3Dx

    subroutine D2c_CPT_ACC_3Dx(f, d2f, n1, n2,n2e, n3,n3e, h, shifted, bc_type)

        implicit none


        integer, intent(in) :: n1, n2,n2e, n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: d2f

        real*8, dimension(n1)                   :: tmp
        real(kind=8), dimension(n1,n2,n3)       :: tmp2

        real(kind=8), dimension(n2,n3)          :: sx

        integer :: i, j, k
        real(kind=8) rhs(n1)


        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n1, shifted, bc_type, h)

        if (bc_type.eq.periodic) then

            do k=1,n3e
                do j=1,n2e
                    tmp2(1, j,k)=     B(0)*f(1, j,k) + B(1)*( f(2, j,k)+f(n1-1, j,k) )  &
                    +   B(2)*( f(3, j,k)+f(n1-2, j,k) )
                    tmp2(2, j,k)=     B(0)*f(2, j,k) + B(1)*( f(3, j,k)+f(1, j,k) )  &
                    +   B(2)*( f(4, j,k)+f(n1-1, j,k) )
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i=3, n1-3
                        tmp2(i, j,k)=     B(0)*f(i,j,k) + B(2)*( f(i+2, j,k)+f(i-2, j,k))   &
                        +   B(1)*( f(i+1, j,k)+f(i-1, j,k) )
                    enddo
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    tmp2(n1-1, j,k)=   B(0)*f(n1-1, j,k) + B(1)*( f(1, j,k)+f(n1-2, j,k) )  &
                    +   B(2)*( f(2, j,k)+f(n1-3, j,k) )
                    tmp2(n1-2, j,k)=   B(0)*f(n1-2, j,k) + B(1)*( f(n1-1, j,k)+f(n1-3, j,k) )  &
                    +   B(2)*( f(1, j,k)+f(n1-4, j,k) )
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i = 2,n1-1
                        tmp2(i, j,k)=tmp2(i, j,k) - m(i)*tmp2(i-1, j,k)
                    enddo
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    tmp2(n1-1, j,k) = tmp2(n1-1, j,k)*bp(n1-1)
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i = n1-2, 1, -1
                        tmp2(i, j,k)=(tmp2(i, j,k) - A3(i)*tmp2(i+1, j,k))*bp(i)
                    enddo
                enddo
            enddo

            do k=1, n3e
                do j=1, n2e
                    sx(j,k) = (v(1)*tmp2(1,j,k)+v(n1-1)*tmp2(n1-1,j,k))/(1+v(1)*q(1)+v(n1-1)*q(n1-1))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i=1, n1-1
                        d2f(i, j,k)=tmp2(i, j,k)-sx(j,k)*q(i) + d2f(i, j,k)
                    enddo
                enddo
            enddo


        end if


        if ((bc_type.eq.Dirichlet)) then

            if (shifted) then

                do k=1,n3e
                    do j=1,n2e

                        tmp2(2, j,k)= Bl(-1)*f(1, j,k) + Bl(0)*f(2, j,k) + Bl(1)*f(3, j,k) + Bl(2)*f(4, j,k)

                        do i=3, n1-3
                            tmp2(i, j,k)=     B(0)*f(i, j,k) + B(2)*( f(i+2, j,k)+f(i-2, j,k))   &
                            +   B(1)*( f(i+1, j,k)+f(i-1, j,k) )
                        enddo

                        tmp2(n1-2, j,k)= Bl(-1)*f(n1-1, j,k) + Bl(0)*f(n1-2, j,k)     &
                        + Bl(1)*f(n1-3, j,k)+ Bl(2)*f(n1-4, j,k)


                        do i=3, n1-3
                            tmp2(i, j,k)= tmp2(i, j,k) - m(i)*tmp2(i-1, j,k)
                        enddo

                        tmp2(n1-2, j,k)= (tmp2(n1-2, j,k) - m(n1-2)*tmp2(n1-3, j,k))*bp(n1-2)


                        do i = n1-3, 2, -1
                            tmp2(i, j,k) = (tmp2(i, j,k) - A3(i)*tmp2(i+1, j,k))*bp(i)
                        enddo

                        do i=2, n1-2
                            d2f(i, j,k)= tmp2(i, j,k) + d2f(i, j,k)
                        enddo

                    enddo
                enddo

            else

                do k=1,n3e
                    do j=1,n2e

                        tmp2(2, j,k)= Bl(-1)*f(1, j,k) + Bl(0)*f(2, j,k) + Bl(1)*f(3, j,k) + Bl(2)*f(4, j,k)

                        do i=3, n1-2
                            tmp2(i, j,k)=     B(0)*f(i, j,k) + B(2)*( f(i+2, j,k)+f(i-2, j,k))   &
                            +   B(1)*( f(i+1, j,k)+f(i-1, j,k) )
                        enddo

                        tmp2(n1-1, j,k)= Bl(-1)*f(n1, j,k) + Bl(0)*f(n1-1, j,k)     &
                        + Bl(1)*f(n1-2, j,k)+ Bl(2)*f(n1-3, j,k)


                        do i=3, n1-2
                            tmp2(i, j,k)= tmp2(i, j,k) - m(i)*tmp2(i-1, j,k)
                        enddo

                        tmp2(n1-1, j,k)= (tmp2(n1-1, j,k) - m(n1-1)*tmp2(n1-2, j,k))*bp(n1-1)


                        do i = n1-2, 2, -1
                            tmp2(i, j,k) = (tmp2(i, j,k) - A3(i)*tmp2(i+1, j,k))*bp(i)
                        enddo

                        do i=2, n1-1
                            d2f(i, j,k)= tmp2(i, j,k) + d2f(i, j,k)
                        enddo

                    enddo
                enddo

            end if
        end if

        return

    end subroutine D2c_CPT_ACC_3Dx



    subroutine D2c_CPT_3Dy(f, d2f, n1,n1e,n2,n3,n3e, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: d2f

        integer :: i, j, k
        real(kind=8) rhs(n2)


        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h



        update_matrix=(last_n/=n2).or.(last_bc_type/=bc_type).or.(last_shifted.neqv.shifted).or.(last_h/=h)


        !write(*,*)'D1_OUCS3'

        ! MATRIX A, B definition -----------------------------------------------------------
         call fill_matrix(n2, shifted, bc_type, h)

        last_n=n2
        last_bc_type=bc_type
        last_shifted=shifted
        last_h=h

        ! SCHEME CALCULATION ---------------------------------------------------------------

        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    d2f(i,1,k)= B(0)*f(i,1,k) + B(1)*( f(i,2,k)+f(i,n2-1,k) )  &
                    +   B(2)*( f(i,3,k)+f(i,n2-2,k) )

                    d2f(i,2,k)= B(0)*f(i,2,k) + B(1)*( f(i,3,k)+f(i,1,k) )  &
                    +   B(2)*( f(i,4,k)+f(i,n2-1,k) )

                    do j=3, n2-3
                        d2f(i,j,k)= B(0)*f(i,j,k) + B(2)*( f(i,j+2,k)+f(i,j-2,k))   &
                        +   B(1)*( f(i,j+1,k)+f(i,j-1,k) )
                    enddo

                    d2f(i,n2-1,k)= B(0)*f(i,n2-1,k) + B(1)*( f(i,1,k)+f(i,n2-2,k) )  &
                    +   B(2)*( f(i,2,k)+f(i,n2-3,k) )


                    d2f(i,n2-2,k)= B(0)*f(i,n2-2,k) + B(1)*( f(i,n2-1,k)+f(i,n2-3,k) )  &
                    +   B(2)*( f(i,1,k)+f(i,n2-4,k) )


                    ! INVERSION
                    do j = 2,n2-1
                        d2f(i,j,k)=d2f(i,j,k) - m(j)*d2f(i,j-1,k)
                    enddo
                    d2f(i,n2-1,k) = d2f(i,n2-1,k)*bp(n2-1)
                    do j = n2-2, 1, -1
                        d2f(i,j,k)=(d2f(i,j,k) - A3(j)*d2f(i,j+1,k))*bp(j)
                    enddo


                    vy=0.d0
                    do j=1, n2-1
                        vy=vy+v(j)*d2f(i,j,k)
                    enddo
                    coef_q = vy/(1+vq)

                    do j=1, n2-1
                        d2f(i,j,k)=d2f(i,j,k)-coef_q*q(j)
                    enddo

                enddo
            enddo

        end if


        if ((bc_type.eq.Dirichlet)) then

            if (shifted) then

                do k=1, n3e
                    do i=1, n1e

                        d2f(i,2,k)= Bl(-1)*f(i,1,k) + Bl(0)*f(i,2,k) + Bl(1)*f(i,3,k) + Bl(2)*f(i,4,k)

                        do j=3, n2-3
                            d2f(i,j,k)= B(0)*f(i,j,k) + B(2)*( f(i,j+2,k)+f(i,j-2,k))   &
                            +   B(1)*( f(i,j+1,k)+f(i,j-1,k) )
                        enddo

                        d2f(i,n2-2,k)= Bl(-1)*f(i,n2-1,k) + Bl(0)*f(i,n2-2,k)     &
                        + Bl(1)*f(i,n2-3,k)+ Bl(2)*f(i,n2-4,k)

                        do j=3, n2-3
                            d2f(i,j,k)= d2f(i,j,k) - m(j)*d2f(i,j-1,k)
                        enddo

                        d2f(i,n2-2,k)= (d2f(i,n2-2,k) - m(n2-2)*d2f(i,n2-3,k))*bp(n2-2)


                        do j = n2-3, 2, -1
                            d2f(i,j,k) = (d2f(i,j,k) - A3(j)*d2f(i,j+1,k))*bp(j)
                        enddo

                    enddo
                enddo

            else

                do k=1, n3e
                    do i=1, n1e

                        d2f(i,2,k)= Bl(-1)*f(i,1,k) + Bl(0)*f(i,2,k) + Bl(1)*f(i,3,k) + Bl(2)*f(i,4,k)

                        do j=3, n2-2
                            d2f(i,j,k)= B(0)*f(i,j,k) + B(2)*( f(i,j+2,k)+f(i,j-2,k))   &
                            +   B(1)*( f(i,j+1,k)+f(i,j-1,k) )
                        enddo

                        d2f(i,n2-1,k)= Bl(-1)*f(i,n2,k) + Bl(0)*f(i,n2-1,k)     &
                        + Bl(1)*f(i,n2-2,k)+ Bl(2)*f(i,n2-3,k)


                        do j=3, n2-2
                            d2f(i,j,k)= d2f(i,j,k) - m(j)*d2f(i,j-1,k)
                        enddo

                        d2f(i,n2-1,k)= (d2f(i,n2-1,k) - m(n2-1)*d2f(i,n2-2,k))*bp(n2-1)


                        do j = n2-2, 2, -1
                            d2f(i,j,k) = (d2f(i,j,k) - A3(j)*d2f(i,j+1,k))*bp(j)
                        enddo

                    enddo
                enddo

            end if

        end if

        return

    end subroutine D2c_CPT_3Dy

    subroutine D2c_CPT_MULT_3Dy(f, d2f, n1,n1e,n2,n3,n3e, h, shifted, bc_type, g)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(:), intent(in)              :: g
        real(kind=8), dimension(n1,n2,n3), intent(in)       :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)      :: d2f

        integer :: i, j, k
        real(kind=8) rhs(n2)


        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h



        update_matrix=(last_n/=n2).or.(last_bc_type/=bc_type).or.(last_shifted.neqv.shifted).or.(last_h/=h)


        !write(*,*)'D1_OUCS3'

        ! MATRIX A, B definition -----------------------------------------------------------
         call fill_matrix(n2, shifted, bc_type, h)

        last_n=n2
        last_bc_type=bc_type
        last_shifted=shifted
        last_h=h

        ! SCHEME CALCULATION ---------------------------------------------------------------

        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    d2f(i,1,k)= B(0)*f(i,1,k) + B(1)*( f(i,2,k)+f(i,n2-1,k) )  &
                    +   B(2)*( f(i,3,k)+f(i,n2-2,k) )

                    d2f(i,2,k)= B(0)*f(i,2,k) + B(1)*( f(i,3,k)+f(i,1,k) )  &
                    +   B(2)*( f(i,4,k)+f(i,n2-1,k) )

                    do j=3, n2-3
                        d2f(i,j,k)= B(0)*f(i,j,k) + B(2)*( f(i,j+2,k)+f(i,j-2,k))   &
                        +   B(1)*( f(i,j+1,k)+f(i,j-1,k) )
                    enddo

                    d2f(i,n2-1,k)= B(0)*f(i,n2-1,k) + B(1)*( f(i,1,k)+f(i,n2-2,k) )  &
                    +   B(2)*( f(i,2,k)+f(i,n2-3,k) )


                    d2f(i,n2-2,k)= B(0)*f(i,n2-2,k) + B(1)*( f(i,n2-1,k)+f(i,n2-3,k) )  &
                    +   B(2)*( f(i,1,k)+f(i,n2-4,k) )


                    ! INVERSION
                    do j = 2,n2-1
                        d2f(i,j,k)=d2f(i,j,k) - m(j)*d2f(i,j-1,k)
                    enddo
                    d2f(i,n2-1,k) = d2f(i,n2-1,k)*bp(n2-1)
                    do j = n2-2, 1, -1
                        d2f(i,j,k)=(d2f(i,j,k) - A3(j)*d2f(i,j+1,k))*bp(j)
                    enddo


                    vy=0.d0
                    do j=1, n2-1
                        vy=vy+v(j)*d2f(i,j,k)
                    enddo
                    coef_q = vy/(1+vq)

                    do j=1, n2-1
                        d2f(i,j,k)=(d2f(i,j,k)-coef_q*q(j))*g(j)
                    enddo

                enddo
            enddo

        end if


        if ((bc_type.eq.Dirichlet)) then

            if (shifted) then

                do k=1, n3e
                    do i=1, n1e

                        d2f(i,2,k)= Bl(-1)*f(i,1,k) + Bl(0)*f(i,2,k) + Bl(1)*f(i,3,k) + Bl(2)*f(i,4,k)

                        do j=3, n2-3
                            d2f(i,j,k)= B(0)*f(i,j,k) + B(2)*( f(i,j+2,k)+f(i,j-2,k))   &
                            +   B(1)*( f(i,j+1,k)+f(i,j-1,k) )
                        enddo

                        d2f(i,n2-2,k)= Bl(-1)*f(i,n2-1,k) + Bl(0)*f(i,n2-2,k)     &
                        + Bl(1)*f(i,n2-3,k)+ Bl(2)*f(i,n2-4,k)

                        do j=3, n2-3
                            d2f(i,j,k)= d2f(i,j,k) - m(j)*d2f(i,j-1,k)
                        enddo

                        d2f(i,n2-2,k)= (d2f(i,n2-2,k) - m(n2-2)*d2f(i,n2-3,k))*bp(n2-2)


                        do j = n2-3, 2, -1
                            d2f(i,j,k) = (d2f(i,j,k) - A3(j)*d2f(i,j+1,k))*bp(j)
                        enddo

                        do j=2, n2-2
                            d2f(i,j,k)= d2f(i,j,k) * g(j)
                        enddo

                    enddo
                enddo

            else

                do k=1, n3e
                    do i=1, n1e

                        d2f(i,2,k)= Bl(-1)*f(i,1,k) + Bl(0)*f(i,2,k) + Bl(1)*f(i,3,k) + Bl(2)*f(i,4,k)

                        do j=3, n2-2
                            d2f(i,j,k)= B(0)*f(i,j,k) + B(2)*( f(i,j+2,k)+f(i,j-2,k))   &
                            +   B(1)*( f(i,j+1,k)+f(i,j-1,k) )
                        enddo

                        d2f(i,n2-1,k)= Bl(-1)*f(i,n2,k) + Bl(0)*f(i,n2-1,k)     &
                        + Bl(1)*f(i,n2-2,k)+ Bl(2)*f(i,n2-3,k)


                        do j=3, n2-2
                            d2f(i,j,k)= d2f(i,j,k) - m(j)*d2f(i,j-1,k)
                        enddo

                        d2f(i,n2-1,k)= (d2f(i,n2-1,k) - m(n2-1)*d2f(i,n2-2,k))*bp(n2-1)


                        do j = n2-2, 2, -1
                            d2f(i,j,k) = (d2f(i,j,k) - A3(j)*d2f(i,j+1,k))*bp(j)
                        enddo

                        do j=2, n2-1
                            d2f(i,j,k)= d2f(i,j,k) * g(j)
                        enddo

                    enddo
                enddo

            end if

        end if

        return

    end subroutine D2c_CPT_MULT_3Dy

    subroutine D2c_CPT_MULTACC_3Dy(f, d2f, n1,n1e,n2,n3,n3e, h, shifted, bc_type, g)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(:), intent(in)              :: g
        real(kind=8), dimension(n1,n2,n3), intent(in)       :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)      :: d2f

        real(kind=8), dimension(n2)             ::  tmp
        real(kind=8), dimension(n1,n2,n3)       :: tmp2

        integer :: i, j, k
        real(kind=8) rhs(n2)


        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n2, shifted, bc_type, h)

        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    tmp(1:n2-1)=d2f(i,1:n2-1,k)

                    d2f(i,1,k)= B(0)*f(i,1,k) + B(1)*( f(i,2,k)+f(i,n2-1,k) )  &
                    +   B(2)*( f(i,3,k)+f(i,n2-2,k) )

                    d2f(i,2,k)= B(0)*f(i,2,k) + B(1)*( f(i,3,k)+f(i,1,k) )  &
                    +   B(2)*( f(i,4,k)+f(i,n2-1,k) )

                    do j=3, n2-3
                        d2f(i,j,k)= B(0)*f(i,j,k) + B(2)*( f(i,j+2,k)+f(i,j-2,k))   &
                        +   B(1)*( f(i,j+1,k)+f(i,j-1,k) )
                    enddo

                    d2f(i,n2-1,k)= B(0)*f(i,n2-1,k) + B(1)*( f(i,1,k)+f(i,n2-2,k) )  &
                    +   B(2)*( f(i,2,k)+f(i,n2-3,k) )


                    d2f(i,n2-2,k)= B(0)*f(i,n2-2,k) + B(1)*( f(i,n2-1,k)+f(i,n2-3,k) )  &
                    +   B(2)*( f(i,1,k)+f(i,n2-4,k) )


                    ! INVERSION
                    do j = 2,n2-1
                        d2f(i,j,k)=d2f(i,j,k) - m(j)*d2f(i,j-1,k)
                    enddo
                    d2f(i,n2-1,k) = d2f(i,n2-1,k)*bp(n2-1)
                    do j = n2-2, 1, -1
                        d2f(i,j,k)=(d2f(i,j,k) - A3(j)*d2f(i,j+1,k))*bp(j)
                    enddo


                    vy=0.d0
                    do j=1, n2-1
                        vy=vy+v(j)*d2f(i,j,k)
                    enddo
                    coef_q = vy/(1+vq)

                    do j=1, n2-1
                        d2f(i,j,k)=(d2f(i,j,k)-coef_q*q(j))*g(j) + tmp(j)
                    enddo

                enddo
            enddo

        end if


        if ((bc_type.eq.Dirichlet)) then

            if (shifted) then

                do k=1, n3e
                    do i=1, n1e

                        tmp2(i,2,k)= Bl(-1)*f(i,1,k) + Bl(0)*f(i,2,k) + Bl(1)*f(i,3,k) + Bl(2)*f(i,4,k)

                        do j=3, n2-3
                            tmp2(i,j,k)= B(0)*f(i,j,k) + B(2)*( f(i,j+2,k)+f(i,j-2,k))   &
                            +   B(1)*( f(i,j+1,k)+f(i,j-1,k) )
                        enddo

                        tmp2(i,n2-2,k)= Bl(-1)*f(i,n2-1,k) + Bl(0)*f(i,n2-2,k)     &
                        + Bl(1)*f(i,n2-3,k)+ Bl(2)*f(i,n2-4,k)

                        do j=3, n2-3
                            tmp2(i,j,k)= tmp2(i,j,k) - m(j)*tmp2(i,j-1,k)
                        enddo

                        tmp2(i,n2-2,k)= (tmp2(i,n2-2,k) - m(n2-2)*tmp2(i,n2-3,k))*bp(n2-2)


                        do j = n2-3, 2, -1
                            tmp2(i,j,k) = (tmp2(i,j,k) - A3(j)*tmp2(i,j+1,k))*bp(j)
                        enddo

                        do j=2, n2-2
                            d2f(i,j,k)= tmp2(i,j,k) * g(j) + d2f(i,j,k)
                        enddo

                    enddo
                enddo

            else

                do k=1, n3e
                    do i=1, n1e

                        tmp2(i,2,k)= Bl(-1)*f(i,1,k) + Bl(0)*f(i,2,k) + Bl(1)*f(i,3,k) + Bl(2)*f(i,4,k)

                        do j=3, n2-2
                            tmp2(i,j,k)= B(0)*f(i,j,k) + B(2)*( f(i,j+2,k)+f(i,j-2,k))   &
                            +   B(1)*( f(i,j+1,k)+f(i,j-1,k) )
                        enddo

                        tmp2(i,n2-1,k)= Bl(-1)*f(i,n2,k) + Bl(0)*f(i,n2-1,k)     &
                        + Bl(1)*f(i,n2-2,k)+ Bl(2)*f(i,n2-3,k)


                        do j=3, n2-2
                            tmp2(i,j,k)= tmp2(i,j,k) - m(j)*tmp2(i,j-1,k)
                        enddo

                        tmp2(i,n2-1,k)= (tmp2(i,n2-1,k) - m(n2-1)*tmp2(i,n2-2,k))*bp(n2-1)


                        do j = n2-2, 2, -1
                            tmp2(i,j,k) = (tmp2(i,j,k) - A3(j)*tmp2(i,j+1,k))*bp(j)
                        enddo

                        do j=2, n2-1
                            d2f(i,j,k)= tmp2(i,j,k) * g(j) + d2f(i,j,k)
                        enddo

                    enddo
                enddo

            end if

        end if

        return

    end subroutine D2c_CPT_MULTACC_3Dy

    subroutine D2c_CPT_3Dz(f, d2f, n1,n1e,n2,n2e,n3, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: d2f

        real(kind=8), dimension(n1,n2) :: sz

        integer :: i,j,k


        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n3, shifted, bc_type, h)

        ! SCHEME CALCULATION ---------------------------------------------------------------

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e
                    d2f(i,j,1)= B(0)*f(i,j,1) + B(1)*( f(i,j,2)+f(i,j,n3-1) )  &
                    +   B(2)*( f(i,j,3)+f(i,j,n3-2) )
                    d2f(i,j,2)= B(0)*f(i,j,2) + B(1)*( f(i,j,3)+f(i,j,1) )  &
                    +   B(2)*( f(i,j,4)+f(i,j,n3-1) )
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k=3, n3-3
                        d2f(i,j,k)= B(0)*f(i,j,k) + B(2)*( f(i,j,k+2)+f(i,j,k-2))   &
                        +   B(1)*( f(i,j,k+1)+f(i,j,k-1) )
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    d2f(i,j,n3-1)= B(0)*f(i,j,n3-1) + B(1)*( f(i,j,1)+f(i,j,n3-2) )  &
                    +   B(2)*( f(i,j,2)+f(i,j,n3-3) )
                    d2f(i,j,n3-2)= B(0)*f(i,j,n3-2) + B(1)*( f(i,j,n3-1)+f(i,j,n3-3) )  &
                    +   B(2)*( f(i,j,1)+f(i,j,n3-4) )
                enddo
            enddo


            ! INVERSION
            do j=1, n2e
                do i=1, n1e
                    do k = 2,n3-1
                        d2f(i,j,k)=d2f(i,j,k) - m(k)*d2f(i,j,k-1)
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    d2f(i,j,n3-1) = d2f(i,j,n3-1)*bp(n3-1)
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k = n3-2, 1, -1
                        d2f(i,j,k)=(d2f(i,j,k) - A3(k)*d2f(i,j,k+1))*bp(k)
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    sz(i,j) = (v(1)*d2f(i,j,1)+v(n3-1)*d2f(i,j,n3-1))/(1+v(1)*q(1)+v(n3-1)*q(n3-1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k=1, n3-1
                        d2f(i,j,k)=d2f(i,j,k)-sz(i,j)*q(k)
                    enddo
                enddo
            enddo

        end if


        if ((bc_type.eq.Dirichlet)) then

            if (shifted) then

                do j=1, n2e
                    do i=1, n1e

                        d2f(i,j,2)= Bl(-1)*f(i,j,1) + Bl(0)*f(i,j,2) + Bl(1)*f(i,j,3) + Bl(2)*f(i,j,4)

                        do k=3, n3-3
                            d2f(i,j,k)= B(0)*f(i,j,k) + B(2)*( f(i,j,k+2)+f(i,j,k-2))   &
                            +   B(1)*( f(i,j,k+1)+f(i,j,k-1) )
                        enddo

                        d2f(i,j,n3-2)= Bl(-1)*f(i,j,n3-1) + Bl(0)*f(i,j,n3-2)     &
                        + Bl(1)*f(i,j,n3-3)+ Bl(2)*f(i,j,n3-4)


                        do k=3, n3-3
                            d2f(i,j,k)= d2f(i,j,k) - m(k)*d2f(i,j,k-1)
                        enddo

                        d2f(i,j,n3-2)= (d2f(i,j,n3-2) - m(n3-2)*d2f(i,j,n3-3))*bp(n3-2)


                        do k = n3-3, 2, -1
                            d2f(i,j,k) = (d2f(i,j,k) - A3(k)*d2f(i,j,k+1))*bp(k)
                        enddo

                    enddo
                enddo

            else

                do j=1, n2e
                    do i=1, n1e

                        d2f(i,j,2)= Bl(-1)*f(i,j,1) + Bl(0)*f(i,j,2) + Bl(1)*f(i,j,3) + Bl(2)*f(i,j,4)

                        do k=3, n3-2
                            d2f(i,j,k)= B(0)*f(i,j,k) + B(2)*( f(i,j,k+2)+f(i,j,k-2))   &
                            +   B(1)*( f(i,j,k+1)+f(i,j,k-1) )
                        enddo

                        d2f(i,j,n3-1)= Bl(-1)*f(i,j,n3) + Bl(0)*f(i,j,n3-1)     &
                        + Bl(1)*f(i,j,n3-2)+ Bl(2)*f(i,j,n3-3)


                        do k=3, n3-2
                            d2f(i,j,k)= d2f(i,j,k) - m(k)*d2f(i,j,k-1)
                        enddo

                        d2f(i,j,n3-1)= (d2f(i,j,n3-1) - m(n3-1)*d2f(i,j,n3-2))*bp(n3-1)


                        do k = n3-2, 2, -1
                            d2f(i,j,k) = (d2f(i,j,k) - A3(k)*d2f(i,j,k+1))*bp(k)
                        enddo

                    enddo
                enddo

            end if
        end if

        return

    end subroutine D2c_CPT_3Dz

end module CPT_D2c




