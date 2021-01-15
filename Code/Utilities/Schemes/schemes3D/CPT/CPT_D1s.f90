module CPT_D1s
    use boundaries_types
    implicit none

    private

    real(kind=8), allocatable, dimension(:) :: q,v
    real(kind=8), allocatable, dimension(:) :: m, bp, cp
    real(kind=8), dimension(:), allocatable :: A1, A2, A3
    real(kind=8)                                 :: vy, vq, coef_q
    real(kind=8)    :: B(2), Bl(-1:1)

    public::D1s_CPT_3Dx,D1s_CPT_ACC_3Dx, D1s_CPT_3Dy,D1s_CPT_MULT_3Dy,D1s_CPT_MULTACC_3Dy, D1s_CPT_3Dz,D1s_CPT_ACC_3Dz

contains

    subroutine fill_matrix(n_max, shifted, bc_type, h)
        use CPT_INIT
        implicit none
        integer     :: n_max
        logical     :: shifted
        integer     :: bc_type
        real*8      :: h
        integer     :: i
        real*8      :: al

        if (allocated(A1)) then
            deallocate(A1)
            deallocate(A2)
            deallocate(A3)
        end if

        allocate(A1(n_max))
        allocate(A2(n_max))
        allocate(A3(n_max))


        al=9.d0/62.d0
        B(1)=(3*(3-2*al))/(8*h)
        B(2)= (22*al-1)/(24*h)


        do i= 1, n_max
            A1(i) = al
            A2(i) = 1.d0
            A3(i) = al
        enddo

        select case (bc_type)

            case (periodic)

                if (allocated(m)) deallocate(m)
                if (allocated(bp)) deallocate(bp)
                if (allocated(cp)) deallocate(cp)
                if (allocated(q)) deallocate(q)
                if (allocated(v)) deallocate(v)

                allocate(q(n_max-1))
                allocate(v(n_max-1))
                allocate(m(n_max-1), bp(n_max-1), cp(n_max-1))

                call TS_initPR(A1(1:n_max-1), A2(1:n_max-1), A3(1:n_max-1), bp, cp, m, q,v,vq, n_max-1)

            case (Dirichlet)

                Bl(-1)=-1.d0               / h
                Bl(0)= 2.d0                / h
                Bl(1)=-1.d0                / h

                A1(n_max-1)=-1.d0
                A2(n_max-1)=1.d0

                A2(1)=1.d0
                A3(1)=-1.d0

                if (allocated(m)) deallocate(m)
                if (allocated(bp)) deallocate(bp)

                allocate(m(n_max-1), bp(n_max-1))

                call TS_initDIR(A1(1:n_max-1), A2(1:n_max-1), A3(1:n_max-1), bp, m, n_max-1)

            case (antisymetric)

                A1(n_max-1)=al
                A2(n_max-1)=1.d0+al

                A2(1)=1.d0+al
                A3(1)=al

                if (allocated(m)) deallocate(m)
                if (allocated(bp)) deallocate(bp)

                allocate(m(n_max-1), bp(n_max-1))

                call TS_initDIR(A1(1:n_max-1), A2(1:n_max-1), A3(1:n_max-1), bp, m, n_max-1)

            case default

        end select

    end subroutine fill_matrix

    subroutine D1s_CPT_3Dx(f, d1f, n1,n2,n2e,n3,n3e, h, shifted, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: d1f

        real(kind=8), dimension(n2,n3)                      :: sx

        integer :: i, j, k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n1, shifted, bc_type, h)

        if (bc_type==periodic) then

            do k=1,n3e
                do j=1,n2e
                    d1f(1,j,k)= B(2)*(f(3,j,k)-f(n1-1,j,k)) + B(1)*(f(2,j,k) - f(1,j,k))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i=2,n1-3
                        d1f(i,j,k)= B(2)*(f(i+2,j,k)-f(i-1,j,k)) + B(1)*(f(i+1,j,k) - f(i,j,k))
                    enddo
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    d1f(n1-2,j,k)= B(2)*(f(1,j,k)-f(n1-3,j,k)) + B(1)*(f(n1-1,j,k) - f(n1-2,j,k))
                    d1f(n1-1,j,k)= B(2)*(f(2,j,k)-f(n1-2,j,k)) + B(1)*(f(1,j,k) - f(n1-1,j,k))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i = 2,n1-1
                        d1f(i, j,k)=d1f(i, j,k) - m(i)*d1f(i-1, j,k)
                    enddo
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    d1f(n1-1, j,k) = d1f(n1-1, j,k)*bp(n1-1)
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i = n1-2, 1, -1
                        d1f(i, j,k)=(d1f(i, j,k) - A3(i)*d1f(i+1, j,k))*bp(i)
                    enddo
                enddo
            enddo

            do k=1, n3e
                do j=1, n2e
                    sx(j,k) = (v(1)*d1f(1,j,k)+v(n1-1)*d1f(n1-1,j,k))/(1+v(1)*q(1)+v(n1-1)*q(n1-1))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i=1, n1-1
                        d1f(i, j,k)=d1f(i, j,k)-sx(j,k)*q(i)
                    enddo
                enddo
            enddo

        end if

        if (bc_type==Dirichlet) then

            do k=1,n3e
                do j=1,n2e

                    do i=2,n1-2
                        d1f(i,j,k)= B(2)*(f(i+2,j,k)-f(i-1,j,k)) + B(1)*(f(i+1,j,k) - f(i,j,k))
                    enddo

                    d1f(1,j,k)= ( Bl(-1)*f(1,j,k) + Bl(0) * f(2,j,k) + Bl(1) * f(3,j,k) )

                    d1f(n1-1,j,k)= -( + Bl(-1)*f(n1,j,k) + Bl(0) * f(n1-1,j,k) + Bl(1)*f(n1-2,j,k) )


                    do i=2, n1-2
                        d1f(i, j,k)= d1f(i, j,k) - m(i)*d1f(i-1, j,k)
                    enddo

                    d1f(n1-1, j,k)= (d1f(n1-1, j,k) - m(n1-1)*d1f(n1-2, j,k))*bp(n1-1)


                    do i = n1-2, 1, -1
                        d1f(i, j,k) = (d1f(i, j,k) - A3(i)*d1f(i+1, j,k))*bp(i)
                    enddo

                enddo
            enddo

        end if

        if (bc_type==antisymetric) then

            do k=1,n3e
                do j=1,n2e

                    do i=2,n1-2
                        d1f(i,j,k)= B(2)*(f(i+2,j,k)-f(i-1,j,k)) + B(1)*(f(i+1,j,k) - f(i,j,k))
                    enddo

                    d1f(1,j,k)= B(2)*(f(3,j,k)+f(2,j,k)) + B(1)*(f(2,j,k) - f(1,j,k))

                    d1f(n1-1,j,k)= B(2)*(-f(n1-1,j,k)-f(n1-2,j,k)) + B(1)*(f(n1,j,k) - f(n1-1,j,k))


                    do i=2, n1-2
                        d1f(i, j,k)= d1f(i, j,k) - m(i)*d1f(i-1, j,k)
                    enddo

                    d1f(n1-1, j,k)= (d1f(n1-1, j,k) - m(n1-1)*d1f(n1-2, j,k))*bp(n1-1)


                    do i = n1-2, 1, -1
                        d1f(i, j,k) = (d1f(i, j,k) - A3(i)*d1f(i+1, j,k))*bp(i)
                    enddo

                enddo
            enddo

        end if


        return

    end subroutine D1s_CPT_3Dx

    subroutine D1s_CPT_ACC_3Dx(f, d1f, n1,n2,n2e,n3,n3e, h, shifted, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: d1f

        real*8, dimension(n1)       :: tmp
        real*8, dimension(n1,n2,n3) :: tmp2

        real(kind=8), dimension(n2,n3)                      :: sx

        integer :: i, j, k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n1, shifted, bc_type, h)

        if (bc_type==periodic) then

            do k=1,n3e
                do j=1,n2e
                    tmp2(1,j,k)= B(2)*(f(3,j,k)-f(n1-1,j,k)) + B(1)*(f(2,j,k) - f(1,j,k))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i=2,n1-3
                        tmp2(i,j,k)= B(2)*(f(i+2,j,k)-f(i-1,j,k)) + B(1)*(f(i+1,j,k) - f(i,j,k))
                    enddo
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    tmp2(n1-2,j,k)= B(2)*(f(1,j,k)-f(n1-3,j,k)) + B(1)*(f(n1-1,j,k) - f(n1-2,j,k))
                    tmp2(n1-1,j,k)= B(2)*(f(2,j,k)-f(n1-2,j,k)) + B(1)*(f(1,j,k) - f(n1-1,j,k))
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
                        d1f(i, j,k)=d1f(i, j,k)+tmp2(i, j,k)-sx(j,k)*q(i)
                    enddo
                enddo
            enddo

        end if

        if (bc_type==Dirichlet) then

            do k=1,n3e
                do j=1,n2e

                    tmp(1:n1-1)=d1f(1:n1-1,j,k)

                    do i=2,n1-2
                        d1f(i,j,k)= B(2)*(f(i+2,j,k)-f(i-1,j,k)) + B(1)*(f(i+1,j,k) - f(i,j,k))
                    enddo

                    d1f(1,j,k)= ( Bl(-1)*f(1,j,k) + Bl(0) * f(2,j,k) + Bl(1) * f(3,j,k) )

                    d1f(n1-1,j,k)= -( + Bl(-1)*f(n1,j,k) + Bl(0) * f(n1-1,j,k) + Bl(1)*f(n1-2,j,k) )


                    do i=2, n1-2
                        d1f(i, j,k)= d1f(i, j,k) - m(i)*d1f(i-1, j,k)
                    enddo

                    d1f(n1-1, j,k)= (d1f(n1-1, j,k) - m(n1-1)*d1f(n1-2, j,k))*bp(n1-1)


                    do i = n1-2, 1, -1
                        d1f(i, j,k) = (d1f(i, j,k) - A3(i)*d1f(i+1, j,k))*bp(i)
                    enddo
                    do i=1, n1-1
                        d1f(i, j,k)= tmp(i)+d1f(i, j,k)
                    enddo

                enddo
            enddo

        end if


        return

    end subroutine D1s_CPT_ACC_3Dx



    subroutine D1s_CPT_3Dy(f, d1f, n1,n1e,n2,n3,n3e, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: d1f

        integer :: i,j,k

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


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=2,n2-2
                        d1f(i,j,k)= B(2)*(f(i,j+2,k)-f(i,j-1,k)) + B(1)*(f(i,j+1,k) - f(i,j,k))
                    enddo

                    d1f(i,1,k)= B(2)*(f(i,3,k)-f(i,n2-1,k)) + B(1)*(f(i,2,k) - f(i,1,k))

                    d1f(i,n2-1,k)= B(2)*(f(i,2,k)-f(i,n2-2,k)) + B(1)*(f(i,1,k) - f(i,n2-1,k))

                    do j = 2,n2-1
                        d1f(i,j,k)=d1f(i,j,k) - m(j)*d1f(i,j-1,k)
                    enddo
                    d1f(i,n2-1,k) = d1f(i,n2-1,k)*bp(n2-1)
                    do j = n2-2, 1, -1
                        d1f(i,j,k)=(d1f(i,j,k) - A3(j)*d1f(i,j+1,k))*bp(j)
                    enddo

                    vy=0.d0
                    do j=1, n2-1
                        vy=vy+v(j)*d1f(i,j,k)
                    enddo
                    coef_q = vy/(1+vq)

                    do j=1, n2-1
                        d1f(i,j,k)=d1f(i,j,k)-coef_q*q(j)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e

                    do j=2,n2-2
                        d1f(i,j,k)= B(2)*(f(i,j+2,k)-f(i,j-1,k)) + B(1)*(f(i,j+1,k) - f(i,j,k))
                    enddo

                    d1f(i,1,k)= ( Bl(-1)*f(i,1,k) + Bl(0) * f(i,2,k) + Bl(1) * f(i,3,k) )

                    d1f(i,n2-1,k)= -( + Bl(-1)*f(i,n2,k) + Bl(0) * f(i,n2-1,k) + Bl(1)*f(i,n2-2,k) )


                    do j=2, n2-2
                        d1f(i,j,k)= d1f(i,j,k) - m(j)*d1f(i,j-1,k)
                    enddo

                    d1f(i,n2-1,k)= (d1f(i,n2-1,k) - m(n2-1)*d1f(i,n2-2,k))*bp(n2-1)


                    do j = n2-2, 1, -1
                        d1f(i,j,k) = (d1f(i,j,k) - A3(j)*d1f(i,j+1,k))*bp(j)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.antisymetric)) then

            do k=1, n3e
                do i=1, n1e

                    do j=2,n2-2
                        d1f(i,j,k)= B(2)*(f(i,j+2,k)-f(i,j-1,k)) + B(1)*(f(i,j+1,k) - f(i,j,k))
                    enddo

                    d1f(i,1,k)= ( B(2)*f(i,3,k) + B(1)*f(i,2,k) - B(1)*f(i,1,k) - B(2)*(2.d0*f(i,1,k)-f(i,2,k)) )

                    d1f(i,n2-1,k)= ( B(2)*(2.d0*f(i,n2,k)-f(i,n2-1,k)) + B(1)*f(i,n2,k) - B(1)*f(i,n2-1,k) - B(2)*f(i,n2-2,k) )


                    do j=2, n2-2
                        d1f(i,j,k)= d1f(i,j,k) - m(j)*d1f(i,j-1,k)
                    enddo

                    d1f(i,n2-1,k)= (d1f(i,n2-1,k) - m(n2-1)*d1f(i,n2-2,k))*bp(n2-1)


                    do j = n2-2, 1, -1
                        d1f(i,j,k) = (d1f(i,j,k) - A3(j)*d1f(i,j+1,k))*bp(j)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D1s_CPT_3Dy

    subroutine D1s_CPT_MULT_3Dy(f, d1f, n1,n1e,n2,n3,n3e, h, shifted, bc_type, g)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(:), intent(in)             :: g
        real(kind=8), dimension(n1,n2,n3), intent(in)       :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)      :: d1f

        real(kind=8), dimension(n1,n2,n3)       :: tmp2

        integer :: i,j,k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n2, shifted, bc_type, h)

        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=2,n2-2
                        d1f(i,j,k)= B(2)*(f(i,j+2,k)-f(i,j-1,k)) + B(1)*(f(i,j+1,k) - f(i,j,k))
                    enddo

                    d1f(i,1,k)= B(2)*(f(i,3,k)-f(i,n2-1,k)) + B(1)*(f(i,2,k) - f(i,1,k))

                    d1f(i,n2-1,k)= B(2)*(f(i,2,k)-f(i,n2-2,k)) + B(1)*(f(i,1,k) - f(i,n2-1,k))

                    do j = 2,n2-1
                        d1f(i,j,k)=d1f(i,j,k) - m(j)*d1f(i,j-1,k)
                    enddo
                    d1f(i,n2-1,k) = d1f(i,n2-1,k)*bp(n2-1)
                    do j = n2-2, 1, -1
                        d1f(i,j,k)=(d1f(i,j,k) - A3(j)*d1f(i,j+1,k))*bp(j)
                    enddo

                    vy=0.d0
                    do j=1, n2-1
                        vy=vy+v(j)*d1f(i,j,k)
                    enddo
                    coef_q = vy/(1+vq)

                    do j=1, n2-1
                        d1f(i,j,k)=(d1f(i,j,k)-coef_q*q(j))*g(j)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e

                    do j=2,n2-2
                        tmp2(i,j,k)= B(2)*(f(i,j+2,k)-f(i,j-1,k)) + B(1)*(f(i,j+1,k) - f(i,j,k))
                    enddo

                    tmp2(i,1,k)= ( Bl(-1)*f(i,1,k) + Bl(0) * f(i,2,k) + Bl(1) * f(i,3,k) )

                    tmp2(i,n2-1,k)= -( + Bl(-1)*f(i,n2,k) + Bl(0) * f(i,n2-1,k) + Bl(1)*f(i,n2-2,k) )


                    do j=2, n2-2
                        tmp2(i,j,k)= tmp2(i,j,k) - m(j)*tmp2(i,j-1,k)
                    enddo

                    tmp2(i,n2-1,k)= (tmp2(i,n2-1,k) - m(n2-1)*tmp2(i,n2-2,k))*bp(n2-1)


                    do j = n2-2, 1, -1
                        tmp2(i,j,k) = (tmp2(i,j,k) - A3(j)*tmp2(i,j+1,k))*bp(j)
                    enddo

                    do j=1, n2-1
                        d1f(i,j,k)= tmp2(i,j,k)*g(j)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.antisymetric)) then

            do k=1, n3e
                do i=1, n1e

                    do j=2,n2-2
                        tmp2(i,j,k)= B(2)*(f(i,j+2,k)-f(i,j-1,k)) + B(1)*(f(i,j+1,k) - f(i,j,k))
                    enddo

                    tmp2(i,1,k)= ( B(2)*f(i,3,k) + B(1)*f(i,2,k) - B(1)*f(i,1,k) - B(2)*(2.d0*f(i,1,k)-f(i,2,k)) )

                    tmp2(i,n2-1,k)= ( B(2)*(2.d0*f(i,n2,k)-f(i,n2-1,k)) + B(1)*f(i,n2,k) - B(1)*f(i,n2-1,k) - B(2)*f(i,n2-2,k) )


                    do j=2, n2-2
                        tmp2(i,j,k)= tmp2(i,j,k) - m(j)*tmp2(i,j-1,k)
                    enddo

                    tmp2(i,n2-1,k)= (tmp2(i,n2-1,k) - m(n2-1)*tmp2(i,n2-2,k))*bp(n2-1)


                    do j = n2-2, 1, -1
                        tmp2(i,j,k) = (tmp2(i,j,k) - A3(j)*tmp2(i,j+1,k))*bp(j)
                    enddo

                    do j=1, n2-1
                        d1f(i,j,k)= tmp2(i,j,k)*g(j)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D1s_CPT_MULT_3Dy

    subroutine D1s_CPT_MULTACC_3Dy(f, d1f, n1,n1e,n2,n3,n3e, h, shifted, bc_type, g)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(:), intent(in)             :: g
        real(kind=8), dimension(n1,n2,n3), intent(in)       :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)      :: d1f

        real(kind=8), dimension(n2)             :: tmp
        real(kind=8), dimension(n1,n2,n3)       :: tmp2

        integer :: i,j,k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n2, shifted, bc_type, h)

        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    tmp(1:n2-1)=d1f(i,1:n2-1,k)

                    do j=2,n2-2
                        d1f(i,j,k)= B(2)*(f(i,j+2,k)-f(i,j-1,k)) + B(1)*(f(i,j+1,k) - f(i,j,k))
                    enddo

                    d1f(i,1,k)= B(2)*(f(i,3,k)-f(i,n2-1,k)) + B(1)*(f(i,2,k) - f(i,1,k))

                    d1f(i,n2-1,k)= B(2)*(f(i,2,k)-f(i,n2-2,k)) + B(1)*(f(i,1,k) - f(i,n2-1,k))

                    do j = 2,n2-1
                        d1f(i,j,k)=d1f(i,j,k) - m(j)*d1f(i,j-1,k)
                    enddo
                    d1f(i,n2-1,k) = d1f(i,n2-1,k)*bp(n2-1)
                    do j = n2-2, 1, -1
                        d1f(i,j,k)=(d1f(i,j,k) - A3(j)*d1f(i,j+1,k))*bp(j)
                    enddo

                    vy=0.d0
                    do j=1, n2-1
                        vy=vy+v(j)*d1f(i,j,k)
                    enddo
                    coef_q = vy/(1+vq)

                    do j=1, n2-1
                        d1f(i,j,k)=(d1f(i,j,k)-coef_q*q(j))*g(j)+tmp(j)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e

                    do j=2,n2-2
                        tmp2(i,j,k)= B(2)*(f(i,j+2,k)-f(i,j-1,k)) + B(1)*(f(i,j+1,k) - f(i,j,k))
                    enddo

                    tmp2(i,1,k)= ( Bl(-1)*f(i,1,k) + Bl(0) * f(i,2,k) + Bl(1) * f(i,3,k) )

                    tmp2(i,n2-1,k)= -( + Bl(-1)*f(i,n2,k) + Bl(0) * f(i,n2-1,k) + Bl(1)*f(i,n2-2,k) )


                    do j=2, n2-2
                        tmp2(i,j,k)= tmp2(i,j,k) - m(j)*tmp2(i,j-1,k)
                    enddo

                    tmp2(i,n2-1,k)= (tmp2(i,n2-1,k) - m(n2-1)*tmp2(i,n2-2,k))*bp(n2-1)


                    do j = n2-2, 1, -1
                        tmp2(i,j,k) = (tmp2(i,j,k) - A3(j)*tmp2(i,j+1,k))*bp(j)
                    enddo

                    do j=1, n2-1
                        d1f(i,j,k)= tmp2(i,j,k)*g(j) + d1f(i,j,k)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.antisymetric)) then

            do k=1, n3e
                do i=1, n1e

                    do j=2,n2-2
                        tmp2(i,j,k)= B(2)*(f(i,j+2,k)-f(i,j-1,k)) + B(1)*(f(i,j+1,k) - f(i,j,k))
                    enddo

                    tmp2(i,1,k)= ( B(2)*f(i,3,k) + B(1)*f(i,2,k) - B(1)*f(i,1,k) - B(2)*(2.d0*f(i,1,k)-f(i,2,k)) )

                    tmp2(i,n2-1,k)= ( B(2)*(2.d0*f(i,n2,k)-f(i,n2-1,k)) + B(1)*f(i,n2,k) - B(1)*f(i,n2-1,k) - B(2)*f(i,n2-2,k) )


                    do j=2, n2-2
                        tmp2(i,j,k)= tmp2(i,j,k) - m(j)*tmp2(i,j-1,k)
                    enddo

                    tmp2(i,n2-1,k)= (tmp2(i,n2-1,k) - m(n2-1)*tmp2(i,n2-2,k))*bp(n2-1)


                    do j = n2-2, 1, -1
                        tmp2(i,j,k) = (tmp2(i,j,k) - A3(j)*tmp2(i,j+1,k))*bp(j)
                    enddo

                    do j=1, n2-1
                        d1f(i,j,k)= tmp2(i,j,k)*g(j) + d1f(i,j,k)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D1s_CPT_MULTACC_3Dy

    subroutine D1s_CPT_3Dz(f, d1f, n1,n1e,n2,n2e,n3, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: d1f

        real(kind=8), dimension(n1,n2)                      :: sz

        integer :: i,j,k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n3, shifted, bc_type, h)

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,1)= B(2)*(f(i,j,3)-f(i,j,n3-1)) + B(1)*(f(i,j,2) - f(i,j,1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k=2,n3-3
                        d1f(i,j,k)= B(2)*(f(i,j,k+2)-f(i,j,k-1)) + B(1)*(f(i,j,k+1) - f(i,j,k))
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,n3-2)= B(2)*(f(i,j,1)-f(i,j,n3-3)) + B(1)*(f(i,j,n3-1) - f(i,j,n3-2))
                    d1f(i,j,n3-1)= B(2)*(f(i,j,2)-f(i,j,n3-2)) + B(1)*(f(i,j,1) - f(i,j,n3-1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k = 2,n3-1
                        d1f(i,j,k)=d1f(i,j,k) - m(k)*d1f(i,j,k-1)
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,n3-1) = d1f(i,j,n3-1)*bp(n3-1)
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k = n3-2, 1, -1
                        d1f(i,j,k)=(d1f(i,j,k) - A3(k)*d1f(i,j,k+1))*bp(k)
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    sz(i,j) = (v(1)*d1f(i,j,1)+v(n3-1)*d1f(i,j,n3-1))/(1+v(1)*q(1)+v(n3-1)*q(n3-1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k=1, n3-1
                        d1f(i,j,k)=d1f(i,j,k)-sz(i,j)*q(k)
                    enddo
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e

                    do k=2,n3-2
                        d1f(i,j,k)= B(2)*(f(i,j,k+2)-f(i,j,k-1)) + B(1)*(f(i,j,k+1) - f(i,j,k))
                    enddo

                    d1f(i,j,1)= ( Bl(-1)*f(i,j,1) + Bl(0) * f(i,j,2) + Bl(1) * f(i,j,3) )

                    d1f(i,j,n3-1)= -( + Bl(-1)*f(i,j,n3) + Bl(0) * f(i,j,n3-1) + Bl(1)*f(i,j,n3-2) )


                    do k=2, n3-2
                        d1f(i,j,k)= d1f(i,j,k) - m(k)*d1f(i,j,k-1)
                    enddo

                    d1f(i,j,n3-1)= (d1f(i,j,n3-1) - m(n3-1)*d1f(i,j,n3-2))*bp(n3-1)


                    do k = n3-2, 1, -1
                        d1f(i,j,k) = (d1f(i,j,k) - A3(k)*d1f(i,j,k+1))*bp(k)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D1s_CPT_3Dz

    subroutine D1s_CPT_ACC_3Dz(f, d1f, n1,n1e,n2,n2e,n3, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: d1f

        real(kind=8), dimension(n3)         :: tmp
        real(kind=8), dimension(n1,n2,n3)   :: tmp2

        real(kind=8), dimension(n1,n2)                      :: sz

        integer :: i,j,k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n3, shifted, bc_type, h)

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e
                    tmp2(i,j,1)= B(2)*(f(i,j,3)-f(i,j,n3-1)) + B(1)*(f(i,j,2) - f(i,j,1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k=2,n3-3
                        tmp2(i,j,k)= B(2)*(f(i,j,k+2)-f(i,j,k-1)) + B(1)*(f(i,j,k+1) - f(i,j,k))
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    tmp2(i,j,n3-2)= B(2)*(f(i,j,1)-f(i,j,n3-3)) + B(1)*(f(i,j,n3-1) - f(i,j,n3-2))
                    tmp2(i,j,n3-1)= B(2)*(f(i,j,2)-f(i,j,n3-2)) + B(1)*(f(i,j,1) - f(i,j,n3-1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k = 2,n3-1
                        tmp2(i,j,k)=tmp2(i,j,k) - m(k)*tmp2(i,j,k-1)
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    tmp2(i,j,n3-1) = tmp2(i,j,n3-1)*bp(n3-1)
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k = n3-2, 1, -1
                        tmp2(i,j,k)=(tmp2(i,j,k) - A3(k)*tmp2(i,j,k+1))*bp(k)
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    sz(i,j) = (v(1)*tmp2(i,j,1)+v(n3-1)*tmp2(i,j,n3-1))/(1+v(1)*q(1)+v(n3-1)*q(n3-1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k=1, n3-1
                        d1f(i,j,k)=tmp2(i,j,k)-sz(i,j)*q(k) + d1f(i,j,k)
                    enddo
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e

                    tmp(1:n3-1)=d1f(i,j,1:n3-1)

                    do k=2,n3-2
                        d1f(i,j,k)= B(2)*(f(i,j,k+2)-f(i,j,k-1)) + B(1)*(f(i,j,k+1) - f(i,j,k))
                    enddo

                    d1f(i,j,1)= ( Bl(-1)*f(i,j,1) + Bl(0) * f(i,j,2) + Bl(1) * f(i,j,3) )

                    d1f(i,j,n3-1)= -( + Bl(-1)*f(i,j,n3) + Bl(0) * f(i,j,n3-1) + Bl(1)*f(i,j,n3-2) )


                    do k=2, n3-2
                        d1f(i,j,k)= d1f(i,j,k) - m(k)*d1f(i,j,k-1)
                    enddo

                    d1f(i,j,n3-1)= (d1f(i,j,n3-1) - m(n3-1)*d1f(i,j,n3-2))*bp(n3-1)


                    do k = n3-2, 1, -1
                        d1f(i,j,k) = (d1f(i,j,k) - A3(k)*d1f(i,j,k+1))*bp(k)
                    enddo

                    do k=1, n3-1
                        d1f(i,j,k)= d1f(i,j,k) +tmp(k)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D1s_CPT_ACC_3Dz

end module CPT_D1s

module CPT_D1s_sh
    use boundaries_types
    implicit none

    private

    real(kind=8), allocatable, dimension(:) :: q,v
    real(kind=8), allocatable, dimension(:) :: m, bp, cp
    real(kind=8), dimension(:), allocatable :: A1, A2, A3
    real(kind=8)                                 :: vy, vq, coef_q
    real(kind=8)    :: B(2), Bl(-1:1)

    public::D1ssh_CPT_3Dx,D1ssh_CPT_ACC_3Dx, D1ssh_CPT_3Dy,D1ssh_CPT_MULT_3Dy,D1ssh_CPT_MULTACC_3Dy, D1ssh_CPT_3Dz,D1ssh_CPT_ACC_3Dz

contains

    subroutine fill_matrix(n_max, shifted, bc_type, h)
        use CPT_INIT
        implicit none
        integer     :: n_max
        logical     :: shifted
        integer     :: bc_type
        real*8      :: h
        integer     :: i
        real*8      :: al

        if (allocated(A1)) then
            deallocate(A1)
            deallocate(A2)
            deallocate(A3)
        end if

        allocate(A1(n_max))
        allocate(A2(n_max))
        allocate(A3(n_max))


        al=9.d0/62.d0
        B(1)=(3*(3-2*al))/(8*h)
        B(2)= (22*al-1)/(24*h)


        do i= 1, n_max
            A1(i) = al
            A2(i) = 1.d0
            A3(i) = al
        enddo

        select case (bc_type)

            case (periodic)

                if (allocated(m)) deallocate(m)
                if (allocated(bp)) deallocate(bp)
                if (allocated(cp)) deallocate(cp)
                if (allocated(q)) deallocate(q)
                if (allocated(v)) deallocate(v)

                allocate(q(n_max-1))
                allocate(v(n_max-1))
                allocate(m(n_max-1), bp(n_max-1), cp(n_max-1))

                call TS_initPR(A1(1:n_max-1), A2(1:n_max-1), A3(1:n_max-1), bp, cp, m, q,v,vq, n_max-1)

            case (Dirichlet)

                Bl(-1)=-1.d0               / h
                Bl(0)= 2.d0                / h
                Bl(1)=-1.d0                / h

                A1(n_max-1)=-1.d0
                A2(n_max-1)=1.d0

                A2(2)=1.d0
                A3(2)=-1.d0

                if (allocated(m)) deallocate(m)
                if (allocated(bp)) deallocate(bp)

                allocate(m(2:n_max-1), bp(2:n_max-1))

                call TS_initDIR(A1(2:n_max-1), A2(2:n_max-1), A3(2:n_max-1), bp(2:n_max-1), m(2:n_max-1), n_max-2)

            case (symetric)

                A1(n_max-1)=al
                A2(n_max-1)=1.d0

                A2(2)=1.d0
                A3(2)=al

                if (allocated(m)) deallocate(m)
                if (allocated(bp)) deallocate(bp)

                allocate(m(2:n_max-1), bp(2:n_max-1))

                call TS_initDIR(A1(2:n_max-1), A2(2:n_max-1), A3(2:n_max-1), bp(2:n_max-1), m(2:n_max-1), n_max-2)

            case default

        end select

    end subroutine fill_matrix

    subroutine D1ssh_CPT_3Dx(f, d1f, n1,n2,n2e,n3,n3e, h, shifted, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: d1f

        real(kind=8), dimension(n2,n3)                      :: sx

        integer :: i,j,k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n1, shifted, bc_type, h)

        if (bc_type==periodic) then

            do k=1,n3e
                do j=1,n2e
                    d1f(1,j,k)= B(2)*(f(2,j,k)-f(n1-2,j,k)) + B(1)*(f(1,j,k) - f(n1-1,j,k))
                    d1f(2,j,k)= B(2)*(f(3,j,k)-f(n1-1,j,k)) + B(1)*(f(2,j,k) - f(1,j,k))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i=3,n1-2
                        d1f(i,j,k)= B(2)*(f(i+1,j,k)-f(i-2,j,k)) + B(1)*(f(i,j,k) - f(i-1,j,k))
                    enddo
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    d1f(n1-1,j,k)= B(2)*(f(1,j,k)-f(n1-3,j,k)) + B(1)*(f(n1-1,j,k) - f(n1-2,j,k))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i = 2,n1-1
                        d1f(i, j,k)=d1f(i, j,k) - m(i)*d1f(i-1, j,k)
                    enddo
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    d1f(n1-1, j,k) = d1f(n1-1, j,k)*bp(n1-1)
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i = n1-2, 1, -1
                        d1f(i, j,k)=(d1f(i, j,k) - A3(i)*d1f(i+1, j,k))*bp(i)
                    enddo
                enddo
            enddo

            do k=1, n3e
                do j=1, n2e
                    sx(j,k) = (v(1)*d1f(1,j,k)+v(n1-1)*d1f(n1-1,j,k))/(1+v(1)*q(1)+v(n1-1)*q(n1-1))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i=1, n1-1
                        d1f(i, j,k)=d1f(i, j,k)-sx(j,k)*q(i)
                    enddo
                enddo
            enddo

        end if

        if (bc_type==Dirichlet) then

            do k=1,n3e
                do j=1,n2e

                    do i=3,n1-2
                        d1f(i,j,k)= B(2)*(f(i+1,j,k)-f(i-2,j,k)) + B(1)*(f(i,j,k) - f(i-1,j,k))
                    enddo

                    d1f(2,j,k)= ( Bl(-1)*f(1,j,k) + Bl(0) * f(2,j,k) + Bl(1) * f(3,j,k) )

                    d1f(n1-1,j,k)= -( + Bl(-1)*f(n1-1,j,k) + Bl(0) * f(n1-2,j,k) + Bl(1)*f(n1-3,j,k) )

                    do i=3, n1-2
                        d1f(i, j,k)= d1f(i, j,k) - m(i)*d1f(i-1, j,k)
                    enddo

                    d1f(n1-1, j,k)= (d1f(n1-1, j,k) - m(n1-1)*d1f(n1-2, j,k))*bp(n1-1)

                    do i = n1-2, 2, -1
                        d1f(i, j,k) = (d1f(i, j,k) - A3(i)*d1f(i+1, j,k))*bp(i)
                    enddo

                enddo
            enddo

        endif

        if (bc_type==symetric) then

            do k=1,n3e
                do j=1,n2e

                    do i=3,n1-2
                        d1f(i,j,k)= B(2)*(f(i+1,j,k)-f(i-2,j,k)) + B(1)*(f(i,j,k) - f(i-1,j,k))
                    enddo
                    d1f(1,j,k)=0.d0
                    d1f(2,j,k)= B(2)*(f(3,j,k)-f(1,j,k)) + B(1)*(f(2,j,k) - f(1,j,k))
                    d1f(n1-1,j,k)= B(2)*(f(n1-1,j,k)-f(n1-3,j,k)) + B(1)*(f(n1-1,j,k) - f(n1-2,j,k))

                    do i=3, n1-2
                        d1f(i, j,k)= d1f(i, j,k) - m(i)*d1f(i-1, j,k)
                    enddo

                    d1f(n1-1, j,k)= (d1f(n1-1, j,k) - m(n1-1)*d1f(n1-2, j,k))*bp(n1-1)

                    do i = n1-2, 2, -1
                        d1f(i, j,k) = (d1f(i, j,k) - A3(i)*d1f(i+1, j,k))*bp(i)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D1ssh_CPT_3Dx

    subroutine D1ssh_CPT_ACC_3Dx(f, d1f, n1,n2,n2e,n3,n3e, h, shifted, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: d1f

        real*8, dimension(n1)           :: tmp
        real*8, dimension(n1,n2,n3)     :: tmp2

        real(kind=8), dimension(n2,n3)                      :: sx

        integer :: i,j,k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n1, shifted, bc_type, h)

        if (bc_type==periodic) then

            do k=1,n3e
                do j=1,n2e
                    tmp2(1,j,k)= B(2)*(f(2,j,k)-f(n1-2,j,k)) + B(1)*(f(1,j,k) - f(n1-1,j,k))
                    tmp2(2,j,k)= B(2)*(f(3,j,k)-f(n1-1,j,k)) + B(1)*(f(2,j,k) - f(1,j,k))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i=3,n1-2
                        tmp2(i,j,k)= B(2)*(f(i+1,j,k)-f(i-2,j,k)) + B(1)*(f(i,j,k) - f(i-1,j,k))
                    enddo
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    tmp2(n1-1,j,k)= B(2)*(f(1,j,k)-f(n1-3,j,k)) + B(1)*(f(n1-1,j,k) - f(n1-2,j,k))
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
                        d1f(i, j,k)=tmp2(i, j,k)-sx(j,k)*q(i) + d1f(i, j,k)
                    enddo
                enddo
            enddo

        end if

        if (bc_type==Dirichlet) then

            do k=1,n3e
                do j=1,n2e

                    tmp(2:n1-1)=d1f(2:n1-1,j,k)

                    do i=3,n1-2
                        d1f(i,j,k)= B(2)*(f(i+1,j,k)-f(i-2,j,k)) + B(1)*(f(i,j,k) - f(i-1,j,k))
                    enddo

                    d1f(2,j,k)= ( Bl(-1)*f(1,j,k) + Bl(0) * f(2,j,k) + Bl(1) * f(3,j,k) )
                    d1f(n1-1,j,k)= -( + Bl(-1)*f(n1-1,j,k) + Bl(0) * f(n1-2,j,k) + Bl(1)*f(n1-3,j,k) )

                    do i=3, n1-2
                        d1f(i, j,k)= d1f(i, j,k) - m(i)*d1f(i-1, j,k)
                    enddo

                    d1f(n1-1, j,k)= (d1f(n1-1, j,k) - m(n1-1)*d1f(n1-2, j,k))*bp(n1-1)

                    do i = n1-2, 2, -1
                        d1f(i, j,k) = (d1f(i, j,k) - A3(i)*d1f(i+1, j,k))*bp(i)
                    enddo

                    do i=2, n1-1
                        d1f(i, j,k)= d1f(i, j,k) + tmp(i)
                    enddo

                enddo
            enddo

        endif

        if (bc_type==symetric) then

            do k=1,n3e
                do j=1,n2e

                    tmp(2:n1-1)=d1f(2:n1-1,j,k)

                    do i=3,n1-2
                        d1f(i,j,k)= B(2)*(f(i+1,j,k)-f(i-2,j,k)) + B(1)*(f(i,j,k) - f(i-1,j,k))
                    enddo
                    d1f(1,j,k)=0.d0
                    d1f(2,j,k)= B(2)*(f(3,j,k)-f(1,j,k)) + B(1)*(f(2,j,k) - f(1,j,k))
                    d1f(n1-1,j,k)= B(2)*(f(n1-1,j,k)-f(n1-3,j,k)) + B(1)*(f(n1-1,j,k) - f(n1-2,j,k))

                    do i=3, n1-2
                        d1f(i, j,k)= d1f(i, j,k) - m(i)*d1f(i-1, j,k)
                    enddo

                    d1f(n1-1, j,k)= (d1f(n1-1, j,k) - m(n1-1)*d1f(n1-2, j,k))*bp(n1-1)

                    do i = n1-2, 2, -1
                        d1f(i, j,k) = (d1f(i, j,k) - A3(i)*d1f(i+1, j,k))*bp(i)
                    enddo

                    do i=2, n1-1
                        d1f(i, j,k)= d1f(i, j,k) + tmp(i)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D1ssh_CPT_ACC_3Dx



    subroutine D1ssh_CPT_3Dy(f, d1f, n1,n1e,n2,n3,n3e, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: d1f

        integer :: i,j,k

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


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=3,n2-2
                        d1f(i,j,k)= B(2)*(f(i,j+1,k)-f(i,j-2,k)) + B(1)*(f(i,j,k) - f(i,j-1,k))
                    enddo

                    d1f(i,1,k)= B(2)*(f(i,2,k)-f(i,n2-2,k)) + B(1)*(f(i,1,k) - f(i,n2-1,k))
                    d1f(i,2,k)= B(2)*(f(i,3,k)-f(i,n2-1,k)) + B(1)*(f(i,2,k) - f(i,1,k))
                    d1f(i,n2-1,k)= B(2)*(f(i,1,k)-f(i,n2-3,k)) + B(1)*(f(i,n2-1,k) - f(i,n2-2,k))

                    do j = 2,n2-1
                        d1f(i,j,k)=d1f(i,j,k) - m(j)*d1f(i,j-1,k)
                    enddo
                    d1f(i,n2-1,k) = d1f(i,n2-1,k)*bp(n2-1)
                    do j = n2-2, 1, -1
                        d1f(i,j,k)=(d1f(i,j,k) - A3(j)*d1f(i,j+1,k))*bp(j)
                    enddo

                    vy=0.d0
                    do j=1, n2-1
                        vy=vy+v(j)*d1f(i,j,k)
                    enddo
                    coef_q = vy/(1+vq)

                    do j=1, n2-1
                        d1f(i,j,k)=d1f(i,j,k)-coef_q*q(j)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.symetric)) then

            do k=1, n3e
                do i=1, n1e

                    do j=3,n2-2
                        d1f(i,j,k)= B(2)*(f(i,j+1,k)-f(i,j-2,k)) + B(1)*(f(i,j,k) - f(i,j-1,k))
                    enddo
                    d1f(i,1,k)=0.d0
                    d1f(i,2,k)= B(2)*(f(i,3,k)-f(i,1,k)) + B(1)*(f(i,2,k) - f(i,1,k))

                    d1f(i,n2-1,k)= B(2)*(f(i,n2-1,k)-f(i,n2-3,k)) + B(1)*(f(i,n2-1,k) - f(i,n2-2,k))


                    do j=3, n2-2
                        d1f(i,j,k)= d1f(i,j,k) - m(j)*d1f(i,j-1,k)
                    enddo

                    d1f(i,n2-1,k)= (d1f(i,n2-1,k) - m(n2-1)*d1f(i,n2-2,k))*bp(n2-1)


                    do j = n2-2, 2, -1
                        d1f(i,j,k) = (d1f(i,j,k) - A3(j)*d1f(i,j+1,k))*bp(j)
                    enddo

                enddo
            enddo

        end if


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e

                    do j=3,n2-2
                        d1f(i,j,k)= B(2)*(f(i,j+1,k)-f(i,j-2,k)) + B(1)*(f(i,j,k) - f(i,j-1,k))
                    enddo

                    d1f(i,2,k)= ( Bl(-1)*f(i,1,k) + Bl(0) * f(i,2,k) + Bl(1) * f(i,3,k) )

                    d1f(i,n2-1,k)= -( + Bl(-1)*f(i,n2-1,k) + Bl(0) * f(i,n2-2,k) + Bl(1)*f(i,n2-3,k) )


                    do j=3, n2-2
                        d1f(i,j,k)= d1f(i,j,k) - m(j)*d1f(i,j-1,k)
                    enddo

                    d1f(i,n2-1,k)= (d1f(i,n2-1,k) - m(n2-1)*d1f(i,n2-2,k))*bp(n2-1)


                    do j = n2-2, 2, -1
                        d1f(i,j,k) = (d1f(i,j,k) - A3(j)*d1f(i,j+1,k))*bp(j)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D1ssh_CPT_3Dy

    subroutine D1ssh_CPT_MULT_3Dy(f, d1f, n1,n1e,n2,n3,n3e, h, shifted, bc_type, g)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(:), intent(in)             :: g
        real(kind=8), dimension(n1,n2,n3), intent(in)       :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)      :: d1f

        integer :: i,j,k

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


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=3,n2-2
                        d1f(i,j,k)= B(2)*(f(i,j+1,k)-f(i,j-2,k)) + B(1)*(f(i,j,k) - f(i,j-1,k))
                    enddo

                    d1f(i,1,k)= B(2)*(f(i,2,k)-f(i,n2-2,k)) + B(1)*(f(i,1,k) - f(i,n2-1,k))
                    d1f(i,2,k)= B(2)*(f(i,3,k)-f(i,n2-1,k)) + B(1)*(f(i,2,k) - f(i,1,k))
                    d1f(i,n2-1,k)= B(2)*(f(i,1,k)-f(i,n2-3,k)) + B(1)*(f(i,n2-1,k) - f(i,n2-2,k))

                    do j = 2,n2-1
                        d1f(i,j,k)=d1f(i,j,k) - m(j)*d1f(i,j-1,k)
                    enddo
                    d1f(i,n2-1,k) = d1f(i,n2-1,k)*bp(n2-1)
                    do j = n2-2, 1, -1
                        d1f(i,j,k)=(d1f(i,j,k) - A3(j)*d1f(i,j+1,k))*bp(j)
                    enddo

                    vy=0.d0
                    do j=1, n2-1
                        vy=vy+v(j)*d1f(i,j,k)
                    enddo
                    coef_q = vy/(1+vq)

                    do j=1, n2-1
                        d1f(i,j,k)=(d1f(i,j,k)-coef_q*q(j))*g(j)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e

                    do j=3,n2-2
                        d1f(i,j,k)= B(2)*(f(i,j+1,k)-f(i,j-2,k)) + B(1)*(f(i,j,k) - f(i,j-1,k))
                    enddo

                    d1f(i,2,k)= ( Bl(-1)*f(i,1,k) + Bl(0) * f(i,2,k) + Bl(1) * f(i,3,k) )

                    d1f(i,n2-1,k)= -( + Bl(-1)*f(i,n2-1,k) + Bl(0) * f(i,n2-2,k) + Bl(1)*f(i,n2-3,k) )


                    do j=3, n2-2
                        d1f(i,j,k)= d1f(i,j,k) - m(j)*d1f(i,j-1,k)
                    enddo

                    d1f(i,n2-1,k)= (d1f(i,n2-1,k) - m(n2-1)*d1f(i,n2-2,k))*bp(n2-1)


                    do j = n2-2, 2, -1
                        d1f(i,j,k) = (d1f(i,j,k) - A3(j)*d1f(i,j+1,k))*bp(j)
                    enddo

                    do j=2, n2-1
                        d1f(i,j,k)= d1f(i,j,k)*g(j)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.symetric)) then

            do k=1, n3e
                do i=1, n1e

                    do j=3,n2-2
                        d1f(i,j,k)= B(2)*(f(i,j+1,k)-f(i,j-2,k)) + B(1)*(f(i,j,k) - f(i,j-1,k))
                    enddo
                    d1f(i,1,k)=0.d0
                    d1f(i,2,k)= B(2)*(f(i,3,k)-f(i,1,k)) + B(1)*(f(i,2,k) - f(i,1,k))

                    d1f(i,n2-1,k)= B(2)*(f(i,n2-1,k)-f(i,n2-3,k)) + B(1)*(f(i,n2-1,k) - f(i,n2-2,k))


                    do j=3, n2-2
                        d1f(i,j,k)= d1f(i,j,k) - m(j)*d1f(i,j-1,k)
                    enddo

                    d1f(i,n2-1,k)= (d1f(i,n2-1,k) - m(n2-1)*d1f(i,n2-2,k))*bp(n2-1)


                    do j = n2-2, 2, -1
                        d1f(i,j,k) = (d1f(i,j,k) - A3(j)*d1f(i,j+1,k))*bp(j)
                    enddo

                    do j=2, n2-1
                        d1f(i,j,k)= d1f(i,j,k)*g(j)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D1ssh_CPT_MULT_3Dy

    subroutine D1ssh_CPT_MULTACC_3Dy(f, d1f, n1,n1e,n2,n3,n3e, h, shifted, bc_type, g)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(:), intent(in)             :: g
        real(kind=8), dimension(n1,n2,n3), intent(in)       :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)      :: d1f

        real(kind=8), dimension(n2)         :: tmp
        real(kind=8), dimension(n1,n2,n3)   :: tmp2

        integer :: i,j,k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n2, shifted, bc_type, h)

        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    tmp(1:n2-1)=d1f(i,1:n2-1,k)

                    do j=3,n2-2
                        d1f(i,j,k)= B(2)*(f(i,j+1,k)-f(i,j-2,k)) + B(1)*(f(i,j,k) - f(i,j-1,k))
                    enddo

                    d1f(i,1,k)= B(2)*(f(i,2,k)-f(i,n2-2,k)) + B(1)*(f(i,1,k) - f(i,n2-1,k))
                    d1f(i,2,k)= B(2)*(f(i,3,k)-f(i,n2-1,k)) + B(1)*(f(i,2,k) - f(i,1,k))
                    d1f(i,n2-1,k)= B(2)*(f(i,1,k)-f(i,n2-3,k)) + B(1)*(f(i,n2-1,k) - f(i,n2-2,k))

                    do j = 2,n2-1
                        d1f(i,j,k)=d1f(i,j,k) - m(j)*d1f(i,j-1,k)
                    enddo
                    d1f(i,n2-1,k) = d1f(i,n2-1,k)*bp(n2-1)
                    do j = n2-2, 1, -1
                        d1f(i,j,k)=(d1f(i,j,k) - A3(j)*d1f(i,j+1,k))*bp(j)
                    enddo

                    vy=0.d0
                    do j=1, n2-1
                        vy=vy+v(j)*d1f(i,j,k)
                    enddo
                    coef_q = vy/(1+vq)

                    do j=1, n2-1
                        d1f(i,j,k)=(d1f(i,j,k)-coef_q*q(j))*g(j) + tmp(j)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e

                    do j=3,n2-2
                        tmp2(i,j,k)= B(2)*(f(i,j+1,k)-f(i,j-2,k)) + B(1)*(f(i,j,k) - f(i,j-1,k))
                    enddo

                    tmp2(i,2,k)= ( Bl(-1)*f(i,1,k) + Bl(0) * f(i,2,k) + Bl(1) * f(i,3,k) )

                    tmp2(i,n2-1,k)= -( + Bl(-1)*f(i,n2-1,k) + Bl(0) * f(i,n2-2,k) + Bl(1)*f(i,n2-3,k) )


                    do j=3, n2-2
                        tmp2(i,j,k)= tmp2(i,j,k) - m(j)*tmp2(i,j-1,k)
                    enddo

                    tmp2(i,n2-1,k)= (tmp2(i,n2-1,k) - m(n2-1)*tmp2(i,n2-2,k))*bp(n2-1)


                    do j = n2-2, 2, -1
                        tmp2(i,j,k) = (tmp2(i,j,k) - A3(j)*tmp2(i,j+1,k))*bp(j)
                    enddo

                    do j=2, n2-1
                        d1f(i,j,k)= tmp2(i,j,k)*g(j) + d1f(i,j,k)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.symetric)) then

            do k=1, n3e
                do i=1, n1e

                    do j=3,n2-2
                        tmp2(i,j,k)= B(2)*(f(i,j+1,k)-f(i,j-2,k)) + B(1)*(f(i,j,k) - f(i,j-1,k))
                    enddo

                    tmp2(i,1,k)=0.d0
                    tmp2(i,2,k)= B(2)*(f(i,3,k)-f(i,1,k)) + B(1)*(f(i,2,k) - f(i,1,k))

                    tmp2(i,n2-1,k)= B(2)*(f(i,n2-1,k)-f(i,n2-3,k)) + B(1)*(f(i,n2-1,k) - f(i,n2-2,k))

                    do j=3, n2-2
                        tmp2(i,j,k)= tmp2(i,j,k) - m(j)*tmp2(i,j-1,k)
                    enddo

                    tmp2(i,n2-1,k)= (tmp2(i,n2-1,k) - m(n2-1)*tmp2(i,n2-2,k))*bp(n2-1)


                    do j = n2-2, 2, -1
                        tmp2(i,j,k) = (tmp2(i,j,k) - A3(j)*tmp2(i,j+1,k))*bp(j)
                    enddo

                    do j=2, n2-1
                        d1f(i,j,k)= tmp2(i,j,k)*g(j) + d1f(i,j,k)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D1ssh_CPT_MULTACC_3Dy

    subroutine D1ssh_CPT_3Dz(f, d1f, n1,n1e,n2,n2e,n3, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: d1f

        real(kind=8), dimension(n1,n2)                      :: sz

        integer :: i,j,k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n3, shifted, bc_type, h)

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,1)= B(2)*(f(i,j,2)-f(i,j,n3-2)) + B(1)*(f(i,j,1) - f(i,j,n3-1))
                    d1f(i,j,2)= B(2)*(f(i,j,3)-f(i,j,n3-1)) + B(1)*(f(i,j,2) - f(i,j,1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k=3,n3-2
                        d1f(i,j,k)= B(2)*(f(i,j,k+1)-f(i,j,k-2)) + B(1)*(f(i,j,k) - f(i,j,k-1))
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,n3-1)= B(2)*(f(i,j,1)-f(i,j,n3-3)) + B(1)*(f(i,j,n3-1) - f(i,j,n3-2))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k = 2,n3-1
                        d1f(i,j,k)=d1f(i,j,k) - m(k)*d1f(i,j,k-1)
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,n3-1) = d1f(i,j,n3-1)*bp(n3-1)
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k = n3-2, 1, -1
                        d1f(i,j,k)=(d1f(i,j,k) - A3(k)*d1f(i,j,k+1))*bp(k)
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    sz(i,j) = (v(1)*d1f(i,j,1)+v(n3-1)*d1f(i,j,n3-1))/(1+v(1)*q(1)+v(n3-1)*q(n3-1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k=1, n3-1
                        d1f(i,j,k)=d1f(i,j,k)-sz(i,j)*q(k)
                    enddo
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e

                    do k=3,n3-2
                        d1f(i,j,k)= B(2)*(f(i,j,k+1)-f(i,j,k-2)) + B(1)*(f(i,j,k) - f(i,j,k-1))
                    enddo

                    d1f(i,j,2)= ( Bl(-1)*f(i,j,1) + Bl(0) * f(i,j,2) + Bl(1) * f(i,j,3) )
                    d1f(i,j,n3-1)= -( + Bl(-1)*f(i,j,n3-1) + Bl(0) * f(i,j,n3-2) + Bl(1)*f(i,j,n3-3) )


                    do k=3, n3-2
                        d1f(i,j,k)= d1f(i,j,k) - m(k)*d1f(i,j,k-1)
                    enddo

                    d1f(i,j,n3-1)= (d1f(i,j,n3-1) - m(n3-1)*d1f(i,j,n3-2))*bp(n3-1)


                    do k = n3-2, 2, -1
                        d1f(i,j,k) = (d1f(i,j,k) - A3(k)*d1f(i,j,k+1))*bp(k)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D1ssh_CPT_3Dz

    subroutine D1ssh_CPT_ACC_3Dz(f, d1f, n1,n1e,n2,n2e,n3, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: d1f

        real(kind=8), dimension(n3)         :: tmp
        real(kind=8), dimension(n1,n2,n3)   :: tmp2

        real(kind=8), dimension(n1,n2)                      :: sz

        integer :: i,j,k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n3, shifted, bc_type, h)

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e
                    tmp2(i,j,1)= B(2)*(f(i,j,2)-f(i,j,n3-2)) + B(1)*(f(i,j,1) - f(i,j,n3-1))
                    tmp2(i,j,2)= B(2)*(f(i,j,3)-f(i,j,n3-1)) + B(1)*(f(i,j,2) - f(i,j,1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k=3,n3-2
                        tmp2(i,j,k)= B(2)*(f(i,j,k+1)-f(i,j,k-2)) + B(1)*(f(i,j,k) - f(i,j,k-1))
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    tmp2(i,j,n3-1)= B(2)*(f(i,j,1)-f(i,j,n3-3)) + B(1)*(f(i,j,n3-1) - f(i,j,n3-2))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k = 2,n3-1
                        tmp2(i,j,k)=tmp2(i,j,k) - m(k)*tmp2(i,j,k-1)
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    tmp2(i,j,n3-1) = tmp2(i,j,n3-1)*bp(n3-1)
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k = n3-2, 1, -1
                        tmp2(i,j,k)=(tmp2(i,j,k) - A3(k)*tmp2(i,j,k+1))*bp(k)
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    sz(i,j) = (v(1)*tmp2(i,j,1)+v(n3-1)*tmp2(i,j,n3-1))/(1+v(1)*q(1)+v(n3-1)*q(n3-1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k=1, n3-1
                        d1f(i,j,k)=tmp2(i,j,k)-sz(i,j)*q(k) + d1f(i,j,k)
                    enddo
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e

                    tmp(2:n3-1)=d1f(i,j,2:n3-1)

                    do k=3,n3-2
                        d1f(i,j,k)= B(2)*(f(i,j,k+1)-f(i,j,k-2)) + B(1)*(f(i,j,k) - f(i,j,k-1))
                    enddo

                    d1f(i,j,2)= ( Bl(-1)*f(i,j,1) + Bl(0) * f(i,j,2) + Bl(1) * f(i,j,3) )
                    d1f(i,j,n3-1)= -( + Bl(-1)*f(i,j,n3-1) + Bl(0) * f(i,j,n3-2) + Bl(1)*f(i,j,n3-3) )


                    do k=3, n3-2
                        d1f(i,j,k)= d1f(i,j,k) - m(k)*d1f(i,j,k-1)
                    enddo

                    d1f(i,j,n3-1)= (d1f(i,j,n3-1) - m(n3-1)*d1f(i,j,n3-2))*bp(n3-1)


                    do k = n3-2, 2, -1
                        d1f(i,j,k) = (d1f(i,j,k) - A3(k)*d1f(i,j,k+1))*bp(k)
                    enddo

                    do k=2, n3-1
                        d1f(i,j,k)= d1f(i,j,k) + tmp(k)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D1ssh_CPT_ACC_3Dz

end module CPT_D1s_sh




