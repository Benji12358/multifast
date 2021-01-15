module CPT_D0s
    use boundaries_types
    implicit none

    private

    real(kind=8), allocatable, dimension(:) :: q,v
    real(kind=8), allocatable, dimension(:) :: m, bp, cp
    real(kind=8), dimension(:), allocatable :: A1, A2, A3
    real(kind=8)                                 :: vy, vq, coef_q
    real(kind=8)    :: B(2), Bl(-1:1)

    public::D0s_CPT_3Dx,D0s_CPT_3Dy, D0s_CPT_3Dz

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

        al=3.00001d0/10.d0
        B(1)=(10.d0*al+9.d0)/16.d0
        B(2)= (6.d0*al-1)/16.d0

        do i= 1, n_max
            A1(i) = al
            A2(i) = 1.d0
            A3(i) = al
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

                Bl(-1)=1.d0/4.d0
                Bl(0)= 3.d0/2.d0
                Bl(1)=1.d0/4.d0

                A1(n_max-1)=1.d0
                A2(n_max-1)=1.d0

                A2(1)=1.d0
                A3(1)=1.d0

                if (allocated(m)) then
                    deallocate(m, bp)
                end if

                allocate(m(n_max-1), bp(n_max-1))

                call TS_initDIR(A1(1:n_max-1), A2(1:n_max-1), A3(1:n_max-1), bp, m, n_max-1)

            case default

        end select

    end subroutine fill_matrix

    subroutine D0s_CPT_3Dx(f, ff, n1,n2,n2e,n3,n3e, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: ff

        real(kind=8), dimension(n2,n3)                      :: sx

        integer :: i, j, k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n1, .false., bc_type, 0.d0)

        if (bc_type==periodic) then

            do k=1,n3e
                do j=1,n2e
                    ff(1,j,k)= B(2)*(f(3,j,k)+f(n1-1,j,k)) + B(1)*(f(2,j,k)+f(1,j,k))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i=2,n1-3
                        ff(i,j,k)= B(2)*(f(i+2,j,k)+f(i-1,j,k)) + B(1)*(f(i+1,j,k) + f(i,j,k))
                    enddo
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    ff(n1-2,j,k)= B(2)*(f(1,j,k)+f(n1-3,j,k)) + B(1)*(f(n1-1,j,k)+f(n1-2,j,k))
                    ff(n1-1,j,k)= B(2)*(f(2,j,k)+f(n1-2,j,k)) + B(1)*(f(1,j,k)+f(n1-1,j,k))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i = 2,n1-1
                        ff(i, j,k)=ff(i, j,k) - m(i)*ff(i-1, j,k)
                    enddo
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    ff(n1-1, j,k) = ff(n1-1, j,k)*bp(n1-1)
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i = n1-2, 1, -1
                        ff(i, j,k)=(ff(i, j,k) - A3(i)*ff(i+1, j,k))*bp(i)
                    enddo
                enddo
            enddo

            do k=1, n3e
                do j=1, n2e
                    sx(j,k) = (v(1)*ff(1,j,k)+v(n1-1)*ff(n1-1,j,k))/(1+v(1)*q(1)+v(n1-1)*q(n1-1))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i=1, n1-1
                        ff(i, j,k)=ff(i, j,k)-sx(j,k)*q(i)
                    enddo
                enddo
            enddo

        end if

        if (bc_type==Dirichlet) then

            do k=1,n3e
                do j=1,n2e

                    do i=2,n1-2
                        ff(i,j,k)= B(2)*(f(i+2,j,k)+f(i-1,j,k)) + B(1)*(f(i+1,j,k)+f(i,j,k))
                    enddo

                    ff(1,j,k)= ( Bl(-1)*f(1,j,k) + Bl(0) * f(2,j,k) + Bl(1) * f(3,j,k) )
                    ff(n1-1,j,k)= Bl(-1)*f(n1,j,k) + Bl(0) * f(n1-1,j,k) + Bl(1)*f(n1-2,j,k)


                    do i=2, n1-2
                        ff(i, j,k)= ff(i, j,k) - m(i)*ff(i-1, j,k)
                    enddo

                    ff(n1-1, j,k)= (ff(n1-1, j,k) - m(n1-1)*ff(n1-2, j,k))*bp(n1-1)


                    do i = n1-2, 1, -1
                        ff(i, j,k) = (ff(i, j,k) - A3(i)*ff(i+1, j,k))*bp(i)
                    enddo

                enddo
            enddo

        end if


        return

    end subroutine D0s_CPT_3Dx



    subroutine D0s_CPT_3Dy(f, ff, n1,n1e,n2,n3,n3e, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: ff

        integer :: i,j,k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: update_matrix



        update_matrix=(last_n/=n2).or.(last_bc_type/=bc_type)


        !write(*,*)'D1_OUCS3'

        ! MATRIX A, B definition -----------------------------------------------------------
        call fill_matrix(n2, .false., bc_type, 0.d0)

        last_n=n2
        last_bc_type=bc_type


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=2,n2-2
                        ff(i,j,k)= B(2)*(f(i,j+2,k)+f(i,j-1,k)) + B(1)*(f(i,j+1,k)+f(i,j,k))
                    enddo

                    ff(i,1,k)= B(2)*(f(i,3,k)+f(i,n2-1,k)) + B(1)*(f(i,2,k)+f(i,1,k))

                    ff(i,n2-1,k)= B(2)*(f(i,2,k)+f(i,n2-2,k)) + B(1)*(f(i,1,k)+f(i,n2-1,k))

                    do j = 2,n2-1
                        ff(i,j,k)=ff(i,j,k) - m(j)*ff(i,j-1,k)
                    enddo
                    ff(i,n2-1,k) = ff(i,n2-1,k)*bp(n2-1)
                    do j = n2-2, 1, -1
                        ff(i,j,k)=(ff(i,j,k) - A3(j)*ff(i,j+1,k))*bp(j)
                    enddo

                    vy=0.d0
                    do j=1, n2-1
                        vy=vy+v(j)*ff(i,j,k)
                    enddo
                    coef_q = vy/(1+vq)

                    do j=1, n2-1
                        ff(i,j,k)=ff(i,j,k)-coef_q*q(j)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e

                    do j=2,n2-2
                        ff(i,j,k)= B(2)*(f(i,j+2,k)+f(i,j-1,k)) + B(1)*(f(i,j+1,k)+f(i,j,k))
                    enddo

                    ff(i,1,k)= ( Bl(-1)*f(i,1,k) + Bl(0) * f(i,2,k) + Bl(1) * f(i,3,k) )
                    ff(i,n2-1,k)= ( + Bl(-1)*f(i,n2,k) + Bl(0) * f(i,n2-1,k) + Bl(1)*f(i,n2-2,k) )


                    do j=2, n2-2
                        ff(i,j,k)= ff(i,j,k) - m(j)*ff(i,j-1,k)
                    enddo

                    ff(i,n2-1,k)= (ff(i,n2-1,k) - m(n2-1)*ff(i,n2-2,k))*bp(n2-1)


                    do j = n2-2, 1, -1
                        ff(i,j,k) = (ff(i,j,k) - A3(j)*ff(i,j+1,k))*bp(j)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D0s_CPT_3Dy

    subroutine D0s_CPT_3Dz(f, ff, n1,n1e,n2,n2e,n3, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: ff

        real(kind=8), dimension(n1,n2)                      :: sz

        integer :: i,j,k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n3, .false., bc_type, 0.d0)

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e
                    ff(i,j,1)= B(2)*(f(i,j,3)+f(i,j,n3-1)) + B(1)*(f(i,j,2)+f(i,j,1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e

                    do k=2,n3-3
                        ff(i,j,k)= B(2)*(f(i,j,k+2)+f(i,j,k-1)) + B(1)*(f(i,j,k+1)+f(i,j,k))
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    ff(i,j,n3-2)= B(2)*(f(i,j,1)+f(i,j,n3-3)) + B(1)*(f(i,j,n3-1)+f(i,j,n3-2))
                    ff(i,j,n3-1)= B(2)*(f(i,j,2)+f(i,j,n3-2)) + B(1)*(f(i,j,1)+f(i,j,n3-1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k = 2,n3-1
                        ff(i,j,k)=ff(i,j,k) - m(k)*ff(i,j,k-1)
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    ff(i,j,n3-1) = ff(i,j,n3-1)*bp(n3-1)
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k = n3-2, 1, -1
                        ff(i,j,k)=(ff(i,j,k) - A3(k)*ff(i,j,k+1))*bp(k)
                    enddo
                enddo
            enddo


            do j=1, n2e
                do i=1, n1e
                    sz(i,j) = (v(1)*ff(i,j,1)+v(n3-1)*ff(i,j,n3-1))/(1+v(1)*q(1)+v(n3-1)*q(n3-1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k=1, n3-1
                        ff(i,j,k)=ff(i,j,k)-sz(i,j)*q(k)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e

                    do k=2,n3-2
                        ff(i,j,k)= B(2)*(f(i,j,k+2)+f(i,j,k-1)) + B(1)*(f(i,j,k+1)+f(i,j,k))
                    enddo

                    ff(i,j,1)= ( Bl(-1)*f(i,j,1) + Bl(0) * f(i,j,2) + Bl(1) * f(i,j,3) )
                    ff(i,j,n3-1)= ( + Bl(-1)*f(i,j,n3) + Bl(0) * f(i,j,n3-1) + Bl(1)*f(i,j,n3-2) )


                    do k=2, n3-2
                        ff(i,j,k)= ff(i,j,k) - m(k)*ff(i,j,k-1)
                    enddo

                    ff(i,j,n3-1)= (ff(i,j,n3-1) - m(n3-1)*ff(i,j,n3-2))*bp(n3-1)


                    do k = n3-2, 1, -1
                        ff(i,j,k) = (ff(i,j,k) - A3(k)*ff(i,j,k+1))*bp(k)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D0s_CPT_3Dz

end module CPT_D0s

module CPT_D0s_sh
    use boundaries_types
    implicit none

    private

    real(kind=8), allocatable, dimension(:) :: q,v
    real(kind=8), allocatable, dimension(:) :: m, bp, cp
    real(kind=8), dimension(:), allocatable :: A1, A2, A3
    real(kind=8)                                 :: vy, vq, coef_q
    real(kind=8)    :: B(2), Bl(-1:1)

    public::D0ssh_CPT_3Dx, D0ssh_CPT_3Dy,D0ssh_CPT_MULT_3Dy, D0ssh_CPT_3Dz,D0ssh_CPT_MULT_3Dz

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

        al=3.00001d0/10.d0
        B(1)=(10.d0*al+9.d0)/16.d0
        B(2)= (6.d0*al-1)/16.d0


        do i= 1, n_max
            A1(i) = al
            A2(i) = 1.d0
            A3(i) = al
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

                Bl(-1)=1.d0/4.d0
                Bl(0)= 3.d0/2.d0
                Bl(1)=1.d0/4.d0

                A1(n_max-1)=1.d0
                A2(n_max-1)=1.d0

                A2(2)=1.d0
                A3(2)=1.d0

                if (allocated(m)) then
                    deallocate(m, bp)
                end if

                allocate(m(2:n_max-1), bp(2:n_max-1))

                call TS_initDIR(A1(2:n_max-1), A2(2:n_max-1), A3(2:n_max-1), bp(2:n_max-1), m(2:n_max-1), n_max-2)

            case default

        end select

    end subroutine fill_matrix



    subroutine D0ssh_CPT_3Dx(f, ff, n1,n2,n2e,n3,n3e, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: ff

        real(kind=8), dimension(n2,n3) :: sz

        integer :: i,j,k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n1, .false., bc_type, 0.d0)

        if (bc_type==periodic) then

            do k=1,n3e
                do j=1,n2e
                    ff(1,j,k)= B(2)*(f(2,j,k)+f(n1-2,j,k)) + B(1)*(f(1,j,k)+f(n1-1,j,k))
                    ff(2,j,k)= B(2)*(f(3,j,k)+f(n1-1,j,k)) + B(1)*(f(2,j,k)+f(1,j,k))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i=3,n1-2
                        ff(i,j,k)= B(2)*(f(i+1,j,k)+f(i-2,j,k)) + B(1)*(f(i,j,k)+f(i-1,j,k))
                    enddo
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    ff(n1-1,j,k)= B(2)*(f(1,j,k)+f(n1-3,j,k)) + B(1)*(f(n1-1,j,k)+f(n1-2,j,k))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i = 2,n1-1
                        ff(i, j,k)=ff(i, j,k) - m(i)*ff(i-1, j,k)
                    enddo
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    ff(n1-1, j,k) = ff(n1-1, j,k)*bp(n1-1)
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i = n1-2, 1, -1
                        ff(i, j,k)=(ff(i, j,k) - A3(i)*ff(i+1, j,k))*bp(i)
                    enddo
                enddo
            enddo

            do k=1, n3e
                do j=1, n2e
                    sz(j,k) = (v(1)*ff(1,j,k)+v(n1-1)*ff(n1-1,j,k))/(1+v(1)*q(1)+v(n1-1)*q(n1-1))
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    do i=1, n1-1
                        ff(i, j,k)=ff(i, j,k)-sz(j,k)*q(i)
                    enddo
                enddo
            enddo

        end if

        if (bc_type==Dirichlet) then

            do k=1,n3e
                do j=1,n2e

                    do i=3,n1-2
                        ff(i,j,k)= B(2)*(f(i+1,j,k)+f(i-2,j,k)) + B(1)*(f(i,j,k)+f(i-1,j,k))
                    enddo

                    ff(2,j,k)= ( Bl(-1)*f(1,j,k) + Bl(0) * f(2,j,k) + Bl(1) * f(3,j,k) )
                    ff(n1-1,j,k)= ( + Bl(-1)*f(n1-1,j,k) + Bl(0) * f(n1-2,j,k) + Bl(1)*f(n1-3,j,k) )

                    do i=3, n1-2
                        ff(i, j,k)= ff(i, j,k) - m(i)*ff(i-1, j,k)
                    enddo

                    ff(n1-1, j,k)= (ff(n1-1, j,k) - m(n1-1)*ff(n1-2, j,k))*bp(n1-1)

                    do i = n1-2, 2, -1
                        ff(i, j,k) = (ff(i, j,k) - A3(i)*ff(i+1, j,k))*bp(i)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D0ssh_CPT_3Dx



    subroutine D0ssh_CPT_3Dy(f, ff, n1,n1e,n2,n3,n3e, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: ff

        integer :: i,j,k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: update_matrix



        update_matrix=(last_n/=n2).or.(last_bc_type/=bc_type)


        !write(*,*)'D1_OUCS3'

        ! MATRIX A, B definition -----------------------------------------------------------
        call fill_matrix(n2, .true., bc_type, 0.d0)

        last_n=n2
        last_bc_type=bc_type

        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=3,n2-2
                        ff(i,j,k)= B(2)*(f(i,j+1,k)+f(i,j-2,k)) + B(1)*(f(i,j,k)+f(i,j-1,k))
                    enddo

                    ff(i,1,k)= B(2)*(f(i,2,k)+f(i,n2-2,k)) + B(1)*(f(i,1,k)+f(i,n2-1,k))
                    ff(i,2,k)= B(2)*(f(i,3,k)+f(i,n2-1,k)) + B(1)*(f(i,2,k)+f(i,1,k))
                    ff(i,n2-1,k)= B(2)*(f(i,1,k)+f(i,n2-3,k)) + B(1)*(f(i,n2-1,k)+f(i,n2-2,k))

                    do j = 2,n2-1
                        ff(i,j,k)=ff(i,j,k) - m(j)*ff(i,j-1,k)
                    enddo
                    ff(i,n2-1,k) = ff(i,n2-1,k)*bp(n2-1)
                    do j = n2-2, 1, -1
                        ff(i,j,k)=(ff(i,j,k) - A3(j)*ff(i,j+1,k))*bp(j)
                    enddo

                    vy=0.d0
                    do j=1, n2-1
                        vy=vy+v(j)*ff(i,j,k)
                    enddo
                    coef_q = vy/(1+vq)

                    do j=1, n2-1
                        ff(i,j,k)=ff(i,j,k)-coef_q*q(j)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e

                    do j=3,n2-2
                        ff(i,j,k)= B(2)*(f(i,j+1,k)+f(i,j-2,k)) + B(1)*(f(i,j,k)+f(i,j-1,k))
                    enddo

                    ff(i,2,k)= ( Bl(-1)*f(i,1,k) + Bl(0) * f(i,2,k) + Bl(1) * f(i,3,k) )
                    ff(i,n2-1,k)= ( + Bl(-1)*f(i,n2-1,k) + Bl(0) * f(i,n2-2,k) + Bl(1)*f(i,n2-3,k) )


                    do j=3, n2-2
                        ff(i,j,k)= ff(i,j,k) - m(j)*ff(i,j-1,k)
                    enddo

                    ff(i,n2-1,k)= (ff(i,n2-1,k) - m(n2-1)*ff(i,n2-2,k))*bp(n2-1)


                    do j = n2-2, 2, -1
                        ff(i,j,k) = (ff(i,j,k) - A3(j)*ff(i,j+1,k))*bp(j)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D0ssh_CPT_3Dy



    subroutine D0ssh_CPT_MULT_3Dy(f, ff, n1,n1e,n2,n3,n3e, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: ff

        real(kind=8), dimension(n2)         :: tmp
        real(kind=8), dimension(n1,n2,n3)   :: tmp2

        integer :: i,j,k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: update_matrix

        call fill_matrix(n2, .true., bc_type, 0.d0)

        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    tmp(1:n2-1)=ff(i,1:n2-1,k)

                    do j=3,n2-2
                        ff(i,j,k)= B(2)*(f(i,j+1,k)+f(i,j-2,k)) + B(1)*(f(i,j,k)+f(i,j-1,k))
                    enddo

                    ff(i,1,k)= B(2)*(f(i,2,k)+f(i,n2-2,k)) + B(1)*(f(i,1,k)+f(i,n2-1,k))
                    ff(i,2,k)= B(2)*(f(i,3,k)+f(i,n2-1,k)) + B(1)*(f(i,2,k)+f(i,1,k))
                    ff(i,n2-1,k)= B(2)*(f(i,1,k)+f(i,n2-3,k)) + B(1)*(f(i,n2-1,k)+f(i,n2-2,k))

                    do j = 2,n2-1
                        ff(i,j,k)=ff(i,j,k) - m(j)*ff(i,j-1,k)
                    enddo
                    ff(i,n2-1,k) = ff(i,n2-1,k)*bp(n2-1)
                    do j = n2-2, 1, -1
                        ff(i,j,k)=(ff(i,j,k) - A3(j)*ff(i,j+1,k))*bp(j)
                    enddo

                    vy=0.d0
                    do j=1, n2-1
                        vy=vy+v(j)*ff(i,j,k)
                    enddo
                    coef_q = vy/(1+vq)

                    do j=1, n2-1
                        ff(i,j,k)=(ff(i,j,k)-coef_q*q(j))*tmp(j)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e

                    do j=3,n2-2
                        tmp2(i,j,k)= B(2)*(f(i,j+1,k)+f(i,j-2,k)) + B(1)*(f(i,j,k)+f(i,j-1,k))
                    enddo

                    tmp2(i,2,k)= ( Bl(-1)*f(i,1,k) + Bl(0) * f(i,2,k) + Bl(1) * f(i,3,k) )
                    tmp2(i,n2-1,k)= ( + Bl(-1)*f(i,n2-1,k) + Bl(0) * f(i,n2-2,k) + Bl(1)*f(i,n2-3,k) )


                    do j=3, n2-2
                        tmp2(i,j,k)= tmp2(i,j,k) - m(j)*tmp2(i,j-1,k)
                    enddo

                    tmp2(i,n2-1,k)= (tmp2(i,n2-1,k) - m(n2-1)*tmp2(i,n2-2,k))*bp(n2-1)


                    do j = n2-2, 2, -1
                        tmp2(i,j,k) = (tmp2(i,j,k) - A3(j)*tmp2(i,j+1,k))*bp(j)
                    enddo

                    do j=2, n2-1
                        ff(i,j,k)= tmp2(i,j,k)*ff(i,j,k)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D0ssh_CPT_MULT_3Dy

    subroutine D0ssh_CPT_3Dz(f, ff, n1,n1e,n2,n2e,n3, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: ff

        real(kind=8), dimension(n1,n2)                      :: sz

        integer :: i,j,k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n3, .true., bc_type, 0.d0)

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e
                    ff(i,j,1)= B(2)*(f(i,j,2)+f(i,j,n3-2)) + B(1)*(f(i,j,1)+f(i,j,n3-1))
                    ff(i,j,2)= B(2)*(f(i,j,3)+f(i,j,n3-1)) + B(1)*(f(i,j,2)+f(i,j,1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k=3,n3-2
                        ff(i,j,k)= B(2)*(f(i,j,k+1)+f(i,j,k-2)) + B(1)*(f(i,j,k)+f(i,j,k-1))
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    ff(i,j,n3-1)= B(2)*(f(i,j,1)+f(i,j,n3-3)) + B(1)*(f(i,j,n3-1)+f(i,j,n3-2))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k = 2,n3-1
                        ff(i,j,k)=ff(i,j,k) - m(k)*ff(i,j,k-1)
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    ff(i,j,n3-1) = ff(i,j,n3-1)*bp(n3-1)
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k = n3-2, 1, -1
                        ff(i,j,k)=(ff(i,j,k) - A3(k)*ff(i,j,k+1))*bp(k)
                    enddo
                enddo
            enddo


            do j=1, n2e
                do i=1, n1e
                    sz(i,j) = (v(1)*ff(i,j,1)+v(n3-1)*ff(i,j,n3-1))/(1+v(1)*q(1)+v(n3-1)*q(n3-1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k=1, n3-1
                        ff(i,j,k)=ff(i,j,k)-sz(i,j)*q(k)
                    enddo
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e

                    do k=3,n3-2
                        ff(i,j,k)= B(2)*(f(i,j,k+1)+f(i,j,k-2)) + B(1)*(f(i,j,k)+f(i,j,k-1))
                    enddo

                    ff(i,j,2)= ( Bl(-1)*f(i,j,1) + Bl(0) * f(i,j,2) + Bl(1) * f(i,j,3) )
                    ff(i,j,n3-1)= ( + Bl(-1)*f(i,j,n3-1) + Bl(0) * f(i,j,n3-2) + Bl(1)*f(i,j,n3-3) )


                    do k=3, n3-2
                        ff(i,j,k)= ff(i,j,k) - m(k)*ff(i,j,k-1)
                    enddo

                    ff(i,j,n3-1)= (ff(i,j,n3-1) - m(n3-1)*ff(i,j,n3-2))*bp(n3-1)


                    do k = n3-2, 2, -1
                        ff(i,j,k) = (ff(i,j,k) - A3(k)*ff(i,j,k+1))*bp(k)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D0ssh_CPT_3Dz

    subroutine D0ssh_CPT_MULT_3Dz(f, ff, n1,n1e,n2,n2e,n3, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: ff

        real(kind=8), dimension(n3)         :: tmp
        real(kind=8), dimension(n1,n2,n3)   :: tmp2

        real(kind=8), dimension(n1,n2)                      :: sz

        integer :: i,j,k

        integer, save :: last_n=0, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h

        call fill_matrix(n3, .true., bc_type, 0.d0)

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e
                    tmp2(i,j,1)= B(2)*(f(i,j,2)+f(i,j,n3-2)) + B(1)*(f(i,j,1)+f(i,j,n3-1))
                    tmp2(i,j,2)= B(2)*(f(i,j,3)+f(i,j,n3-1)) + B(1)*(f(i,j,2)+f(i,j,1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    do k=3,n3-2
                        tmp2(i,j,k)= B(2)*(f(i,j,k+1)+f(i,j,k-2)) + B(1)*(f(i,j,k)+f(i,j,k-1))
                    enddo
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    tmp2(i,j,n3-1)= B(2)*(f(i,j,1)+f(i,j,n3-3)) + B(1)*(f(i,j,n3-1)+f(i,j,n3-2))
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
                        ff(i,j,k)=(tmp2(i,j,k)-sz(i,j)*q(k))*ff(i,j,k)
                    enddo
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e

                    tmp(2:n3-1)=ff(i,j,2:n3-1)

                    do k=3,n3-2
                        ff(i,j,k)= B(2)*(f(i,j,k+1)+f(i,j,k-2)) + B(1)*(f(i,j,k)+f(i,j,k-1))
                    enddo

                    ff(i,j,2)= ( Bl(-1)*f(i,j,1) + Bl(0) * f(i,j,2) + Bl(1) * f(i,j,3) )
                    ff(i,j,n3-1)= ( + Bl(-1)*f(i,j,n3-1) + Bl(0) * f(i,j,n3-2) + Bl(1)*f(i,j,n3-3) )


                    do k=3, n3-2
                        ff(i,j,k)= ff(i,j,k) - m(k)*ff(i,j,k-1)
                    enddo

                    ff(i,j,n3-1)= (ff(i,j,n3-1) - m(n3-1)*ff(i,j,n3-2))*bp(n3-1)


                    do k = n3-2, 2, -1
                        ff(i,j,k) = (ff(i,j,k) - A3(k)*ff(i,j,k+1))*bp(k)
                    enddo
                    do k=2, n3-1
                        ff(i,j,k)= ff(i,j,k)*tmp(k)
                    enddo

                enddo
            enddo

        endif

        return

    end subroutine D0ssh_CPT_MULT_3Dz

end module CPT_D0s_sh




