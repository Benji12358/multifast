
module d0s_schemes

    use boundaries_types
    use schemes_settings

    implicit none

contains

    subroutine D0s_ExpCtr_O2Fp0(f, d0f, n_max, h, shifted, bc_type)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d0f(:)

        integer :: i
        real(kind=8) A

        A=0.5d0

        if (shifted) then

            do i=2,n_max-1
                d0f(i)= ( f(i) + f(i-1) ) * A
            enddo

            if (bc_type.eq.periodic) then
                d0f(1)= ( f(1) + f(n_max-1) ) * A
            end if


        else

            do i=1,n_max-2
                d0f(i)= ( f(i+1) + f(i) ) * A
            enddo

            if (bc_type.eq.periodic) then
                d0f(n_max-1)= ( f(1) + f(n_max-1) ) * A
            end if

            if (bc_type.eq.Dirichlet) then
                d0f(n_max-1)= ( f(n_max) + f(n_max-1) ) * A
            end if

        end if

        return

    end subroutine

    subroutine D0s_ExpCtr_O2Fp0_MULT(f, d0f, n_max, h, shifted, bc_type)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(inout)   :: d0f(:)

        integer :: i
        real(kind=8) A

        A=0.5d0

        if (shifted) then

            do i=2,n_max-1
                d0f(i)= ( f(i) + f(i-1) ) * A * d0f(i)
            enddo

            if (bc_type.eq.periodic) then
                d0f(1)= ( f(1) + f(n_max-1) ) * A * d0f(1)
            end if


        else

            do i=1,n_max-2
                d0f(i)= ( f(i+1) + f(i) ) * A * d0f(i)
            enddo

            if (bc_type.eq.periodic) then
                d0f(n_max-1)= ( f(1) + f(n_max-1) ) * A * d0f(n_max-1)
            end if

            if (bc_type.eq.Dirichlet) then
                d0f(n_max-1)= ( f(n_max) + f(n_max-1) ) * A * d0f(n_max-1)
            end if

        end if

        return

    end subroutine

    subroutine D0s_ExpCtr_O0Fp5(f, d0f, n_max, h, shifted, bc_type)
        ! This schemes has been designed in the file
        ! FILTER_EXP_6pts.wxm. The different cases taken into
        ! account here are for different values of cutting freq. wc
        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d0f(:)

        integer :: i
        real(kind=8) A(5)

        A(1) =  1.234324263495304d0     /2.d0
        A(2) =  -0.31902734931011d0     /2.d0
        A(3) =  0.11097825729549d0      /2.d0
        A(4) =  -0.031121869885979d0    /2.d0
        A(5) =  0.0048483383805266d0    /2.d0

        if (shifted) then

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=6,n_max-5
                d0f(i)=   A(5)*(f(i+4) + f(i-5))    &
                + A(4)*(f(i+3) + f(i-4))            &
                + A(3)*(f(i+2) + f(i-3))            &
                + A(2)*(f(i+1) + f(i-2))            &
                + A(1)*(f(i) + f(i-1))
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 5
                    d0f(i)=   A(5)*(f(i+4) + f(mod(n_max+i-7,n_max-1)+1))   &
                    + A(4)*(f(i+3) + f(mod(n_max+i-6,n_max-1)+1))           &
                    + A(3)*(f(i+2) + f(mod(n_max+i-5,n_max-1)+1))           &
                    + A(2)*(f(i+1) + f(mod(n_max+i-4,n_max-1)+1))           &
                    + A(1)*(f(i) + f(mod(n_max+i-3,n_max-1)+1))
                enddo

                do i=n_max-4,n_max-1
                    d0f(i)=   A(5)*(f(mod(i+3,n_max-1)+1) + f(i-5))     &
                    + A(4)*(f(mod(i+2,n_max-1)+1) + f(i-4))             &
                    + A(3)*(f(mod(i+1,n_max-1)+1) + f(i-3))             &
                    + A(2)*(f(mod(i,n_max-1)+1) + f(i-2))             &
                    + A(1)*(f(i) + f(i-1))
                enddo
            end if

            if (bc_type.eq.Dirichlet) then
                call D0s_CTR_CLOSURE_INF(f, d0f, h, 5, shifted)
            end if

            if (bc_type.eq.Dirichlet) then
                call D0s_CTR_CLOSURE_SUP(f, d0f, n_max, h, 5, shifted)
            end if


        else

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=5,n_max-6
                d0f(i)=   A(5)*(f(i+5) + f(i-4))    &
                + A(4)*(f(i+4) + f(i-3))            &
                + A(3)*(f(i+3) + f(i-2))            &
                + A(2)*(f(i+2) + f(i-1))            &
                + A(1)*(f(i+1) + f(i))
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 4
                    d0f(i)=   A(5)*(f(i+5) + f(mod(n_max+i-6,n_max-1)+1))   &
                    + A(4)*(f(i+4) + f(mod(n_max+i-5,n_max-1)+1))           &
                    + A(3)*(f(i+3) + f(mod(n_max+i-4,n_max-1)+1))           &
                    + A(2)*(f(i+2) + f(mod(n_max+i-3,n_max-1)+1))           &
                    + A(1)*(f(i+1) + f(i))
                enddo

                do i=n_max-5,n_max-1
                    d0f(i)=   A(5)*(f(mod(i+4,n_max-1)+1) + f(i-4))     &
                    + A(4)*(f(mod(i+3,n_max-1)+1) + f(i-3))             &
                    + A(3)*(f(mod(i+2,n_max-1)+1) + f(i-2))             &
                    + A(2)*(f(mod(i+1,n_max-1)+1) + f(i-1))             &
                    + A(1)*(f(mod(i,n_max-1)+1) + f(i))
                enddo
            end if

            if (bc_type.eq.Dirichlet) then
                call D0s_CTR_CLOSURE_INF(f, d0f, h, 4, shifted)
            end if

            if (bc_type.eq.Dirichlet) then
                call D0s_CTR_CLOSURE_SUP(f, d0f, n_max, h, 5, shifted)
            end if

        end if

        return

    end subroutine

    subroutine D0s_ExpCtr_O0Fp5_MULT(f, d0f, n_max, h, shifted, bc_type)
        ! This schemes has been designed in the file
        ! FILTER_EXP_6pts.wxm. The different cases taken into
        ! account here are for different values of cutting freq. wc
        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(inout)   :: d0f(:)

        integer :: i
        real(kind=8) A(5)

        A(1) =  1.234324263495304d0     /2.d0
        A(2) =  -0.31902734931011d0     /2.d0
        A(3) =  0.11097825729549d0      /2.d0
        A(4) =  -0.031121869885979d0    /2.d0
        A(5) =  0.0048483383805266d0    /2.d0

        if (shifted) then

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=6,n_max-5
                d0f(i)= ( A(5)*(f(i+4) + f(i-5))    &
                + A(4)*(f(i+3) + f(i-4))            &
                + A(3)*(f(i+2) + f(i-3))            &
                + A(2)*(f(i+1) + f(i-2))            &
                + A(1)*(f(i) + f(i-1)) )*d0f(i)
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 5
                    d0f(i)= ( A(5)*(f(i+4) + f(mod(n_max+i-7,n_max-1)+1))   &
                    + A(4)*(f(i+3) + f(mod(n_max+i-6,n_max-1)+1))           &
                    + A(3)*(f(i+2) + f(mod(n_max+i-5,n_max-1)+1))           &
                    + A(2)*(f(i+1) + f(mod(n_max+i-4,n_max-1)+1))           &
                    + A(1)*(f(i) + f(mod(n_max+i-3,n_max-1)+1)) )*d0f(i)
                enddo

                do i=n_max-4,n_max-1
                    d0f(i)=   ( A(5)*(f(mod(i+3,n_max-1)+1) + f(i-5))       &
                    + A(4)*(f(mod(i+2,n_max-1)+1) + f(i-4))                 &
                    + A(3)*(f(mod(i+1,n_max-1)+1) + f(i-3))                 &
                    + A(2)*(f(mod(i,n_max-1)+1) + f(i-2))                   &
                    + A(1)*(f(i) + f(i-1)) )*d0f(i)
                enddo
            end if

            if (bc_type.eq.Dirichlet) then
                call D0s_CTR_CLOSURE_INF_MULT(f, d0f, h, 5, shifted)
            end if

            if (bc_type.eq.Dirichlet) then
                call D0s_CTR_CLOSURE_SUP_MULT(f, d0f, n_max, h, 5, shifted)
            end if


        else

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=5,n_max-6
                d0f(i)= ( A(5)*(f(i+5) + f(i-4))    &
                + A(4)*(f(i+4) + f(i-3))            &
                + A(3)*(f(i+3) + f(i-2))            &
                + A(2)*(f(i+2) + f(i-1))            &
                + A(1)*(f(i+1) + f(i)) )*d0f(i)
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 4
                    d0f(i)= ( A(5)*(f(i+5) + f(mod(n_max+i-6,n_max-1)+1))   &
                    + A(4)*(f(i+4) + f(mod(n_max+i-5,n_max-1)+1))           &
                    + A(3)*(f(i+3) + f(mod(n_max+i-4,n_max-1)+1))           &
                    + A(2)*(f(i+2) + f(mod(n_max+i-3,n_max-1)+1))           &
                    + A(1)*(f(i+1) + f(i)) )*d0f(i)
                enddo

                do i=n_max-5,n_max-1
                    d0f(i)= ( A(5)*(f(mod(i+4,n_max-1)+1) + f(i-4))     &
                    + A(4)*(f(mod(i+3,n_max-1)+1) + f(i-3))             &
                    + A(3)*(f(mod(i+2,n_max-1)+1) + f(i-2))             &
                    + A(2)*(f(mod(i+1,n_max-1)+1) + f(i-1))             &
                    + A(1)*(f(mod(i,n_max-1)+1) + f(i)) )*d0f(i)
                enddo
            end if

            if (bc_type.eq.Dirichlet) then
                call D0s_CTR_CLOSURE_INF_MULT(f, d0f, h, 4, shifted)
            end if

            if (bc_type.eq.Dirichlet) then
                call D0s_CTR_CLOSURE_SUP_MULT(f, d0f, n_max, h, 5, shifted)
            end if

        end if

        return

    end subroutine

    subroutine D0s_CptCtr_O6Fp0(f, d0f, n_max, h, shifted, bc_type)

        ! f     : [1..n_max]
        ! d0f   : [1..n_max+1]
        ! fs ---d0f(1)---f(1)---d0f(2)--...f(nmax)--d0f(nmax)---fn
        use Thomas_solver

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d0f(:)

        integer :: i
        real(kind=8)  :: al

        real(kind=8)  ::rhs(n_max)

        real(kind=8), dimension(:,:), allocatable, save :: A
        real(kind=8), dimension(:), save, allocatable   :: u, v
        real(kind=8), save                              :: vq
        real*8                                          :: coef_q, vy
        real(kind=8), save    :: B(-2:2), B_lim(-1:3)

        integer, save :: last_nmax, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h



        update_matrix=(last_nmax/=n_max).or.(last_bc_type/=bc_type)    &
        .or.(last_shifted.neqv.shifted).or.(last_h/=h)

        !if (update_matrix) then
        call fill_matrix
        !end if

        last_nmax=n_max
        last_bc_type=bc_type
        last_shifted=shifted
        last_h=h

        if (shifted) then

            do i=3,n_max-2
                rhs(i)= ( B(2)*f(i+1) + B(1)*f(i) + B(-1)*f(i-1) + B(-2)*f(i-2) )
            enddo

            select case (bc_type)

                case (periodic)

                    rhs(1)= ( B(2)*f(2) + B(1)*f(1) + B(-1)*f(n_max-1) + B(-2)*f(n_max-2) )
                    rhs(2)= ( B(2)*f(3) + B(1)*f(2) + B(-1)*f(1) + B(-2)*f(n_max-1) )

                    rhs(n_max-1)= ( B(2)*f(1) + B(1)*f(n_max-1) + B(-1)*f(n_max-2) + B(-2)*f(n_max-3) )

                    call TS_NPr(rhs, d0f, n_max-1)

                    vy=0.d0
                    do i=1, n_max-1
                        vy=vy+v(i)*d0f(i)
                    enddo
                    coef_q = vy/(1+vq)

                    do i=1, n_max-1
                        d0f(i)=d0f(i)-coef_q*u(i)
                    enddo

                case (Dirichlet)

                    rhs(2)= ( B_lim(-1)*f(1) + B_lim(0) * f(2) + B_lim(1) * f(3) )

                    rhs(n_max-1)= ( + B_lim(-1)*f(n_max-1) + B_lim(0) * f(n_max-2) + B_lim(1)*f(n_max-3) )

                    call TS_NPr(rhs(2:n_max-1), d0f(2:n_max-1), n_max-2)


            end select

        else

            do i=2,n_max-3
                rhs(i)= ( B(2)*f(i+2) + B(1)*f(i+1) + B(-1)*f(i) + B(-2)*f(i-1) )
            enddo

            select case (bc_type)

                case (periodic)

                    rhs(1)= ( B(2)*f(3) + B(1)*f(2) + B(-1)*f(1) + B(-2)*f(n_max-1) )

                    rhs(n_max-1)= ( B(2)*f(2) + B(1)*f(1) + B(-1)*f(n_max-1) + B(-2)*f(n_max-2) )
                    rhs(n_max-2)= ( B(2)*f(1) + B(1)*f(n_max-1) + B(-1)*f(n_max-2) + B(-2)*f(n_max-3) )

                    call TS_Pr_opt(A(-1,1:n_max-1), A(0,1:n_max-1), A(1,1:n_max-1), rhs(1:n_max-1), d0f(1:n_max-1), u,v,vq,n_max-1)

                case (Dirichlet)

                    rhs(1)= ( B_lim(-1)*f(1) + B_lim(0) * f(2) + B_lim(1) * f(3) )

                    rhs(n_max-2)= ( B(2)*f(n_max) + B(1)*f(n_max-1) + B(-1)*f(n_max-2) + B(-2)*f(n_max-3) )
                    rhs(n_max-1)= ( + B_lim(-1)*f(n_max) + B_lim(0) * f(n_max-1) + B_lim(1)*f(n_max-2) )

                    call TS_NPr(rhs(1:n_max-1), d0f(1:n_max-1), n_max-1)

            end select

        end if

        return

    contains


        subroutine fill_matrix()

            if (allocated(A)) then
                deallocate(A)
                deallocate(u)
                deallocate(v)
            end if

            allocate(A(-1:1, n_max))
            allocate(u(n_max), v(n_max))

            al=3.00001d0/10.d0
            B(1)=(10.d0*al+9.d0)/16.d0
            B(2)= (6.d0*al-1)/16.d0
            B(-1)=B(1)
            B(-2)=B(2)

            B_lim(-1)=1.d0/4.d0
            B_lim(0)= 3.d0/2.d0
            B_lim(1)=1.d0/4.d0

            if (shifted) then

                do i= 3, n_max-2
                    A(-1,i) = al
                    A(0,i) = 1.d0
                    A(1,i) = al
                enddo

                if (bc_type.eq.periodic) then

                    A(-1,n_max-2:n_max-1) = al
                    A(0,n_max-2:n_max-1) = 1.d0
                    A(1,n_max-2:n_max-1) = al

                    A(-1,1:2) = al
                    A(0,1:2) = 1.d0
                    A(1,1:2) = al

                    call TS_init_Pr(A(-1,1:n_max-1), A(0,1:n_max-1), A(1,1:n_max-1), u, v, vq,n_max-1)

                end if

                if (bc_type.eq.Dirichlet) then
                    A(-1,2)=0.d0
                    A(0,2)=1.d0
                    A(1,2)=1.d0

                    A(-1,n_max-1)=1.d0
                    A(0,n_max-1)=1.d0
                    A(1,n_max-1)=0.d0

                    call TS_init(A(-1,2:n_max-1), A(0,2:n_max-1), A(1,2:n_max-1), n_max-2)

                end if

            else

                do i= 2, n_max-2
                    A(-1,i) = al
                    A(0,i) = 1.d0
                    A(1,i) = al
                enddo

                if (bc_type.eq.periodic) then

                    A(-1,n_max-1) = al
                    A(0,n_max-1) = 1.d0
                    A(1,n_max-1) = al

                    A(-1,1) = al
                    A(0,1) = 1.d0
                    A(1,1) = al

                    call TS_init_Pr(A(-1,1:n_max-1), A(0,1:n_max-1), A(1,1:n_max-1), u, v, vq,n_max-1)

                end if

                if (bc_type.eq.Dirichlet) then
                    A(-1,n_max-1)=1.d0
                    A(0,n_max-1)=1.d0

                    A(0,1)=1.d0
                    A(1,1)=1.d0

                    call TS_init(A(-1,1:n_max-1), A(0,1:n_max-1), A(1,1:n_max-1), n_max-1)

                end if

            end if

        end subroutine fill_matrix


    end subroutine

    subroutine D0s_CptCtr_O6Fp0_MULT(f, d0f, n_max, h, shifted, bc_type)

        implicit none
        logical, intent(in) :: shifted
        integer, intent(in) :: n_max, bc_type
        real(kind=8) , intent(in) :: h, f(:)
        real(kind=8), intent(inout) :: d0f(:)

        real(kind=8), dimension(n_max) :: d0f_tmp
        integer :: i

        select case (bc_type)

            case (periodic)

                d0f_tmp(1:n_max-1)=d0f(1:n_max-1)
                call D0s_CptCtr_O6Fp0(f, d0f_tmp, n_max, h, shifted, bc_type)
                d0f(1:n_max-1) = d0f(1:n_max-1) * d0f_tmp(1:n_max-1)

            case (Dirichlet)

                if (shifted) then

                    d0f_tmp(2:n_max-1)=d0f(2:n_max-1)
                    call D0s_CptCtr_O6Fp0(f, d0f_tmp, n_max, h, shifted, bc_type)
                    d0f(2:n_max-1) = d0f(2:n_max-1) * d0f_tmp(2:n_max-1)

                else

                    d0f_tmp(1:n_max-1)=d0f(1:n_max-1)
                    call D0s_CptCtr_O6Fp0(f, d0f_tmp, n_max, h, shifted, bc_type)
                    d0f(1:n_max-1) = d0f(1:n_max-1) * d0f_tmp(1:n_max-1)
                end if

            case default

        end select

        return

    end subroutine

    subroutine D0s_CTR_CLOSURE_INF(f, d0f, h, nb_ghost, shifted)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(nb_ghost)

        implicit none

        integer :: i, n_max, nb_ghost
        logical ::  shifted
        real(kind=8) h, fs
        real(kind=8) d0f(:)
        real(kind=8) f(:)

        real(kind=8) d0f_temp(CPT_MIN)

        if (nb_ghost.gt.5) then
            call exit(1)
        end if

        if (shifted) then


            if (treat_boundary_by_cpt_scheme) then

                call D0s_CptCtr_O6Fp0(f(1:CPT_MIN), d0f_temp(1:CPT_MIN), CPT_MIN, 0.d0, .true., Dirichlet)

                d0f(2)=d0f_temp(2)
                if (nb_ghost.eq.2) return

                d0f(3)=d0f_temp(3)
                if (nb_ghost.eq.3) return

            else

                ! i=2---------------------------------------------------
                i=2
                d0f(i)=0.5d0*(f(i)+f(i-1))
                ! i=3---------------------------------------------------
                i=3
                d0f(i)  =   9.0d0/16.d0  * (f(i)+f(i-1))        &
                -   0.0625d0  * (f(i+1)+f(i-2))


                if (nb_ghost.eq.i) return

            endif

            ! i=4---------------------------------------------------
            i=4
            d0f(i)  =   0.59395104312381d0  *   (f(i)+f(i-1))       &
            -   0.10967656468571d0  *   (f(i+1)+f(i-2))             &
            +   0.015725521561903d0 *   (f(i+2) + f(i-3))

            if (nb_ghost.eq.i) return

            ! i=5---------------------------------------------------
            i=5
            d0f(i)  =   0.59395104312381d0  *   (f(i)+f(i-1))       &
            -   0.10967656468571d0  *   (f(i+1)+f(i-2))             &
            +   0.015725521561903d0 *   (f(i+2) + f(i-3))

            if (nb_ghost.eq.i) return

        else


            if (treat_boundary_by_cpt_scheme) then

                call D0s_CptCtr_O6Fp0(f(1:CPT_MIN), d0f_temp(1:CPT_MIN), CPT_MIN, 0.d0, .false., Dirichlet)

                d0f(1)=d0f_temp(1)
                if (nb_ghost.eq.1) return

                d0f(2)=d0f_temp(2)
                if (nb_ghost.eq.2) return

            else

                ! i=1---------------------------------------------------
                i=1
                d0f(i)=0.5d0*(f(i+1)+f(i))
                ! i=2---------------------------------------------------
                i=2
                d0f(i)  =   9.0d0/16.d0  * (f(i+1)+f(i))        &
                -   0.0625d0  * (f(i+2)+f(i-1))


                if (nb_ghost.eq.i) return

            endif

            ! i=3---------------------------------------------------
            i=3
            d0f(i)  =   0.59395104312381d0  *   (f(i+1)+f(i))       &
            -   0.10967656468571d0  *   (f(i+2)+f(i-1))             &
            +   0.015725521561903d0 *   (f(i+3) + f(i-2))

            if (nb_ghost.eq.i) return

            ! i=4---------------------------------------------------
            i=4
            d0f(i)  =   0.59395104312381d0  *   (f(i+1)+f(i))       &
            -   0.10967656468571d0  *   (f(i+2)+f(i-1))             &
            +   0.015725521561903d0 *   (f(i+3) + f(i-2))

            if (nb_ghost.eq.i) return

            ! i=5---SEE D0s_5pts.wxm-------wc=1.8-------------------
            i=5
            d0f(i)  =   0.61716213174765d0  *   (f(i+1)+f(i))       &
            -   0.15951367465505d0  *   (f(i+2)+f(i-1))     &
            +   0.055489128647743d0 *   (f(i+3) + f(i-2))   &
            -   0.015560934942989d0 *   (f(i+4) + f(i-3))   &
            +   0.0024241691902633d0*   (f(i+5) + f(i-4))

        end if

        return

    end subroutine

    subroutine D0s_CTR_CLOSURE_INF_MULT(f, d0f, h, nb_ghost, shifted)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(nb_ghost)

        implicit none

        integer :: i, n_max, nb_ghost
        logical ::  shifted
        real(kind=8) h, fs
        real(kind=8) d0f(:)
        real(kind=8) f(:)

        real(kind=8) d0f_temp(CPT_MIN)

        if (nb_ghost.gt.5) then
            call exit(1)
        end if

        if (shifted) then


            if (treat_boundary_by_cpt_scheme) then

                call D0s_CptCtr_O6Fp0(f(1:CPT_MIN), d0f_temp(1:CPT_MIN), CPT_MIN, 0.d0, .true., Dirichlet)

                d0f(2)=d0f_temp(2)*d0f(2)
                if (nb_ghost.eq.2) return

                d0f(3)=d0f_temp(3)*d0f(3)
                if (nb_ghost.eq.3) return

            else

                ! i=2---------------------------------------------------
                i=2
                d0f(i)=( 0.5d0*(f(i)+f(i-1)) )*d0f(i)
                ! i=3---------------------------------------------------
                i=3
                d0f(i)  = ( 9.0d0/16.d0  * (f(i)+f(i-1))        &
                -   0.0625d0  * (f(i+1)+f(i-2)) )*d0f(i)


                if (nb_ghost.eq.i) return

            endif

            ! i=4---------------------------------------------------
            i=4
            d0f(i)  = ( 0.59395104312381d0  *   (f(i)+f(i-1))       &
            -   0.10967656468571d0  *   (f(i+1)+f(i-2))             &
            +   0.015725521561903d0 *   (f(i+2) + f(i-3)) )*d0f(i)

            if (nb_ghost.eq.i) return

            ! i=5---------------------------------------------------
            i=5
            d0f(i)  = ( 0.59395104312381d0  *   (f(i)+f(i-1))       &
            -   0.10967656468571d0  *   (f(i+1)+f(i-2))             &
            +   0.015725521561903d0 *   (f(i+2) + f(i-3)) )*d0f(i)

            if (nb_ghost.eq.i) return

        else


            if (treat_boundary_by_cpt_scheme) then

                call D0s_CptCtr_O6Fp0(f(1:CPT_MIN), d0f_temp(1:CPT_MIN), CPT_MIN, 0.d0, .false., Dirichlet)

                d0f(1)=d0f_temp(1)*d0f(1)
                if (nb_ghost.eq.1) return

                d0f(2)=d0f_temp(2)*d0f(2)
                if (nb_ghost.eq.2) return

            else

                ! i=1---------------------------------------------------
                i=1
                d0f(i)=( 0.5d0*(f(i+1)+f(i)) )*d0f(i)
                ! i=2---------------------------------------------------
                i=2
                d0f(i)  = ( 9.0d0/16.d0  * (f(i+1)+f(i))        &
                -   0.0625d0  * (f(i+2)+f(i-1)) )*d0f(i)


                if (nb_ghost.eq.i) return

            endif

            ! i=3---------------------------------------------------
            i=3
            d0f(i)  = ( 0.59395104312381d0  *   (f(i+1)+f(i))       &
            -   0.10967656468571d0  *   (f(i+2)+f(i-1))             &
            +   0.015725521561903d0 *   (f(i+3) + f(i-2)) )*d0f(i)

            if (nb_ghost.eq.i) return

            ! i=4---------------------------------------------------
            i=4
            d0f(i)  = ( 0.59395104312381d0  *   (f(i+1)+f(i))       &
            -   0.10967656468571d0  *   (f(i+2)+f(i-1))             &
            +   0.015725521561903d0 *   (f(i+3) + f(i-2)) )*d0f(i)

            if (nb_ghost.eq.i) return

            ! i=5---SEE D0s_5pts.wxm-------wc=1.8-------------------
            i=5
            d0f(i)  = ( 0.61716213174765d0  *   (f(i+1)+f(i))       &
            -   0.15951367465505d0  *   (f(i+2)+f(i-1))     &
            +   0.055489128647743d0 *   (f(i+3) + f(i-2))   &
            -   0.015560934942989d0 *   (f(i+4) + f(i-3))   &
            +   0.0024241691902633d0*   (f(i+5) + f(i-4)) )*d0f(i)

        end if

        return

    end subroutine

    subroutine D0s_CTR_CLOSURE_SUP(f, d0f, n_max, h, nb_ghost, shifted)

        implicit none

        integer :: i, n_max, nb_ghost, j
        real(kind=8) h, fn
        real(kind=8) d0f(:)
        real(kind=8) f(:)
        logical :: shifted

        real(kind=8) d0f_temp(n_max-CPT_MIN:n_max-1)

        if (nb_ghost.gt.5) then
            call exit(1)
        end if

        if (shifted) then

            if (treat_boundary_by_cpt_scheme) then

                !   ---f------d0f-----f------d0f-------f
                call D0s_CptCtr_O6Fp0(f(n_max-CPT_MIN-1:n_max-1), d0f_temp, CPT_MIN+1, h, .false., Dirichlet)

                d0f(n_max-1)=d0f_temp(n_max-1)

                d0f(n_max-2)=d0f_temp(n_max-2)
                if (nb_ghost.eq.2) return

            else

                ! First point-------------------------------------------
                i=n_max-1
                d0f(i)=0.5d0*(f(i)+f(i-1))
                ! 2nd point----------------------------------------------
                i=n_max-2
                d0f(i)  =   9.0d0/16.d0 * (f(i)+f(i-1))     &
                -   0.0625d0    * (f(i+1)+f(i-2))


                if (nb_ghost.eq.2) return

            endif

            ! 3rd point---------------------------------------------
            i=n_max-3
            d0f(i)  =   0.59395104312381d0  *   (f(i)+f(i-1))       &
            -   0.10967656468571d0  *   (f(i+1)+f(i-2))             &
            +   0.015725521561903d0 *   (f(i+2) + f(i-3))

            if (nb_ghost.eq.3) return

            ! 4th point---------------------------------------------
            i=n_max-4
            d0f(i)  =   0.59395104312381d0  *   (f(i)+f(i-1))       &
            -   0.10967656468571d0  *   (f(i+1)+f(i-2))             &
            +   0.015725521561903d0 *   (f(i+2) + f(i-3))

            if (nb_ghost.eq.4) return

        else

            if (treat_boundary_by_cpt_scheme) then

                !   ---f------d0f-----f------d0f-------f
                call D0s_CptCtr_O6Fp0(f(n_max-CPT_MIN:n_max), d0f_temp, CPT_MIN+1, h, .false., Dirichlet)

                d0f(n_max-1)=d0f_temp(n_max-1)

                d0f(n_max-2)=d0f_temp(n_max-2)
                if (nb_ghost.eq.2) return

            else

                ! First point-------------------------------------------
                i=n_max-1
                d0f(i)=0.5d0*(f(i+1)+f(i))
                ! 2nd point----------------------------------------------
                i=n_max-2
                d0f(i)  =   9.0d0/16.d0 * (f(i+1)+f(i))     &
                -   0.0625d0    * (f(i+2)+f(i-1))


                if (nb_ghost.eq.2) return

            endif

            ! 3rd point---------------------------------------------
            i=n_max-3
            d0f(i)  =   0.59395104312381d0  *   (f(i+1)+f(i))       &
            -   0.10967656468571d0  *   (f(i+2)+f(i-1))             &
            +   0.015725521561903d0 *   (f(i+3) + f(i-2))

            if (nb_ghost.eq.3) return

            ! 4th point---------------------------------------------
            i=n_max-4
            d0f(i)  =   0.59395104312381d0  *   (f(i+1)+f(i))       &
            -   0.10967656468571d0  *   (f(i+2)+f(i-1))             &
            +   0.015725521561903d0 *   (f(i+3) + f(i-2))

            if (nb_ghost.eq.4) return

            ! 5th point ----D0s_5pts.wxm-------wc=1.8---------------
            i=n_max-5
            d0f(i)  =   0.61716213174765d0  *   (f(i+1)+f(i))       &
            -   0.15951367465505d0  *   (f(i+2)+f(i-1))             &
            +   0.055489128647743d0 *   (f(i+3) + f(i-2))           &
            -   0.015560934942989d0 *   (f(i+4) + f(i-3))           &
            +   0.0024241691902633d0*   (f(i+5) + f(i-4))

        end if

        return

    end subroutine

    subroutine D0s_CTR_CLOSURE_SUP_MULT(f, d0f, n_max, h, nb_ghost, shifted)

        implicit none

        integer :: i, n_max, nb_ghost, j
        real(kind=8) h, fn
        real(kind=8) d0f(:)
        real(kind=8) f(:)
        logical :: shifted

        real(kind=8) d0f_temp(n_max-CPT_MIN:n_max-1)

        if (nb_ghost.gt.5) then
            call exit(1)
        end if

        if (shifted) then

            if (treat_boundary_by_cpt_scheme) then

                !   ---f------d0f-----f------d0f-------f
                call D0s_CptCtr_O6Fp0(f(n_max-CPT_MIN-1:n_max-1), d0f_temp, CPT_MIN+1, h, .false., Dirichlet)

                d0f(n_max-1)=d0f_temp(n_max-1)*d0f(n_max-1)

                d0f(n_max-2)=d0f_temp(n_max-2)*d0f(n_max-2)
                if (nb_ghost.eq.2) return

            else

                ! First point-------------------------------------------
                i=n_max-1
                d0f(i)= ( 0.5d0*(f(i)+f(i-1)) )*d0f(i)
                ! 2nd point----------------------------------------------
                i=n_max-2
                d0f(i)  =  ( 9.0d0/16.d0 * (f(i)+f(i-1))     &
                -   0.0625d0    * (f(i+1)+f(i-2)) )*d0f(i)


                if (nb_ghost.eq.2) return

            endif

            ! 3rd point---------------------------------------------
            i=n_max-3
            d0f(i)  = ( 0.59395104312381d0  *   (f(i)+f(i-1))       &
            -   0.10967656468571d0  *   (f(i+1)+f(i-2))             &
            +   0.015725521561903d0 *   (f(i+2) + f(i-3)) )*d0f(i)

            if (nb_ghost.eq.3) return

            ! 4th point---------------------------------------------
            i=n_max-4
            d0f(i)  = ( 0.59395104312381d0  *   (f(i)+f(i-1))       &
            -   0.10967656468571d0  *   (f(i+1)+f(i-2))             &
            +   0.015725521561903d0 *   (f(i+2) + f(i-3)) )*d0f(i)

            if (nb_ghost.eq.4) return

        else

            if (treat_boundary_by_cpt_scheme) then

                !   ---f------d0f-----f------d0f-------f
                call D0s_CptCtr_O6Fp0(f(n_max-CPT_MIN:n_max), d0f_temp, CPT_MIN+1, h, .false., Dirichlet)

                d0f(n_max-1)=d0f_temp(n_max-1)*d0f(n_max-1)

                d0f(n_max-2)=d0f_temp(n_max-2)*d0f(n_max-2)
                if (nb_ghost.eq.2) return

            else

                ! First point-------------------------------------------
                i=n_max-1
                d0f(i)= ( 0.5d0*(f(i+1)+f(i)) )*d0f(i)
                ! 2nd point----------------------------------------------
                i=n_max-2
                d0f(i)  = ( 9.0d0/16.d0 * (f(i+1)+f(i))     &
                -   0.0625d0    * (f(i+2)+f(i-1)) )*d0f(i)


                if (nb_ghost.eq.2) return

            endif

            ! 3rd point---------------------------------------------
            i=n_max-3
            d0f(i)  = ( 0.59395104312381d0  *   (f(i+1)+f(i))       &
            -   0.10967656468571d0  *   (f(i+2)+f(i-1))             &
            +   0.015725521561903d0 *   (f(i+3) + f(i-2)) )*d0f(i)

            if (nb_ghost.eq.3) return

            ! 4th point---------------------------------------------
            i=n_max-4
            d0f(i)  = ( 0.59395104312381d0  *   (f(i+1)+f(i))       &
            -   0.10967656468571d0  *   (f(i+2)+f(i-1))             &
            +   0.015725521561903d0 *   (f(i+3) + f(i-2)) )*d0f(i)

            if (nb_ghost.eq.4) return

            ! 5th point ----D0s_5pts.wxm-------wc=1.8---------------
            i=n_max-5
            d0f(i)  = ( 0.61716213174765d0  *   (f(i+1)+f(i))       &
            -   0.15951367465505d0  *   (f(i+2)+f(i-1))             &
            +   0.055489128647743d0 *   (f(i+3) + f(i-2))           &
            -   0.015560934942989d0 *   (f(i+4) + f(i-3))           &
            +   0.0024241691902633d0*   (f(i+5) + f(i-4)) )*d0f(i)

        end if

        return

    end subroutine

end module

module d1s_schemes
    use boundaries_types
    use schemes_settings

    implicit none

contains

    subroutine D1s_ExpCtr_O2Fp0(f, d1f, n_max, h, shifted, bc_type)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)   :: d1f(:)

        integer :: i
        real(kind=8) A

        A=1.d0/h

        if (shifted) then

            do i=2,n_max-1
                d1f(i)= ( f(i) - f(i-1) ) * A
            enddo

            if (bc_type.eq.periodic) then
                d1f(1)= ( f(1) - f(n_max-1) ) * A
            end if

        else

            do i=1,n_max-2
                d1f(i)= ( f(i+1) - f(i) ) * A
            enddo

            if (bc_type.eq.periodic) then
                d1f(n_max-1)= ( f(1) - f(n_max-1) ) * A
            end if

            if ((bc_type==Dirichlet).or.(bc_type==antisymetric)) then
                d1f(n_max-1)= ( f(n_max) - f(n_max-1) ) * A
            end if

        end if

        return

    end subroutine

    subroutine D1s_ExpCtr_O2Fp0_ACC(f, d1f, n_max, h, shifted, bc_type)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(inout)   :: d1f(:)

        integer :: i
        real(kind=8) A

        A=1.d0/h

        if (shifted) then

            do i=2,n_max-1
                d1f(i)= d1f(i) + ( f(i) - f(i-1) ) * A
            enddo

            if (bc_type.eq.periodic) then
                d1f(1)= d1f(1) + ( f(1) - f(n_max-1) ) * A
            end if

        else

            do i=1,n_max-2
                d1f(i)= d1f(i) + ( f(i+1) - f(i) ) * A
            enddo

            if (bc_type.eq.periodic) then
                d1f(n_max-1)= d1f(n_max-1) + ( f(1) - f(n_max-1) ) * A
            end if

            if ((bc_type==Dirichlet).or.(bc_type==antisymetric)) then
                d1f(n_max-1)= d1f(n_max-1) + ( f(n_max) - f(n_max-1) ) * A
            end if

        end if

        return

    end subroutine

    subroutine D1s_ExpCtr_O2Fp0_MULT(f, d1f, n_max, h, shifted, bc_type, geom_coef)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, geom_coef(:), f(:)
        real(kind=8), intent(out)   :: d1f(:)

        integer :: i
        real(kind=8) A

        A=1.d0/h

        if (shifted) then

            do i=2,n_max-1
                d1f(i)= ( f(i) - f(i-1) ) * A * geom_coef(i)
            enddo

            if (bc_type.eq.periodic) then
                d1f(1)= ( f(1) - f(n_max-1) ) * A * geom_coef(1)
            end if

        else

            do i=1,n_max-2
                d1f(i)= ( f(i+1) - f(i) ) * A * geom_coef(i)
            enddo

            if (bc_type.eq.periodic) then
                d1f(n_max-1)= ( f(1) - f(n_max-1) ) * A * geom_coef(n_max-1)
            end if

            if ((bc_type==Dirichlet).or.(bc_type==antisymetric)) then
                d1f(n_max-1)= ( f(n_max) - f(n_max-1) ) * A * geom_coef(n_max-1)
            end if

        end if

        return

    end subroutine

    subroutine D1s_ExpCtr_O2Fp0_MULT_ACC(f, d1f, n_max, h, shifted, bc_type, geom_coef)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, geom_coef(:), f(:)
        real(kind=8), intent(inout)   :: d1f(:)

        integer :: i
        real(kind=8) A

        A=1.d0/h

        if (shifted) then

            do i=2,n_max-1
                d1f(i)= d1f(i) + ( f(i) - f(i-1) ) * A * geom_coef(i)
            enddo

            if (bc_type.eq.periodic) then
                d1f(1)= d1f(1) + ( f(1) - f(n_max-1) ) * A * geom_coef(1)
            end if

        else

            do i=1,n_max-2
                d1f(i)= d1f(i) + ( f(i+1) - f(i) ) * A * geom_coef(i)
            enddo

            if (bc_type.eq.periodic) then
                d1f(n_max-1)= d1f(n_max-1) + ( f(1) - f(n_max-1) ) * A * geom_coef(n_max-1)
            end if

            if ((bc_type==Dirichlet).or.(bc_type==antisymetric)) then
                d1f(n_max-1)= d1f(n_max-1) + ( f(n_max) - f(n_max-1) ) * A * geom_coef(n_max-1)
            end if

        end if

        return

    end subroutine

    subroutine D1s_Tamm(f, d1f, n_max, h, shifted, bc_type)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)   :: d1f(:)

        integer :: i
        real(kind=8) A(3)

        A(1) =  1.189101951200031d0      /   h
        A(2) =  -0.073717642266682d0     /   h
        A(3) =  0.006410195120003d0      /   h


        if (shifted) then

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=4,n_max-3
                d1f(i)= A(3)*(f(i+2) - f(i-3))     &
                + A(2)*(f(i+1) - f(i-2))                    &
                + A(1)*(f(i) - f(i-1))
            enddo

            if (bc_type==periodic) then

                do i=1, 3
                    d1f(i)= A(3)*(f(i+2) - f(mod(n_max+i-5,n_max-1)+1))           &
                    + A(2)*(f(i+1) - f(mod(n_max+i-4,n_max-1)+1))           &
                    + A(1)*(f(i) - f(mod(n_max+i-3,n_max-1)+1))
                enddo

                do i=n_max-2,n_max-1
                    d1f(i)=   A(3)*(f(mod(i+1,n_max-1)+1) - f(i-3))             &
                    + A(2)*(f(mod(i,n_max-1)+1) - f(i-2))             &
                    + A(1)*(f(i) - f(i-1))
                enddo

            end if

            if (bc_type==symetric) then

                do i=1, 3
                    d1f(i)=   A(3)*(f(i+2) - f(max(i-3,abs(i-4))))           &
                    + A(2)*(f(i+1) - f(max(i-2,abs(i-3))))           &
                    + A(1)*(f(i) - f(max(i-1,abs(i-2))))
                enddo

                do i=n_max-2,n_max-1
                    d1f(i)=   A(3)*(f(min(i+2, 2*n_max-(i+3))) - f(i-3))             &
                    + A(2)*(f(min(i+1, 2*n_max-(i+2))) - f(i-2))             &
                    + A(1)*(f(i) - f(i-1))
                enddo

            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_INF(f, d1f, h, 3, shifted)
            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_SUP(f, d1f, n_max, h, 2, shifted)
            end if


        else

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=3,n_max-4
                d1f(i)= A(3)*(f(i+3) - f(i-2))            &
                + A(2)*(f(i+2) - f(i-1))                    &
                + A(1)*(f(i+1) - f(i))
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 2
                    d1f(i)= A(3)*(f(i+3) - f(mod(n_max+i-4,n_max-1)+1))   &
                    + A(2)*(f(i+2) - f(mod(n_max+i-3,n_max-1)+1))           &
                    + A(1)*(f(i+1) - f(i))
                enddo

                do i=n_max-3,n_max-1
                    d1f(i)=A(3)*(f(mod(i+2,n_max-1)+1) - f(i-2))     &
                    + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-1))             &
                    + A(1)*(f(mod(i,n_max-1)+1) - f(i))
                enddo

            end if

            if (bc_type.eq.antisymetric) then

                d1f(1)=   A(3)*(f(4) + f(3))           &
                + A(2)*(f(3) + f(2))           &
                + A(1)*(f(2) - f(1))

                d1f(1)=d1f(1)-2.d0*(A(2)+A(3))*f(1)

                d1f(2)=   A(3)*(f(5) + f(2))           &
                + A(2)*(f(4) - f(1))           &
                + A(1)*(f(3) - f(2))

                d1f(2)=d1f(2)-2.d0*(A(3))*f(1)


                d1f(n_max-1)=   A(3)*(-f(n_max-2) - f(n_max-3))             &
                + A(2)*(-f(n_max-1) - f(n_max-2))             &
                + A(1)*(f(n_max) - f(n_max-1))

                d1f(n_max-1)=d1f(n_max-1)+2.d0*(A(2)+A(3))*f(n_max)


                d1f(n_max-2)=   A(3)*(-f(n_max-1) - f(n_max-4))             &
                + A(2)*(f(n_max) - f(n_max-3))             &
                + A(1)*(f(n_max-1) - f(n_max-2))

                d1f(n_max-2)=d1f(n_max-2)+2.d0*(A(3))*f(n_max)


                d1f(n_max-3)=   A(3)*(f(n_max) - f(n_max-5))             &
                + A(2)*(f(n_max-1) - f(n_max-4))             &
                + A(1)*(f(n_max-2) - f(n_max-3))

            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_INF(f, d1f, h, 2, shifted)
            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_SUP(f, d1f, n_max, h, 3, shifted)
            end if

        end if

        return
    end subroutine

    subroutine D1s_Tamm_ACC(f, d1f, n_max, h, shifted, bc_type)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(inout)   :: d1f(:)

        integer :: i
        real(kind=8) A(3)

        A(1) =  1.189101951200031d0      /   h
        A(2) =  -0.073717642266682d0     /   h
        A(3) =  0.006410195120003d0      /   h


        if (shifted) then

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=4,n_max-3
                d1f(i)= d1f(i) + A(3)*(f(i+2) - f(i-3))     &
                + A(2)*(f(i+1) - f(i-2))                    &
                + A(1)*(f(i) - f(i-1))
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 3
                    d1f(i)= d1f(i) + A(3)*(f(i+2) - f(mod(n_max+i-5,n_max-1)+1))           &
                    + A(2)*(f(i+1) - f(mod(n_max+i-4,n_max-1)+1))           &
                    + A(1)*(f(i) - f(mod(n_max+i-3,n_max-1)+1))
                enddo

                do i=n_max-2,n_max-1
                    d1f(i)=   d1f(i) + A(3)*(f(mod(i+1,n_max-1)+1) - f(i-3))             &
                    + A(2)*(f(mod(i,n_max-1)+1) - f(i-2))             &
                    + A(1)*(f(i) - f(i-1))
                enddo

            end if

            if (bc_type==symetric) then

                do i=1, 3
                    d1f(i)= d1f(i) + A(3)*(f(i+2) - f(max(i-3,abs(i-4))))           &
                    + A(2)*(f(i+1) - f(max(i-2,abs(i-3))))           &
                    + A(1)*(f(i) - f(max(i-1,abs(i-2))))
                enddo

                do i=n_max-2,n_max-1
                    d1f(i)= d1f(i) + A(3)*(f(min(i+2, 2*n_max-(i+3))) - f(i-3))             &
                    + A(2)*(f(min(i+1, 2*n_max-(i+2))) - f(i-2))             &
                    + A(1)*(f(i) - f(i-1))
                enddo

            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_INF_ACC(f, d1f, h, 3, shifted)
            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_SUP_ACC(f, d1f, n_max, h, 2, shifted)
            end if


        else

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=3,n_max-4
                d1f(i)= d1f(i) + A(3)*(f(i+3) - f(i-2))            &
                + A(2)*(f(i+2) - f(i-1))                    &
                + A(1)*(f(i+1) - f(i))
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 2
                    d1f(i)= d1f(i) + A(3)*(f(i+3) - f(mod(n_max+i-4,n_max-1)+1))   &
                    + A(2)*(f(i+2) - f(mod(n_max+i-3,n_max-1)+1))           &
                    + A(1)*(f(i+1) - f(i))
                enddo

                do i=n_max-3,n_max-1
                    d1f(i)=d1f(i) + A(3)*(f(mod(i+2,n_max-1)+1) - f(i-2))     &
                    + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-1))             &
                    + A(1)*(f(mod(i,n_max-1)+1) - f(i))
                enddo

            end if

            if (bc_type.eq.antisymetric) then

                d1f(1)= d1f(1)+A(3)*(f(4) + f(3))           &
                + A(2)*(f(3) + f(2))           &
                + A(1)*(f(2) - f(1))

                d1f(2)= d1f(2)+A(3)*(f(5) + f(2))           &
                + A(2)*(f(4) - f(1))           &
                + A(1)*(f(3) - f(2))


                d1f(n_max-1)= d1f(n_max-1)+A(3)*(-f(n_max-2) - f(n_max-3))             &
                + A(2)*(-f(n_max-1) - f(n_max-2))             &
                + A(1)*(f(n_max) - f(n_max-1))


                d1f(n_max-2)= d1f(n_max-2) + A(3)*(-f(n_max-1) - f(n_max-4))             &
                + A(2)*(f(n_max) - f(n_max-3))             &
                + A(1)*(f(n_max-1) - f(n_max-2))


                d1f(n_max-3)= d1f(n_max-3)+A(3)*(f(n_max) - f(n_max-5))             &
                + A(2)*(f(n_max-1) - f(n_max-4))             &
                + A(1)*(f(n_max-2) - f(n_max-3))

            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_INF_ACC(f, d1f, h, 2, shifted)
            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_SUP_ACC(f, d1f, n_max, h, 3, shifted)
            end if

        end if

        return
    end subroutine

    subroutine D1s_Tamm_MULT(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, geom_coefs(:), f(:)
        real(kind=8), intent(out)   :: d1f(:)

        integer :: i
        real(kind=8) A(3)

        A(1) =  1.189101951200031d0      /   h
        A(2) =  -0.073717642266682d0     /   h
        A(3) =  0.006410195120003d0      /   h


        if (shifted) then

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=4,n_max-3
                d1f(i)= ( A(3)*(f(i+2) - f(i-3))     &
                + A(2)*(f(i+1) - f(i-2))                    &
                + A(1)*(f(i) - f(i-1)) ) * geom_coefs(i)
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 3
                    d1f(i)= ( A(3)*(f(i+2) - f(mod(n_max+i-5,n_max-1)+1))           &
                    + A(2)*(f(i+1) - f(mod(n_max+i-4,n_max-1)+1))           &
                    + A(1)*(f(i) - f(mod(n_max+i-3,n_max-1)+1)) ) * geom_coefs(i)
                enddo

                do i=n_max-2,n_max-1
                    d1f(i)=   ( A(3)*(f(mod(i+1,n_max-1)+1) - f(i-3))             &
                    + A(2)*(f(mod(i,n_max-1)+1) - f(i-2))             &
                    + A(1)*(f(i) - f(i-1)) ) * geom_coefs(i)
                enddo

            end if

            if (bc_type==symetric) then

                do i=1, 3
                    d1f(i)= ( A(3)*(f(i+2) - f(max(i-3,abs(i-4))))           &
                    + A(2)*(f(i+1) - f(max(i-2,abs(i-3))))           &
                    + A(1)*(f(i) - f(max(i-1,abs(i-2)))) ) * geom_coefs(i)
                enddo

                do i=n_max-2,n_max-1
                    d1f(i)=  ( A(3)*(f(min(i+2, 2*n_max-(i+3))) - f(i-3))             &
                    + A(2)*(f(min(i+1, 2*n_max-(i+2))) - f(i-2))             &
                    + A(1)*(f(i) - f(i-1)) ) * geom_coefs(i)
                enddo

            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_INF_MULT(f, d1f, h, 3, shifted, geom_coefs)
            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_SUP_MULT(f, d1f, n_max, h, 2, shifted, geom_coefs)
            end if


        else

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=3,n_max-4
                d1f(i)= ( A(3)*(f(i+3) - f(i-2))            &
                + A(2)*(f(i+2) - f(i-1))                    &
                + A(1)*(f(i+1) - f(i)) ) * geom_coefs(i)
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 2
                    d1f(i)= ( A(3)*(f(i+3) - f(mod(n_max+i-4,n_max-1)+1))   &
                    + A(2)*(f(i+2) - f(mod(n_max+i-3,n_max-1)+1))           &
                    + A(1)*(f(i+1) - f(i)) ) * geom_coefs(i)
                enddo

                do i=n_max-3,n_max-1
                    d1f(i)=( A(3)*(f(mod(i+2,n_max-1)+1) - f(i-2))     &
                    + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-1))             &
                    + A(1)*(f(mod(i,n_max-1)+1) - f(i)) ) * geom_coefs(i)
                enddo

            end if

            if (bc_type.eq.antisymetric) then

                d1f(1)=   (A(3)*(f(4) + f(3))           &
                + A(2)*(f(3) + f(2))           &
                + A(1)*(f(2) - f(1)) ) * geom_coefs(1)

                d1f(2)=   (A(3)*(f(5) + f(2))           &
                + A(2)*(f(4) - f(1))           &
                + A(1)*(f(3) - f(2)) ) * geom_coefs(2)


                d1f(n_max-1)=   (A(3)*(-f(n_max-2) - f(n_max-3))             &
                + A(2)*(-f(n_max-1) - f(n_max-2))             &
                + A(1)*(f(n_max) - f(n_max-1)) ) * geom_coefs(n_max-1)


                d1f(n_max-2)=   (A(3)*(-f(n_max-1) - f(n_max-4))             &
                + A(2)*(f(n_max) - f(n_max-3))             &
                + A(1)*(f(n_max-1) - f(n_max-2)) ) * geom_coefs(n_max-2)


                d1f(n_max-3)=   (A(3)*(f(n_max) - f(n_max-5))             &
                + A(2)*(f(n_max-1) - f(n_max-4))             &
                + A(1)*(f(n_max-2) - f(n_max-3)) ) * geom_coefs(n_max-3)

            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_INF_MULT(f, d1f, h, 2, shifted, geom_coefs)
            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_SUP_MULT(f, d1f, n_max, h, 3, shifted, geom_coefs)
            end if

        end if

        return

    end subroutine

    subroutine D1s_Tamm_MULT_ACC(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, geom_coefs(:), f(:)
        real(kind=8), intent(inout)   :: d1f(:)

        integer :: i
        real(kind=8) A(3)

        A(1) =  1.189101951200031d0      /   h
        A(2) =  -0.073717642266682d0     /   h
        A(3) =  0.006410195120003d0      /   h


        if (shifted) then

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=4,n_max-3
                d1f(i)= d1f(i) + ( A(3)*(f(i+2) - f(i-3))     &
                + A(2)*(f(i+1) - f(i-2))                    &
                + A(1)*(f(i) - f(i-1)) ) * geom_coefs(i)
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 3
                    d1f(i)= d1f(i) + ( A(3)*(f(i+2) - f(mod(n_max+i-5,n_max-1)+1))           &
                    + A(2)*(f(i+1) - f(mod(n_max+i-4,n_max-1)+1))           &
                    + A(1)*(f(i) - f(mod(n_max+i-3,n_max-1)+1)) ) * geom_coefs(i)
                enddo

                do i=n_max-2,n_max-1
                    d1f(i)=   d1f(i) + ( A(3)*(f(mod(i+1,n_max-1)+1) - f(i-3))             &
                    + A(2)*(f(mod(i,n_max-1)+1) - f(i-2))             &
                    + A(1)*(f(i) - f(i-1)) ) * geom_coefs(i)
                enddo

            end if

            if (bc_type==symetric) then

                do i=1, 3
                    d1f(i)= d1f(i) + ( A(3)*(f(i+2) - f(max(i-3,abs(i-4))))           &
                    + A(2)*(f(i+1) - f(max(i-2,abs(i-3))))           &
                    + A(1)*(f(i) - f(max(i-1,abs(i-2)))) ) * geom_coefs(i)
                enddo

                do i=n_max-2,n_max-1
                    d1f(i)=  d1f(i) + ( A(3)*(f(min(i+2, 2*n_max-(i+3))) - f(i-3))             &
                    + A(2)*(f(min(i+1, 2*n_max-(i+2))) - f(i-2))             &
                    + A(1)*(f(i) - f(i-1)) ) * geom_coefs(i)
                enddo

            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_INF_MULT_ACC(f, d1f, h, 3, shifted, geom_coefs)
            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_SUP_MULT_ACC(f, d1f, n_max, h, 2, shifted, geom_coefs)
            end if


        else

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=3,n_max-4
                d1f(i)= d1f(i) + ( A(3)*(f(i+3) - f(i-2))            &
                + A(2)*(f(i+2) - f(i-1))                    &
                + A(1)*(f(i+1) - f(i)) ) * geom_coefs(i)
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 2
                    d1f(i)= d1f(i) + ( A(3)*(f(i+3) - f(mod(n_max+i-4,n_max-1)+1))   &
                    + A(2)*(f(i+2) - f(mod(n_max+i-3,n_max-1)+1))           &
                    + A(1)*(f(i+1) - f(i)) ) * geom_coefs(i)
                enddo

                do i=n_max-3,n_max-1
                    d1f(i)=d1f(i) + ( A(3)*(f(mod(i+2,n_max-1)+1) - f(i-2))     &
                    + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-1))             &
                    + A(1)*(f(mod(i,n_max-1)+1) - f(i)) ) * geom_coefs(i)
                enddo
            end if

            if (bc_type.eq.antisymetric) then

                d1f(1)=  d1f(1) + (A(3)*(f(4) + f(3))           &
                + A(2)*(f(3) + f(2))           &
                + A(1)*(f(2) - f(1)) ) * geom_coefs(1)

                d1f(2)=  d1f(2) + (A(3)*(f(5) + f(2))           &
                + A(2)*(f(4) - f(1))           &
                + A(1)*(f(3) - f(2)) ) * geom_coefs(2)


                d1f(n_max-1)=  d1f(n_max-1) + (A(3)*(-f(n_max-2) - f(n_max-3))             &
                + A(2)*(-f(n_max-1) - f(n_max-2))             &
                + A(1)*(f(n_max) - f(n_max-1)) ) * geom_coefs(n_max-1)


                d1f(n_max-2)=  d1f(n_max-2) + (A(3)*(-f(n_max-1) - f(n_max-4))             &
                + A(2)*(f(n_max) - f(n_max-3))             &
                + A(1)*(f(n_max-1) - f(n_max-2)) ) * geom_coefs(n_max-2)


                d1f(n_max-3)=  d1f(n_max-3) + (A(3)*(f(n_max) - f(n_max-5))             &
                + A(2)*(f(n_max-1) - f(n_max-4))             &
                + A(1)*(f(n_max-2) - f(n_max-3)) ) * geom_coefs(n_max-3)

            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_INF_MULT_ACC(f, d1f, h, 2, shifted, geom_coefs)
            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_SUP_MULT_ACC(f, d1f, n_max, h, 3, shifted, geom_coefs)
            end if

        end if

        return
    end subroutine

    subroutine D1s_ExpCtr_O0Fp5(f, d1f, n_max, h, shifted, bc_type)
        ! This schemes has been designed in the file
        ! FILTER_EXP_6pts.wxm. The different cases taken into
        ! account here are for different values of cutting freq. wc
        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)   :: d1f(:)

        integer :: i
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        if (shifted) then

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=6,n_max-5
                d1f(i)=   A(5)*(f(i+4) - f(i-5))    &
                + A(4)*(f(i+3) - f(i-4))            &
                + A(3)*(f(i+2) - f(i-3))            &
                + A(2)*(f(i+1) - f(i-2))            &
                + A(1)*(f(i) - f(i-1))
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 5
                    d1f(i)=   A(5)*(f(i+4) - f(mod(n_max+i-7,n_max-1)+1))   &
                    + A(4)*(f(i+3) - f(mod(n_max+i-6,n_max-1)+1))           &
                    + A(3)*(f(i+2) - f(mod(n_max+i-5,n_max-1)+1))           &
                    + A(2)*(f(i+1) - f(mod(n_max+i-4,n_max-1)+1))           &
                    + A(1)*(f(i) - f(mod(n_max+i-3,n_max-1)+1))
                enddo

                do i=n_max-4,n_max-1
                    d1f(i)=   A(5)*(f(mod(i+3,n_max-1)+1) - f(i-5))     &
                    + A(4)*(f(mod(i+2,n_max-1)+1) - f(i-4))             &
                    + A(3)*(f(mod(i+1,n_max-1)+1) - f(i-3))             &
                    + A(2)*(f(mod(i,n_max-1)+1) - f(i-2))             &
                    + A(1)*(f(i) - f(i-1))
                enddo

            end if

            if (bc_type==symetric) then

                do i=1, 5
                    d1f(i)=   A(5)*(f(i+4) - f(max(i-5,abs(i-6))))   &
                    + A(4)*(f(i+3) - f(max(i-4,abs(i-5))))           &
                    + A(3)*(f(i+2) - f(max(i-3,abs(i-4))))           &
                    + A(2)*(f(i+1) - f(max(i-2,abs(i-3))))           &
                    + A(1)*(f(i) - f(max(i-1,abs(i-2))))
                enddo

                do i=n_max-4,n_max-1
                    d1f(i)=   A(5)*(f(min(i+4, 2*n_max-(i+5))) - f(i-5))     &
                    + A(4)*(f(min(i+3, 2*n_max-(i+4))) - f(i-4))             &
                    + A(3)*(f(min(i+2, 2*n_max-(i+3))) - f(i-3))             &
                    + A(2)*(f(min(i+1, 2*n_max-(i+2))) - f(i-2))             &
                    + A(1)*(f(i) - f(i-1))
                enddo

            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_INF(f, d1f, h, 5, shifted)
            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_SUP(f, d1f, n_max, h, 4, shifted)
            end if


        else

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=5,n_max-6
                d1f(i)=   A(5)*(f(i+5) - f(i-4))    &
                + A(4)*(f(i+4) - f(i-3))            &
                + A(3)*(f(i+3) - f(i-2))            &
                + A(2)*(f(i+2) - f(i-1))            &
                + A(1)*(f(i+1) - f(i))
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 4
                    d1f(i)=   A(5)*(f(i+5) - f(mod(n_max+i-6,n_max-1)+1))   &
                    + A(4)*(f(i+4) - f(mod(n_max+i-5,n_max-1)+1))           &
                    + A(3)*(f(i+3) - f(mod(n_max+i-4,n_max-1)+1))           &
                    + A(2)*(f(i+2) - f(mod(n_max+i-3,n_max-1)+1))           &
                    + A(1)*(f(i+1) - f(i))
                enddo

                do i=n_max-5,n_max-1
                    d1f(i)=   A(5)*(f(mod(i+4,n_max-1)+1) - f(i-4))     &
                    + A(4)*(f(mod(i+3,n_max-1)+1) - f(i-3))             &
                    + A(3)*(f(mod(i+2,n_max-1)+1) - f(i-2))             &
                    + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-1))             &
                    + A(1)*(f(mod(i,n_max-1)+1) - f(i))
                enddo
            end if

            if (bc_type.eq.antisymetric) then

                d1f(1)=   A(5)*(f(6) + f(5))   &
                + A(4)*(f(5) + f(4))           &
                + A(3)*(f(4) + f(3))           &
                + A(2)*(f(3) + f(2))           &
                + A(1)*(f(2) - f(1))

                d1f(1)=d1f(1)-2.d0*(A(2)+A(3)+A(4)+A(5))*f(1)

                d1f(2)=   A(5)*(f(7) + f(4))   &
                + A(4)*(f(6) + f(3))           &
                + A(3)*(f(5) + f(2))           &
                + A(2)*(f(4) - f(1))           &
                + A(1)*(f(3) - f(2))

                d1f(2)=d1f(2)-2.d0*(A(3)+A(4)+A(5))*f(1)

                d1f(3)=   A(5)*(f(8) + f(3))   &
                + A(4)*(f(7) + f(2))           &
                + A(3)*(f(6) - f(1))           &
                + A(2)*(f(5) - f(2))           &
                + A(1)*(f(4) - f(3))

                d1f(3)=d1f(3)-2.d0*(A(4)+A(5))*f(1)

                d1f(4)=   A(5)*(f(9) + f(2))   &
                + A(4)*(f(8) - f(1))           &
                + A(3)*(f(7) - f(2))           &
                + A(2)*(f(6) - f(3))           &
                + A(1)*(f(5) - f(4))

                d1f(4)=d1f(4)-2.d0*(A(5))*f(1)


                d1f(n_max-1)=   A(5)*(-f(n_max-4) - f(n_max-5))     &
                + A(4)*(-f(n_max-3) - f(n_max-4))             &
                + A(3)*(-f(n_max-2) - f(n_max-3))             &
                + A(2)*(-f(n_max-1) - f(n_max-2))             &
                + A(1)*(f(n_max) - f(n_max-1))

                d1f(n_max-1)=d1f(n_max-1)+2.d0*(A(2)+A(3)+A(4)+A(5))*f(n_max)


                d1f(n_max-2)=   A(5)*(-f(n_max-3) - f(n_max-6))     &
                + A(4)*(-f(n_max-2) - f(n_max-5))             &
                + A(3)*(-f(n_max-1) - f(n_max-4))             &
                + A(2)*(f(n_max) - f(n_max-3))             &
                + A(1)*(f(n_max-1) - f(n_max-2))

                d1f(n_max-2)=d1f(n_max-2)+2.d0*(A(3)+A(4)+A(5))*f(n_max)


                d1f(n_max-3)=   A(5)*(-f(n_max-2) - f(n_max-7))     &
                + A(4)*(-f(n_max-1) - f(n_max-6))             &
                + A(3)*(f(n_max) - f(n_max-5))             &
                + A(2)*(f(n_max-1) - f(n_max-4))             &
                + A(1)*(f(n_max-2) - f(n_max-3))

                d1f(n_max-3)=d1f(n_max-3)+2.d0*(A(4)+A(5))*f(n_max)


                d1f(n_max-4)=   A(5)*(-f(n_max-1) - f(n_max-8))     &
                + A(4)*(f(n_max) - f(n_max-7))             &
                + A(3)*(f(n_max-1) - f(n_max-6))             &
                + A(2)*(f(n_max-2) - f(n_max-5))             &
                + A(1)*(f(n_max-3) - f(n_max-4))

                d1f(n_max-4)=d1f(n_max-4)+2.d0*(A(5))*f(n_max)


                d1f(n_max-5)=   A(5)*(f(n_max) - f(n_max-9))     &
                + A(4)*(f(n_max-1) - f(n_max-8))             &
                + A(3)*(f(n_max-2) - f(n_max-7))             &
                + A(2)*(f(n_max-3) - f(n_max-6))             &
                + A(1)*(f(n_max-4) - f(n_max-5))

            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_INF(f, d1f, h, 4, shifted)
            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_SUP(f, d1f, n_max, h, 5, shifted)
            end if

        end if

        return

    end subroutine

    subroutine D1s_ExpCtr_O0Fp5_ACC(f, d1f, n_max, h, shifted, bc_type)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(inout)   :: d1f(:)

        integer :: i
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        if (shifted) then

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=6,n_max-5
                d1f(i)=d1f(i) + A(5)*(f(i+4) - f(i-5))    &
                + A(4)*(f(i+3) - f(i-4))            &
                + A(3)*(f(i+2) - f(i-3))            &
                + A(2)*(f(i+1) - f(i-2))            &
                + A(1)*(f(i) - f(i-1))
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 5
                    d1f(i)= d1f(i) + A(5)*(f(i+4) - f(mod(n_max+i-7,n_max-1)+1))   &
                    + A(4)*(f(i+3) - f(mod(n_max+i-6,n_max-1)+1))           &
                    + A(3)*(f(i+2) - f(mod(n_max+i-5,n_max-1)+1))           &
                    + A(2)*(f(i+1) - f(mod(n_max+i-4,n_max-1)+1))           &
                    + A(1)*(f(i) - f(mod(n_max+i-3,n_max-1)+1))
                enddo

                do i=n_max-4,n_max-1
                    d1f(i)= d1f(i) + A(5)*(f(mod(i+3,n_max-1)+1) - f(i-5))     &
                    + A(4)*(f(mod(i+2,n_max-1)+1) - f(i-4))             &
                    + A(3)*(f(mod(i+1,n_max-1)+1) - f(i-3))             &
                    + A(2)*(f(mod(i,n_max-1)+1) - f(i-2))             &
                    + A(1)*(f(i) - f(i-1))
                enddo
            end if

            if (bc_type==symetric) then

                do i=1, 5
                    d1f(i)= d1f(i) + A(5)*(f(i+4) - f(max(i-5,abs(i-6))))   &
                    + A(4)*(f(i+3) - f(max(i-4,abs(i-5))))           &
                    + A(3)*(f(i+2) - f(max(i-3,abs(i-4))))           &
                    + A(2)*(f(i+1) - f(max(i-2,abs(i-3))))           &
                    + A(1)*(f(i) - f(max(i-1,abs(i-2))))
                enddo

                do i=n_max-4,n_max-1
                    d1f(i)= d1f(i) + A(5)*(f(min(i+4, 2*n_max-(i+5))) - f(i-5))     &
                    + A(4)*(f(min(i+3, 2*n_max-(i+4))) - f(i-4))             &
                    + A(3)*(f(min(i+2, 2*n_max-(i+3))) - f(i-3))             &
                    + A(2)*(f(min(i+1, 2*n_max-(i+2))) - f(i-2))             &
                    + A(1)*(f(i) - f(i-1))
                enddo

            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_INF_ACC(f, d1f, h, 5, shifted)
            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_SUP_ACC(f, d1f, n_max, h, 4, shifted)
            end if


        else

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=5,n_max-6
                d1f(i)= d1f(i) + A(5)*(f(i+5) - f(i-4))    &
                + A(4)*(f(i+4) - f(i-3))            &
                + A(3)*(f(i+3) - f(i-2))            &
                + A(2)*(f(i+2) - f(i-1))            &
                + A(1)*(f(i+1) - f(i))
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 4
                    d1f(i)= d1f(i) + A(5)*(f(i+5) - f(mod(n_max+i-6,n_max-1)+1))   &
                    + A(4)*(f(i+4) - f(mod(n_max+i-5,n_max-1)+1))           &
                    + A(3)*(f(i+3) - f(mod(n_max+i-4,n_max-1)+1))           &
                    + A(2)*(f(i+2) - f(mod(n_max+i-3,n_max-1)+1))           &
                    + A(1)*(f(i+1) - f(i))
                enddo

                do i=n_max-5,n_max-1
                    d1f(i)= d1f(i) + A(5)*(f(mod(i+4,n_max-1)+1) - f(i-4))     &
                    + A(4)*(f(mod(i+3,n_max-1)+1) - f(i-3))             &
                    + A(3)*(f(mod(i+2,n_max-1)+1) - f(i-2))             &
                    + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-1))             &
                    + A(1)*(f(mod(i,n_max-1)+1) - f(i))
                enddo

            end if

            if (bc_type.eq.antisymetric) then

                d1f(1)= d1f(1)+A(5)*(f(6) + f(5))   &
                + A(4)*(f(5) + f(4))           &
                + A(3)*(f(4) + f(3))           &
                + A(2)*(f(3) + f(2))           &
                + A(1)*(f(2) - f(1))

                d1f(1)=d1f(1)-2.d0*(A(2)+A(3)+A(4)+A(5))*f(1)

                d1f(2)= d1f(2)+A(5)*(f(7) + f(4))   &
                + A(4)*(f(6) + f(3))           &
                + A(3)*(f(5) + f(2))           &
                + A(2)*(f(4) - f(1))           &
                + A(1)*(f(3) - f(2))

                d1f(2)=d1f(2)-2.d0*(A(3)+A(4)+A(5))*f(1)

                d1f(3)= d1f(3)+A(5)*(f(8) + f(3))   &
                + A(4)*(f(7) + f(2))           &
                + A(3)*(f(6) - f(1))           &
                + A(2)*(f(5) - f(2))           &
                + A(1)*(f(4) - f(3))

                d1f(3)=d1f(3)-2.d0*(A(4)+A(5))*f(1)

                d1f(4)= d1f(4)+A(5)*(f(9) + f(2))   &
                + A(4)*(f(8) - f(1))           &
                + A(3)*(f(7) - f(2))           &
                + A(2)*(f(6) - f(3))           &
                + A(1)*(f(5) - f(4))

                d1f(4)=d1f(4)-2.d0*(A(5))*f(1)


                d1f(n_max-1)= d1f(n_max-1)+A(5)*(-f(n_max-4) - f(n_max-5))     &
                + A(4)*(-f(n_max-3) - f(n_max-4))             &
                + A(3)*(-f(n_max-2) - f(n_max-3))             &
                + A(2)*(-f(n_max-1) - f(n_max-2))             &
                + A(1)*(f(n_max) - f(n_max-1))

                d1f(n_max-1)=d1f(n_max-1)+2.d0*(A(2)+A(3)+A(4)+A(5))*f(n_max)


                d1f(n_max-2)= d1f(n_max-2)+A(5)*(-f(n_max-3) - f(n_max-6))     &
                + A(4)*(-f(n_max-2) - f(n_max-5))             &
                + A(3)*(-f(n_max-1) - f(n_max-4))             &
                + A(2)*(f(n_max) - f(n_max-3))             &
                + A(1)*(f(n_max-1) - f(n_max-2))

                d1f(n_max-2)=d1f(n_max-2)+2.d0*(A(3)+A(4)+A(5))*f(n_max)


                d1f(n_max-3)= d1f(n_max-3)+A(5)*(-f(n_max-2) - f(n_max-7))     &
                + A(4)*(-f(n_max-1) - f(n_max-6))             &
                + A(3)*(f(n_max) - f(n_max-5))             &
                + A(2)*(f(n_max-1) - f(n_max-4))             &
                + A(1)*(f(n_max-2) - f(n_max-3))

                d1f(n_max-3)=d1f(n_max-3)+2.d0*(A(4)+A(5))*f(n_max)


                d1f(n_max-4)= d1f(n_max-4)+A(5)*(-f(n_max-1) - f(n_max-8))     &
                + A(4)*(f(n_max) - f(n_max-7))             &
                + A(3)*(f(n_max-1) - f(n_max-6))             &
                + A(2)*(f(n_max-2) - f(n_max-5))             &
                + A(1)*(f(n_max-3) - f(n_max-4))

                d1f(n_max-4)=d1f(n_max-4)+2.d0*(A(5))*f(n_max)


                d1f(n_max-5)= d1f(n_max-5)+A(5)*(f(n_max) - f(n_max-9))     &
                + A(4)*(f(n_max-1) - f(n_max-8))             &
                + A(3)*(f(n_max-2) - f(n_max-7))             &
                + A(2)*(f(n_max-3) - f(n_max-6))             &
                + A(1)*(f(n_max-4) - f(n_max-5))

            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_INF_ACC(f, d1f, h, 4, shifted)
            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_SUP_ACC(f, d1f, n_max, h, 5, shifted)
            end if

        end if

        return

    end subroutine

    subroutine D1s_ExpCtr_O0Fp5_MULT(f, d1f, n_max, h, shifted, bc_type, geom_coef)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, geom_coef(:), f(:)
        real(kind=8), intent(out)   :: d1f(:)

        integer :: i
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        if (shifted) then

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=6,n_max-5
                d1f(i)=(A(5)*(f(i+4) - f(i-5))    &
                + A(4)*(f(i+3) - f(i-4))            &
                + A(3)*(f(i+2) - f(i-3))            &
                + A(2)*(f(i+1) - f(i-2))            &
                + A(1)*(f(i) - f(i-1))) * geom_coef(i)
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 5
                    d1f(i)= ( A(5)*(f(i+4) - f(mod(n_max+i-7,n_max-1)+1))   &
                    + A(4)*(f(i+3) - f(mod(n_max+i-6,n_max-1)+1))           &
                    + A(3)*(f(i+2) - f(mod(n_max+i-5,n_max-1)+1))           &
                    + A(2)*(f(i+1) - f(mod(n_max+i-4,n_max-1)+1))           &
                    + A(1)*(f(i) - f(mod(n_max+i-3,n_max-1)+1)) ) * geom_coef(i)
                enddo

                do i=n_max-4,n_max-1
                    d1f(i)= ( A(5)*(f(mod(i+3,n_max-1)+1) - f(i-5))     &
                    + A(4)*(f(mod(i+2,n_max-1)+1) - f(i-4))             &
                    + A(3)*(f(mod(i+1,n_max-1)+1) - f(i-3))             &
                    + A(2)*(f(mod(i,n_max-1)+1) - f(i-2))               &
                    + A(1)*(f(i) - f(i-1)) ) * geom_coef(i)
                enddo

            end if

            if (bc_type==symetric) then

                do i=1, 5
                    d1f(i)= ( A(5)*(f(i+4) - f(max(i-5,abs(i-6))))   &
                    + A(4)*(f(i+3) - f(max(i-4,abs(i-5))))           &
                    + A(3)*(f(i+2) - f(max(i-3,abs(i-4))))           &
                    + A(2)*(f(i+1) - f(max(i-2,abs(i-3))))           &
                    + A(1)*(f(i) - f(max(i-1,abs(i-2)))) ) * geom_coef(i)
                enddo

                do i=n_max-4,n_max-1
                    d1f(i)=  ( A(5)*(f(min(i+4, 2*n_max-(i+5))) - f(i-5))     &
                    + A(4)*(f(min(i+3, 2*n_max-(i+4))) - f(i-4))             &
                    + A(3)*(f(min(i+2, 2*n_max-(i+3))) - f(i-3))             &
                    + A(2)*(f(min(i+1, 2*n_max-(i+2))) - f(i-2))             &
                    + A(1)*(f(i) - f(i-1)) ) * geom_coef(i)
                enddo

            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_INF_MULT(f, d1f, h, 5, shifted, geom_coef)
            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_SUP_MULT(f, d1f, n_max, h, 4, shifted, geom_coef)
            end if


        else

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=5,n_max-6
                d1f(i)= ( A(5)*(f(i+5) - f(i-4))    &
                + A(4)*(f(i+4) - f(i-3))            &
                + A(3)*(f(i+3) - f(i-2))            &
                + A(2)*(f(i+2) - f(i-1))            &
                + A(1)*(f(i+1) - f(i)) ) * geom_coef(i)
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 4
                    d1f(i)= (A(5)*(f(i+5) - f(mod(n_max+i-6,n_max-1)+1))   &
                    + A(4)*(f(i+4) - f(mod(n_max+i-5,n_max-1)+1))           &
                    + A(3)*(f(i+3) - f(mod(n_max+i-4,n_max-1)+1))           &
                    + A(2)*(f(i+2) - f(mod(n_max+i-3,n_max-1)+1))           &
                    + A(1)*(f(i+1) - f(i)) ) * geom_coef(i)
                enddo

                do i=n_max-5,n_max-1
                    d1f(i)= ( A(5)*(f(mod(i+4,n_max-1)+1) - f(i-4))     &
                    + A(4)*(f(mod(i+3,n_max-1)+1) - f(i-3))             &
                    + A(3)*(f(mod(i+2,n_max-1)+1) - f(i-2))             &
                    + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-1))             &
                    + A(1)*(f(mod(i,n_max-1)+1) - f(i)) ) * geom_coef(i)
                enddo

            end if

            if (bc_type.eq.antisymetric) then

                d1f(1)=   (A(5)*(f(6) + f(5))   &
                + A(4)*(f(5) + f(4))           &
                + A(3)*(f(4) + f(3))           &
                + A(2)*(f(3) + f(2))           &
                + A(1)*(f(2) - f(1)) ) * geom_coef(1)

                d1f(1)=d1f(1)-2.d0*(A(2)+A(3)+A(4)+A(5))*f(1)* geom_coef(1)

                d1f(2)=   (A(5)*(f(7) + f(4))   &
                + A(4)*(f(6) + f(3))           &
                + A(3)*(f(5) + f(2))           &
                + A(2)*(f(4) - f(1))           &
                + A(1)*(f(3) - f(2)) ) * geom_coef(2)

                d1f(2)=d1f(2)-2.d0*(A(3)+A(4)+A(5))*f(1)* geom_coef(2)

                d1f(3)=   (A(5)*(f(8) + f(3))   &
                + A(4)*(f(7) + f(2))           &
                + A(3)*(f(6) - f(1))           &
                + A(2)*(f(5) - f(2))           &
                + A(1)*(f(4) - f(3)) ) * geom_coef(3)

                d1f(3)=d1f(3)-2.d0*(A(4)+A(5))*f(1)* geom_coef(3)

                d1f(4)=   (A(5)*(f(9) + f(2))   &
                + A(4)*(f(8) - f(1))           &
                + A(3)*(f(7) - f(2))           &
                + A(2)*(f(6) - f(3))           &
                + A(1)*(f(5) - f(4)) ) * geom_coef(4)

                d1f(4)=d1f(4)-2.d0*(A(5))*f(1)* geom_coef(4)


                d1f(n_max-1)=   (A(5)*(-f(n_max-4) - f(n_max-5))     &
                + A(4)*(-f(n_max-3) - f(n_max-4))             &
                + A(3)*(-f(n_max-2) - f(n_max-3))             &
                + A(2)*(-f(n_max-1) - f(n_max-2))             &
                + A(1)*(f(n_max) - f(n_max-1)) ) * geom_coef(n_max-1)

                d1f(n_max-1)=d1f(n_max-1)+2.d0*(A(2)+A(3)+A(4)+A(5))*f(n_max)* geom_coef(n_max-1)


                d1f(n_max-2)=   (A(5)*(-f(n_max-3) - f(n_max-6))     &
                + A(4)*(-f(n_max-2) - f(n_max-5))             &
                + A(3)*(-f(n_max-1) - f(n_max-4))             &
                + A(2)*(f(n_max) - f(n_max-3))             &
                + A(1)*(f(n_max-1) - f(n_max-2)) ) * geom_coef(n_max-2)

                d1f(n_max-2)=d1f(n_max-2)+2.d0*(A(3)+A(4)+A(5))*f(n_max)* geom_coef(n_max-2)


                d1f(n_max-3)=   (A(5)*(-f(n_max-2) - f(n_max-7))     &
                + A(4)*(-f(n_max-1) - f(n_max-6))             &
                + A(3)*(f(n_max) - f(n_max-5))             &
                + A(2)*(f(n_max-1) - f(n_max-4))             &
                + A(1)*(f(n_max-2) - f(n_max-3)) ) * geom_coef(n_max-3)

                d1f(n_max-3)=d1f(n_max-3)+2.d0*(A(4)+A(5))*f(n_max)* geom_coef(n_max-3)


                d1f(n_max-4)=   (A(5)*(-f(n_max-1) - f(n_max-8))     &
                + A(4)*(f(n_max) - f(n_max-7))             &
                + A(3)*(f(n_max-1) - f(n_max-6))             &
                + A(2)*(f(n_max-2) - f(n_max-5))             &
                + A(1)*(f(n_max-3) - f(n_max-4)) ) * geom_coef(n_max-4)

                d1f(n_max-4)=d1f(n_max-4)+2.d0*(A(5))*f(n_max)* geom_coef(n_max-4)


                d1f(n_max-5)=   (A(5)*(f(n_max) - f(n_max-9))     &
                + A(4)*(f(n_max-1) - f(n_max-8))             &
                + A(3)*(f(n_max-2) - f(n_max-7))             &
                + A(2)*(f(n_max-3) - f(n_max-6))             &
                + A(1)*(f(n_max-4) - f(n_max-5)) ) * geom_coef(n_max-5)

            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_INF_MULT(f, d1f, h, 4, shifted, geom_coef)
            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_SUP_MULT(f, d1f, n_max, h, 5, shifted, geom_coef)
            end if

        end if

        return

    end subroutine

    subroutine D1s_ExpCtr_O0Fp5_MULT_ACC(f, d1f, n_max, h, shifted, bc_type, geom_coef)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, geom_coef(:), f(:)
        real(kind=8), intent(inout)   :: d1f(:)

        integer :: i
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        if (shifted) then

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=6,n_max-5
                d1f(i)=d1f(i) + (A(5)*(f(i+4) - f(i-5))    &
                + A(4)*(f(i+3) - f(i-4))            &
                + A(3)*(f(i+2) - f(i-3))            &
                + A(2)*(f(i+1) - f(i-2))            &
                + A(1)*(f(i) - f(i-1))) * geom_coef(i)
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 5
                    d1f(i)= d1f(i) + ( A(5)*(f(i+4) - f(mod(n_max+i-7,n_max-1)+1))   &
                    + A(4)*(f(i+3) - f(mod(n_max+i-6,n_max-1)+1))           &
                    + A(3)*(f(i+2) - f(mod(n_max+i-5,n_max-1)+1))           &
                    + A(2)*(f(i+1) - f(mod(n_max+i-4,n_max-1)+1))           &
                    + A(1)*(f(i) - f(mod(n_max+i-3,n_max-1)+1)) ) * geom_coef(i)
                enddo

                do i=n_max-4,n_max-1
                    d1f(i)= d1f(i) + ( A(5)*(f(mod(i+3,n_max-1)+1) - f(i-5))     &
                    + A(4)*(f(mod(i+2,n_max-1)+1) - f(i-4))             &
                    + A(3)*(f(mod(i+1,n_max-1)+1) - f(i-3))             &
                    + A(2)*(f(mod(i,n_max-1)+1) - f(i-2))               &
                    + A(1)*(f(i) - f(i-1)) ) * geom_coef(i)
                enddo

            end if

            if (bc_type==symetric) then

                do i=1, 5
                    d1f(i)= d1f(i) + ( A(5)*(f(i+4) - f(max(i-5,abs(i-6))))   &
                    + A(4)*(f(i+3) - f(max(i-4,abs(i-5))))           &
                    + A(3)*(f(i+2) - f(max(i-3,abs(i-4))))           &
                    + A(2)*(f(i+1) - f(max(i-2,abs(i-3))))           &
                    + A(1)*(f(i) - f(max(i-1,abs(i-2)))) ) * geom_coef(i)
                enddo

                do i=n_max-4,n_max-1
                    d1f(i)=  d1f(i) + ( A(5)*(f(min(i+4, 2*n_max-(i+5))) - f(i-5))     &
                    + A(4)*(f(min(i+3, 2*n_max-(i+4))) - f(i-4))             &
                    + A(3)*(f(min(i+2, 2*n_max-(i+3))) - f(i-3))             &
                    + A(2)*(f(min(i+1, 2*n_max-(i+2))) - f(i-2))             &
                    + A(1)*(f(i) - f(i-1)) ) * geom_coef(i)
                enddo

            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_INF_MULT_ACC(f, d1f, h, 5, shifted, geom_coef)
            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_SUP_MULT_ACC(f, d1f, n_max, h, 4, shifted, geom_coef)
            end if


        else

            ! do i=N, n_max-N-1 (N size of stencil)
            do i=5,n_max-6
                d1f(i)= d1f(i) + ( A(5)*(f(i+5) - f(i-4))    &
                + A(4)*(f(i+4) - f(i-3))            &
                + A(3)*(f(i+3) - f(i-2))            &
                + A(2)*(f(i+2) - f(i-1))            &
                + A(1)*(f(i+1) - f(i)) ) * geom_coef(i)
            enddo

            if (bc_type.eq.periodic) then

                do i=1, 4
                    d1f(i)= d1f(i) + (A(5)*(f(i+5) - f(mod(n_max+i-6,n_max-1)+1))   &
                    + A(4)*(f(i+4) - f(mod(n_max+i-5,n_max-1)+1))           &
                    + A(3)*(f(i+3) - f(mod(n_max+i-4,n_max-1)+1))           &
                    + A(2)*(f(i+2) - f(mod(n_max+i-3,n_max-1)+1))           &
                    + A(1)*(f(i+1) - f(i)) ) * geom_coef(i)
                enddo

                do i=n_max-5,n_max-1
                    d1f(i)= d1f(i) + ( A(5)*(f(mod(i+4,n_max-1)+1) - f(i-4))     &
                    + A(4)*(f(mod(i+3,n_max-1)+1) - f(i-3))             &
                    + A(3)*(f(mod(i+2,n_max-1)+1) - f(i-2))             &
                    + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-1))             &
                    + A(1)*(f(mod(i,n_max-1)+1) - f(i)) ) * geom_coef(i)
                enddo

            end if

            if (bc_type.eq.antisymetric) then

                d1f(1)=  d1f(1) + (A(5)*(f(6) + f(5))   &
                + A(4)*(f(5) + f(4))           &
                + A(3)*(f(4) + f(3))           &
                + A(2)*(f(3) + f(2))           &
                + A(1)*(f(2) - f(1)) ) * geom_coef(1)

                d1f(1)=d1f(1)-2.d0*(A(2)+A(3)+A(4)+A(5))*f(1)* geom_coef(1)

                d1f(2)=  d1f(2) + (A(5)*(f(7) + f(4))   &
                + A(4)*(f(6) + f(3))           &
                + A(3)*(f(5) + f(2))           &
                + A(2)*(f(4) - f(1))           &
                + A(1)*(f(3) - f(2)) ) * geom_coef(2)

                d1f(2)=d1f(2)-2.d0*(A(3)+A(4)+A(5))*f(1)* geom_coef(2)

                d1f(3)=  d1f(3) + (A(5)*(f(8) + f(3))   &
                + A(4)*(f(7) + f(2))           &
                + A(3)*(f(6) - f(1))           &
                + A(2)*(f(5) - f(2))           &
                + A(1)*(f(4) - f(3)) ) * geom_coef(3)

                d1f(3)=d1f(3)-2.d0*(A(4)+A(5))*f(1)* geom_coef(3)

                d1f(4)=  d1f(4) + (A(5)*(f(9) + f(2))   &
                + A(4)*(f(8) - f(1))           &
                + A(3)*(f(7) - f(2))           &
                + A(2)*(f(6) - f(3))           &
                + A(1)*(f(5) - f(4)) ) * geom_coef(4)

                d1f(4)=d1f(4)-2.d0*(A(5))*f(1)* geom_coef(4)


                d1f(n_max-1)=  d1f(n_max-1) + (A(5)*(-f(n_max-4) - f(n_max-5))     &
                + A(4)*(-f(n_max-3) - f(n_max-4))             &
                + A(3)*(-f(n_max-2) - f(n_max-3))             &
                + A(2)*(-f(n_max-1) - f(n_max-2))             &
                + A(1)*(f(n_max) - f(n_max-1)) ) * geom_coef(n_max-1)

                d1f(n_max-1)=d1f(n_max-1)+2.d0*(A(2)+A(3)+A(4)+A(5))*f(n_max)* geom_coef(n_max-1)


                d1f(n_max-2)=  d1f(n_max-2) + (A(5)*(-f(n_max-3) - f(n_max-6))     &
                + A(4)*(-f(n_max-2) - f(n_max-5))             &
                + A(3)*(-f(n_max-1) - f(n_max-4))             &
                + A(2)*(f(n_max) - f(n_max-3))             &
                + A(1)*(f(n_max-1) - f(n_max-2)) ) * geom_coef(n_max-2)

                d1f(n_max-2)=d1f(n_max-2)+2.d0*(A(3)+A(4)+A(5))*f(n_max)* geom_coef(n_max-2)


                d1f(n_max-3)=  d1f(n_max-3) + (A(5)*(-f(n_max-2) - f(n_max-7))     &
                + A(4)*(-f(n_max-1) - f(n_max-6))             &
                + A(3)*(f(n_max) - f(n_max-5))             &
                + A(2)*(f(n_max-1) - f(n_max-4))             &
                + A(1)*(f(n_max-2) - f(n_max-3)) ) * geom_coef(n_max-3)

                d1f(n_max-3)=d1f(n_max-3)+2.d0*(A(4)+A(5))*f(n_max)* geom_coef(n_max-3)


                d1f(n_max-4)=  d1f(n_max-4) + (A(5)*(-f(n_max-1) - f(n_max-8))     &
                + A(4)*(f(n_max) - f(n_max-7))             &
                + A(3)*(f(n_max-1) - f(n_max-6))             &
                + A(2)*(f(n_max-2) - f(n_max-5))             &
                + A(1)*(f(n_max-3) - f(n_max-4)) ) * geom_coef(n_max-4)

                d1f(n_max-4)=d1f(n_max-4)+2.d0*(A(5))*f(n_max)* geom_coef(n_max-4)


                d1f(n_max-5)=  d1f(n_max-5) + (A(5)*(f(n_max) - f(n_max-9))     &
                + A(4)*(f(n_max-1) - f(n_max-8))             &
                + A(3)*(f(n_max-2) - f(n_max-7))             &
                + A(2)*(f(n_max-3) - f(n_max-6))             &
                + A(1)*(f(n_max-4) - f(n_max-5)) ) * geom_coef(n_max-5)

            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_INF_MULT_ACC(f, d1f, h, 4, shifted, geom_coef)
            end if

            if (bc_type.eq.Dirichlet) then
                call D1s_CTR_CLOSURE_SUP_MULT_ACC(f, d1f, n_max, h, 5, shifted, geom_coef)
            end if

        end if

        return

    end subroutine

    subroutine D1s_CptCtr_O6Fp0(f, d1f, n_max, h, shifted, bc_type)

        ! Values localisation:
        ! |------x------o------x------o------x------o------x------|
        ! fs    df1     f1    df2     f2          fnmax  fnmax+1  fn
        use Thomas_solver

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)   :: d1f(:)

        integer :: i
        real(kind=8)  :: al
        real(kind=8)  :: rhs(n_max)

        ! A, B MATRIX
        real(kind=8), save:: B(-2:2), B_lim(-1:3)
        real(kind=8), dimension(:,:), save, allocatable   :: A

        integer, save :: last_nmax, last_bc_type
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h



        update_matrix=(last_nmax/=n_max).or.(last_bc_type/=bc_type)    &
        .or.(last_shifted.neqv.shifted).or.(last_h/=h)

        !if (update_matrix) then
        call fill_matrix
        !end if

        last_nmax=n_max
        last_bc_type=bc_type
        last_shifted=shifted
        last_h=h

        if (shifted) then

            do i=3,n_max-2
                rhs(i)= ( B(2)*f(i+1) + B(1)*f(i) + B(-1)*f(i-1) + B(-2)*f(i-2) )
            enddo

            select case (bc_type)

                case (periodic)

                    rhs(1)= ( B(2)*f(2) + B(1)*f(1) + B(-1)*f(n_max-1) + B(-2)*f(n_max-2) )
                    rhs(2)= ( B(2)*f(3) + B(1)*f(2) + B(-1)*f(1) + B(-2)*f(n_max-1) )

                    rhs(n_max-1)= ( B(2)*f(1) + B(1)*f(n_max-1) + B(-1)*f(n_max-2) + B(-2)*f(n_max-3) )

                    call TS_Pr(A(-1,1:n_max-1), A(0,1:n_max-1), A(1,1:n_max-1), rhs(1:n_max-1), d1f(1:n_max-1), n_max-1)

                case (Dirichlet)

                    rhs(2)= ( B_lim(-1)*f(1) + B_lim(0) * f(2) + B_lim(1) * f(3) )

                    rhs(n_max-1)= -( + B_lim(-1)*f(n_max-1) + B_lim(0) * f(n_max-2) + B_lim(1)*f(n_max-3) )

                    call TS_NPr(rhs(2:n_max-1), d1f(2:n_max-1), n_max-2)

                case (symetric)
                    rhs(2)= ( B(2)*f(3) + B(1)*f(2) + B(-1)*f(1) + B(-2)*f(1) )

                    rhs(n_max-1)= ( B(2)*f(n_max-1) + B(1)*f(n_max-1) + B(-1)*f(n_max-2) + B(-2)*f(n_max-3) )


                    call TS_NPr(rhs(2:n_max-1), d1f(2:n_max-1), n_max-2)

            end select

        else

            do i=2,n_max-3
                rhs(i)= ( B(2)*f(i+2) + B(1)*f(i+1) + B(-1)*f(i) + B(-2)*f(i-1) )
            enddo

            select case (bc_type)

                case (periodic)

                    rhs(1)= ( B(2)*f(3) + B(1)*f(2) + B(-1)*f(1) + B(-2)*f(n_max-1) )

                    rhs(n_max-2)= ( B(2)*f(1) + B(1)*f(n_max-1) + B(-1)*f(n_max-2) + B(-2)*f(n_max-3) )
                    rhs(n_max-1)= ( B(2)*f(2) + B(1)*f(1) + B(-1)*f(n_max-1) + B(-2)*f(n_max-2) )

                    call TS_Pr(A(-1,1:n_max-1), A(0,1:n_max-1), A(1,1:n_max-1), rhs(1:n_max-1), d1f(1:n_max-1), n_max-1)

                case (Dirichlet)

                    rhs(1)= ( B_lim(-1)*f(1) + B_lim(0) * f(2) + B_lim(1) * f(3) )

                    rhs(n_max-2)= ( B(2)*f(n_max) + B(1)*f(n_max-1) + B(-1)*f(n_max-2) + B(-2)*f(n_max-3) )
                    rhs(n_max-1)= -( + B_lim(-1)*f(n_max) + B_lim(0) * f(n_max-1) + B_lim(1)*f(n_max-2) )

                    call TS_NPr(rhs(1:n_max-1), d1f(1:n_max-1), n_max-1)

                case (antisymetric)

                    ! TOCOMPLETE

                    rhs(1)= ( B(2)*f(3) + B(1)*f(2) + B(-1)*f(1) + B(-2)*(2.d0*f(1)-f(2)) )

                    rhs(n_max-2)= ( B(2)*f(n_max) + B(1)*f(n_max-1) + B(-1)*f(n_max-2) + B(-2)*f(n_max-3) )
                    rhs(n_max-1)= ( B(2)*(2.d0*f(n_max)-f(n_max-1)) + B(1)*f(n_max) + B(-1)*f(n_max-1) + B(-2)*f(n_max-2) )

                    call TS_NPr(rhs(1:n_max-1), d1f(1:n_max-1), n_max-1)

            end select

        end if

        return

    contains

        subroutine fill_matrix()
            implicit none

            if (allocated(A)) then
                deallocate(A)
            end if

            allocate(A(-1:1, n_max))


            al=9.d0/62.d0
            B(1)=(3*(3-2*al))/(8*h)
            B(2)= (22*al-1)/(24*h)
            B(-1)=-B(1)
            B(-2)=-B(2)

            B_lim(-1)=-1.d0               / h
            B_lim(0)= 2.d0                / h
            B_lim(1)=-1.d0                / h

            if (shifted) then

                ! On détermine la matrice A ainsi que les RHS

                do i= 3, n_max-2
                    A(-1,i) = al
                    A(0,i) = 1.d0
                    A(1,i) = al
                enddo

                select case (bc_type)
                    case (periodic)

                        A(-1,n_max-2:n_max-1) = al
                        A(0,n_max-2:n_max-1) = 1.d0
                        A(1,n_max-2:n_max-1) = al

                        A(-1,1:2) = al
                        A(0,1:2) = 1.d0
                        A(1,1:2) = al

                    case (Dirichlet)
                        A(-1,2)=0.d0
                        A(0,2)=1.d0
                        A(1,2)=-1.d0

                        A(-1,n_max-1)=-1.d0
                        A(0,n_max-1)=1.d0
                        A(1,n_max-1)=0.d0

                        call TS_init(A(-1,2:n_max-1), A(0,2:n_max-1), A(1,2:n_max-1), n_max-2)

                    case (symetric)
                        A(-1,2)=0.d0
                        A(0,2)=1.d0
                        A(1,2)=al

                        A(-1,n_max-1)=al
                        A(0,n_max-1)=1.d0
                        A(1,n_max-1)=0.d0

                        call TS_init(A(-1,2:n_max-1), A(0,2:n_max-1), A(1,2:n_max-1), n_max-2)

                end select

            else

                ! On détermine la matrice A ainsi que les RHS

                do i= 2, n_max-2
                    A(-1,i) = al
                    A(0,i) = 1.d0
                    A(1,i) = al
                enddo

                select case (bc_type)

                    case (periodic)

                        A(-1,n_max-1) = al
                        A(0,n_max-1) = 1.d0
                        A(1,n_max-1) = al

                        A(-1,1) = al
                        A(0,1) = 1.d0
                        A(1,1) = al

                    case (Dirichlet)
                        A(-1,n_max-1)=-1.d0
                        A(0,n_max-1)=1.d0

                        A(0,1)=1.d0
                        A(1,1)=-1.d0

                        call TS_init(A(-1,1:n_max-1), A(0,1:n_max-1), A(1,1:n_max-1), n_max-1)

                    case (antisymetric)
                        A(-1,n_max-1)=al
                        A(0,n_max-1)=1.d0+al

                        A(0,1)=1.d0+al
                        A(1,1)=al

                        call TS_init(A(-1,1:n_max-1), A(0,1:n_max-1), A(1,1:n_max-1), n_max-1)

                end select

            end if

        end subroutine fill_matrix

    end subroutine



    subroutine D1s_CptCtr_O6Fp0_ACC(f, d1f, n_max, h, shifted, bc_type)

        implicit none
        logical, intent(in) :: shifted
        integer, intent(in) :: n_max, bc_type
        real(kind=8) , intent(in) :: h, f(:)
        real(kind=8), intent(inout) :: d1f(:)

        real(kind=8), dimension(n_max) :: d1f_tmp
        integer :: i

        select case (bc_type)

            case (periodic)

                d1f_tmp(1:n_max-1)=d1f(1:n_max-1)
                call D1s_CptCtr_O6Fp0(f, d1f_tmp, n_max, h, shifted, bc_type)
                d1f(1:n_max-1) = d1f(1:n_max-1) + d1f_tmp(1:n_max-1)

            case (Dirichlet)

                if (shifted) then

                    d1f_tmp(2:n_max-1)=d1f(2:n_max-1)
                    call D1s_CptCtr_O6Fp0(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(2:n_max-1) = d1f(2:n_max-1) + d1f_tmp(2:n_max-1)

                else

                    d1f_tmp(1:n_max-1)=d1f(1:n_max-1)
                    call D1s_CptCtr_O6Fp0(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(1:n_max-1) = d1f(1:n_max-1) + d1f_tmp(1:n_max-1)
                end if

            case (symetric)

                if (shifted) then

                    d1f_tmp(2:n_max-1)=d1f(2:n_max-1)
                    call D1s_CptCtr_O6Fp0(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(2:n_max-1) = d1f(2:n_max-1) + d1f_tmp(2:n_max-1)
                end if

            case (antisymetric)

                if (.not. shifted) then

                    d1f_tmp(1:n_max-1)=d1f(1:n_max-1)
                    call D1s_CptCtr_O6Fp0(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(1:n_max-1) = d1f(1:n_max-1) + d1f_tmp(1:n_max-1)
                end if

            case default

        end select

        return

    end subroutine

    subroutine D1s_CptCtr_O6Fp0_MULT(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        logical, intent(in) :: shifted
        integer, intent(in) :: n_max, bc_type
        real(kind=8) , intent(in) :: h, geom_coefs(:), f(:)
        real(kind=8), intent(out) :: d1f(:)

        integer :: i

        call D1s_CptCtr_O6Fp0(f, d1f, n_max, h, shifted, bc_type)

        select case (bc_type)

            case (periodic)

                do i = 1, n_max-1
                    d1f(i)=d1f(i)*geom_coefs(i)
                end do

            case (Dirichlet)
                if (shifted) then
                    d1f(2:n_max-1)=d1f(2:n_max-1)*geom_coefs(2:n_max-1)
                else
                    d1f(1:n_max-1)=d1f(1:n_max-1)*geom_coefs(1:n_max-1)
                end if

            case (symetric)
                if (shifted) then
                    d1f(2:n_max-1)=d1f(2:n_max-1)*geom_coefs(2:n_max-1)
                end if

            case (antisymetric)
                if (.not. shifted) then
                    d1f(1:n_max-1)=d1f(1:n_max-1)*geom_coefs(1:n_max-1)
                end if

            case default

        end select

        return

    end subroutine
    subroutine D1s_CptCtr_O6Fp0_MULT_ACC(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        logical, intent(in) :: shifted
        integer, intent(in) :: n_max, bc_type
        real(kind=8) , intent(in) :: h, geom_coefs(:), f(:)
        real(kind=8), intent(inout) :: d1f(:)

        real(kind=8), dimension(n_max) :: d1f_tmp
        integer :: i

        select case (bc_type)

            case (periodic)

                d1f_tmp(1:n_max-1)=d1f(1:n_max-1)
                call D1s_CptCtr_O6Fp0(f, d1f_tmp, n_max, h, shifted, bc_type)
                d1f(1:n_max-1) = d1f(1:n_max-1) + d1f_tmp(1:n_max-1)*geom_coefs(1:n_max-1)

            case (Dirichlet)

                if (shifted) then

                    d1f_tmp(2:n_max-1)=d1f(2:n_max-1)
                    call D1s_CptCtr_O6Fp0(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(2:n_max-1) = d1f(2:n_max-1) + d1f_tmp(2:n_max-1)*geom_coefs(2:n_max-1)

                else

                    d1f_tmp(1:n_max-1)=d1f(1:n_max-1)
                    call D1s_CptCtr_O6Fp0(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(1:n_max-1) = d1f(1:n_max-1) + d1f_tmp(1:n_max-1)*geom_coefs(1:n_max-1)
                end if


            case (symetric)

                if (shifted) then

                    d1f_tmp(2:n_max-1)=d1f(2:n_max-1)
                    call D1s_CptCtr_O6Fp0(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(2:n_max-1) = d1f(2:n_max-1) + d1f_tmp(2:n_max-1)*geom_coefs(2:n_max-1)

                end if

            case (antisymetric)

                if (.not. shifted) then

                    d1f_tmp(1:n_max-1)=d1f(1:n_max-1)
                    call D1s_CptCtr_O6Fp0(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(1:n_max-1) = d1f(1:n_max-1) + d1f_tmp(1:n_max-1)*geom_coefs(1:n_max-1)
                end if

            case default

        end select

        return

    end subroutine

    subroutine D1s_CTR_CLOSURE_INF(f, d1f, h, nb_ghost, shifted)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(nb_ghost)

        implicit none

        integer :: i, n_max, nb_ghost
        logical ::  shifted
        real(kind=8) h
        real(kind=8) d1f(:)
        real(kind=8) f(:)

        real(kind=8) d1f_temp(CPT_MIN)

        if (nb_ghost.gt.5) then
            call exit(1)
        end if

        if (shifted) then


            if (treat_boundary_by_cpt_scheme) then

                call D1s_CptCtr_O6Fp0(f(1:CPT_MIN), d1f_temp(1:CPT_MIN), CPT_MIN, h, .true., Dirichlet)

                d1f(2)=d1f_temp(2)
                if (nb_ghost.eq.2) return

                d1f(3)=d1f_temp(3)
                if (nb_ghost.eq.3) return

            else

                ! i=2---------------------------------------------------
                i=2
                d1f(i)= (f(i)-f(i-1))/h
                ! i=3---------------------------------------------------
                i=3
                d1f(i)  = 9.0d0/8.d0  * (f(i)-f(i-1))     / h       &
                -   1.d0/24.d0  * (f(i+1)-f(i-2))           / h

                if (nb_ghost.eq.i) return

            endif

            ! i=4---------------------------------------------------
            i=4
            d1f(i)  = 1.189101951200031d0  *   (f(i)-f(i-1))  / h       &
            -   0.073717642266682d0  *   (f(i+1)-f(i-2))        / h       &
            +   0.006410195120003d0 *   (f(i+2) - f(i-3))       / h

            if (nb_ghost.eq.i) return

            ! i=5---------------------------------------------------
            i=5
            d1f(i)  = 1.189101951200031d0  *   (f(i)-f(i-1))  / h         &
            -   0.073717642266682d0  *   (f(i+1)-f(i-2))        / h         &
            +   0.006410195120003d0 *   (f(i+2) - f(i-3))       / h

            if (nb_ghost.eq.i) return

        else


            if (treat_boundary_by_cpt_scheme) then

                call D1s_CptCtr_O6Fp0(f(1:CPT_MIN), d1f_temp(1:CPT_MIN), CPT_MIN, h, .false., Dirichlet)

                d1f(1)=d1f_temp(1)
                if (nb_ghost.eq.1) return

                d1f(2)=d1f_temp(2)
                if (nb_ghost.eq.2) return

            else

                ! i=1---------------------------------------------------
                i=1
                d1f(i)= (f(i+1)-f(i))/h
                ! i=2---------------------------------------------------
                i=2
                d1f(i)  =  9.0d0/8.d0  * (f(i+1)-f(i))     / h     &
                -   1.d0/24.d0  * (f(i+2)-f(i-1))           / h


                if (nb_ghost.eq.i) return

            endif

            ! i=3---------------------------------------------------
            i=3
            d1f(i)  = 1.189101951200031d0  *   (f(i+1)-f(i))  / h         &
            -   0.073717642266682d0  *   (f(i+2)-f(i-1))        / h         &
            +   0.006410195120003d0 *   (f(i+3) - f(i-2))       / h

            if (nb_ghost.eq.i) return

            ! i=4---------------------------------------------------
            i=4
            d1f(i)  = 1.189101951200031d0  *   (f(i+1)-f(i))  / h         &
            -   0.073717642266682d0  *   (f(i+2)-f(i-1))        / h         &
            +   0.006410195120003d0 *   (f(i+3) - f(i-2))       / h

            if (nb_ghost.eq.i) return

            ! i=5---SEE D0s_5pts.wxm-------wc=1.8-------------------
            i=5
            d1f(i)  = 1.23519718551682d0  *   (f(i+1)-f(i))   / h     &
            -   0.10695865844795d0  *   (f(i+2)-f(i-1))         / h     &
            +   0.022474471128913d0 *   (f(i+3) - f(i-2))       / h     &
            -   0.0045027791193291d0 *   (f(i+4) - f(i-3))      / h     &
            +   5.3624403051412786d-4*   (f(i+5) - f(i-4))      / h

        end if

        return

    end subroutine

    subroutine D1s_CTR_CLOSURE_INF_ACC(f, d1f, h, nb_ghost, shifted)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(nb_ghost)

        implicit none

        integer :: i, n_max, nb_ghost
        logical ::  shifted
        real(kind=8) h
        real(kind=8) d1f(:)
        real(kind=8) f(:)

        real(kind=8) d1f_temp(CPT_MIN)

        if (nb_ghost.gt.5) then
            call exit(1)
        end if

        if (shifted) then


            if (treat_boundary_by_cpt_scheme) then

                call D1s_CptCtr_O6Fp0(f(1:CPT_MIN), d1f_temp(1:CPT_MIN), CPT_MIN, h, .true., Dirichlet)

                d1f(2)=d1f(2) + d1f_temp(2)
                if (nb_ghost.eq.2) return

                d1f(3)=d1f(3) + d1f_temp(3)
                if (nb_ghost.eq.3) return

            else

                ! i=2---------------------------------------------------
                i=2
                d1f(i)=d1f(i) + (f(i)-f(i-1))/h
                ! i=3---------------------------------------------------
                i=3
                d1f(i)  = d1f(i) + 9.0d0/8.d0  * (f(i)-f(i-1))     / h       &
                -   1.d0/24.d0  * (f(i+1)-f(i-2))           / h

                if (nb_ghost.eq.i) return

            endif

            ! i=4---------------------------------------------------
            i=4
            d1f(i)  = d1f(i) + 1.189101951200031d0  *   (f(i)-f(i-1))  / h       &
            -   0.073717642266682d0  *   (f(i+1)-f(i-2))        / h       &
            +   0.006410195120003d0 *   (f(i+2) - f(i-3))       / h

            if (nb_ghost.eq.i) return

            ! i=5---------------------------------------------------
            i=5
            d1f(i)  = d1f(i) + 1.189101951200031d0  *   (f(i)-f(i-1))  / h         &
            -   0.073717642266682d0  *   (f(i+1)-f(i-2))        / h         &
            +   0.006410195120003d0 *   (f(i+2) - f(i-3))       / h

            if (nb_ghost.eq.i) return

        else


            if (treat_boundary_by_cpt_scheme) then

                call D1s_CptCtr_O6Fp0(f(1:CPT_MIN), d1f_temp(1:CPT_MIN), CPT_MIN, h, .false., Dirichlet)

                d1f(1)=d1f(1) + d1f_temp(1)
                if (nb_ghost.eq.1) return

                d1f(2)=d1f(2) + d1f_temp(2)
                if (nb_ghost.eq.2) return

            else

                ! i=1---------------------------------------------------
                i=1
                d1f(i)= d1f(i) + (f(i+1)-f(i))/h
                ! i=2---------------------------------------------------
                i=2
                d1f(i)  = d1f(i) + 9.0d0/8.d0  * (f(i+1)-f(i))     / h     &
                -   1.d0/24.d0  * (f(i+2)-f(i-1))           / h


                if (nb_ghost.eq.i) return

            endif

            ! i=3---------------------------------------------------
            i=3
            d1f(i)  = d1f(i) + 1.189101951200031d0  *   (f(i+1)-f(i))  / h         &
            -   0.073717642266682d0  *   (f(i+2)-f(i-1))        / h         &
            +   0.006410195120003d0 *   (f(i+3) - f(i-2))       / h

            if (nb_ghost.eq.i) return

            ! i=4---------------------------------------------------
            i=4
            d1f(i)  = d1f(i) + 1.189101951200031d0  *   (f(i+1)-f(i))  / h         &
            -   0.073717642266682d0  *   (f(i+2)-f(i-1))        / h         &
            +   0.006410195120003d0 *   (f(i+3) - f(i-2))       / h

            if (nb_ghost.eq.i) return

            ! i=5---SEE D0s_5pts.wxm-------wc=1.8-------------------
            i=5
            d1f(i)  = d1f(i) + 1.23519718551682d0  *   (f(i+1)-f(i))   / h     &
            -   0.10695865844795d0  *   (f(i+2)-f(i-1))         / h     &
            +   0.022474471128913d0 *   (f(i+3) - f(i-2))       / h     &
            -   0.0045027791193291d0 *   (f(i+4) - f(i-3))      / h     &
            +   5.3624403051412786d-4*   (f(i+5) - f(i-4))      / h

        end if

        return

    end subroutine

    subroutine D1s_CTR_CLOSURE_INF_MULT(f, d1f, h, nb_ghost, shifted, geom_coef)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(nb_ghost)

        implicit none

        integer :: i, n_max, nb_ghost
        logical ::  shifted
        real(kind=8) h
        real(kind=8) d1f(:), geom_coef(:)
        real(kind=8) f(:)

        real(kind=8) d1f_temp(CPT_MIN)

        if (nb_ghost.gt.5) then
            call exit(1)
        end if

        if (shifted) then


            if (treat_boundary_by_cpt_scheme) then

                call D1s_CptCtr_O6Fp0(f(1:CPT_MIN), d1f_temp(1:CPT_MIN), CPT_MIN, h, .true., Dirichlet)

                d1f(2)=d1f_temp(2) * geom_coef(2)
                if (nb_ghost.eq.2) return

                d1f(3)=d1f_temp(3) * geom_coef(3)
                if (nb_ghost.eq.3) return

            else

                ! i=2---------------------------------------------------
                i=2
                d1f(i)=( (f(i)-f(i-1))/h ) * geom_coef(i)
                ! i=3---------------------------------------------------
                i=3
                d1f(i)  = ( 9.0d0/8.d0  * (f(i)-f(i-1))             / h       &
                -   1.d0/24.d0  * (f(i+1)-f(i-2))                   / h ) * geom_coef(i)

                if (nb_ghost.eq.i) return

            endif

            ! i=4---------------------------------------------------
            i=4
            d1f(i)  = ( 1.189101951200031d0  *   (f(i)-f(i-1))          / h       &
            -   0.073717642266682d0  *   (f(i+1)-f(i-2))                / h       &
            +   0.006410195120003d0 *   (f(i+2) - f(i-3))               / h ) * geom_coef(i)

            if (nb_ghost.eq.i) return

            ! i=5---------------------------------------------------
            i=5
            d1f(i)  = ( 1.189101951200031d0  *   (f(i)-f(i-1))              / h         &
            -   0.073717642266682d0  *   (f(i+1)-f(i-2))                    / h         &
            +   0.006410195120003d0 *   (f(i+2) - f(i-3))                   / h ) * geom_coef(i)

            if (nb_ghost.eq.i) return

        else


            if (treat_boundary_by_cpt_scheme) then

                call D1s_CptCtr_O6Fp0(f(1:CPT_MIN), d1f_temp(1:CPT_MIN), CPT_MIN, h, .false., Dirichlet)

                d1f(1)=d1f_temp(1) * geom_coef(1)
                if (nb_ghost.eq.1) return

                d1f(2)=d1f_temp(2) * geom_coef(2)
                if (nb_ghost.eq.2) return

            else

                ! i=1---------------------------------------------------
                i=1
                d1f(i)= ( (f(i+1)-f(i))/h ) * geom_coef(i)
                ! i=2---------------------------------------------------
                i=2
                d1f(i)  = ( 9.0d0/8.d0  * (f(i+1)-f(i))             / h     &
                -   1.d0/24.d0  * (f(i+2)-f(i-1))                   / h ) * geom_coef(i)


                if (nb_ghost.eq.i) return

            endif

            ! i=3---------------------------------------------------
            i=3
            d1f(i)  = ( 1.189101951200031d0  *   (f(i+1)-f(i))     / h         &
            -   0.073717642266682d0  *   (f(i+2)-f(i-1))                    / h         &
            +   0.006410195120003d0 *   (f(i+3) - f(i-2))                   / h ) * geom_coef(i)

            if (nb_ghost.eq.i) return

            ! i=4---------------------------------------------------
            i=4
            d1f(i)  = ( 1.189101951200031d0  *   (f(i+1)-f(i)) / h         &
            -   0.073717642266682d0  *   (f(i+2)-f(i-1))                / h         &
            +   0.006410195120003d0 *   (f(i+3) - f(i-2))               / h ) * geom_coef(i)

            if (nb_ghost.eq.i) return

            ! i=5---SEE D0s_5pts.wxm-------wc=1.8-------------------
            i=5
            d1f(i)  = ( 1.23519718551682d0  *   (f(i+1)-f(i))           / h     &
            -   0.10695865844795d0  *   (f(i+2)-f(i-1))                 / h     &
            +   0.022474471128913d0 *   (f(i+3) - f(i-2))               / h     &
            -   0.0045027791193291d0 *   (f(i+4) - f(i-3))              / h     &
            +   5.3624403051412786d-4*   (f(i+5) - f(i-4))              / h ) * geom_coef(i)

        end if

        return

    end subroutine

    subroutine D1s_CTR_CLOSURE_INF_MULT_ACC(f, d1f, h, nb_ghost, shifted, geom_coef)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(nb_ghost)

        implicit none

        integer :: i, n_max, nb_ghost
        logical ::  shifted
        real(kind=8) h
        real(kind=8) d1f(:), geom_coef(:)
        real(kind=8) f(:)

        real(kind=8) d1f_temp(CPT_MIN)

        if (nb_ghost.gt.5) then
            call exit(1)
        end if

        if (shifted) then


            if (treat_boundary_by_cpt_scheme) then

                call D1s_CptCtr_O6Fp0(f(1:CPT_MIN), d1f_temp(1:CPT_MIN), CPT_MIN, h, .true., Dirichlet)

                d1f(2)=d1f(2) + d1f_temp(2) * geom_coef(2)
                if (nb_ghost.eq.2) return

                d1f(3)=d1f(3) + d1f_temp(3) * geom_coef(3)
                if (nb_ghost.eq.3) return

            else

                ! i=2---------------------------------------------------
                i=2
                d1f(i)=d1f(i) + ( (f(i)-f(i-1))/h ) * geom_coef(i)
                ! i=3---------------------------------------------------
                i=3
                d1f(i)  = d1f(i) + ( 9.0d0/8.d0  * (f(i)-f(i-1))    / h       &
                -   1.d0/24.d0  * (f(i+1)-f(i-2))                   / h ) * geom_coef(i)

                if (nb_ghost.eq.i) return

            endif

            ! i=4---------------------------------------------------
            i=4
            d1f(i)  = d1f(i) + ( 1.189101951200031d0  *   (f(i)-f(i-1))  / h       &
            -   0.073717642266682d0  *   (f(i+1)-f(i-2))                / h       &
            +   0.006410195120003d0 *   (f(i+2) - f(i-3))               / h ) * geom_coef(i)

            if (nb_ghost.eq.i) return

            ! i=5---------------------------------------------------
            i=5
            d1f(i)  = d1f(i) + ( 1.189101951200031d0  *   (f(i)-f(i-1))     / h         &
            -   0.073717642266682d0  *   (f(i+1)-f(i-2))                    / h         &
            +   0.006410195120003d0 *   (f(i+2) - f(i-3))                   / h ) * geom_coef(i)

            if (nb_ghost.eq.i) return

        else


            if (treat_boundary_by_cpt_scheme) then

                call D1s_CptCtr_O6Fp0(f(1:CPT_MIN), d1f_temp(1:CPT_MIN), CPT_MIN, h, .false., Dirichlet)

                d1f(1)=d1f(1) + d1f_temp(1) * geom_coef(1)
                if (nb_ghost.eq.1) return

                d1f(2)=d1f(2) + d1f_temp(2) * geom_coef(2)
                if (nb_ghost.eq.2) return

            else

                ! i=1---------------------------------------------------
                i=1
                d1f(i)= d1f(i) + ( (f(i+1)-f(i))/h ) * geom_coef(i)
                ! i=2---------------------------------------------------
                i=2
                d1f(i)  = d1f(i) + ( 9.0d0/8.d0  * (f(i+1)-f(i))    / h     &
                -   1.d0/24.d0  * (f(i+2)-f(i-1))                   / h ) * geom_coef(i)


                if (nb_ghost.eq.i) return

            endif

            ! i=3---------------------------------------------------
            i=3
            d1f(i)  = d1f(i) + ( 1.189101951200031d0  *   (f(i+1)-f(i))     / h         &
            -   0.073717642266682d0  *   (f(i+2)-f(i-1))                    / h         &
            +   0.006410195120003d0 *   (f(i+3) - f(i-2))                   / h ) * geom_coef(i)

            if (nb_ghost.eq.i) return

            ! i=4---------------------------------------------------
            i=4
            d1f(i)  = d1f(i) + ( 1.189101951200031d0  *   (f(i+1)-f(i)) / h         &
            -   0.073717642266682d0  *   (f(i+2)-f(i-1))                / h         &
            +   0.006410195120003d0 *   (f(i+3) - f(i-2))               / h ) * geom_coef(i)

            if (nb_ghost.eq.i) return

            ! i=5---SEE D0s_5pts.wxm-------wc=1.8-------------------
            i=5
            d1f(i)  = d1f(i) + ( 1.23519718551682d0  *   (f(i+1)-f(i))  / h     &
            -   0.10695865844795d0  *   (f(i+2)-f(i-1))                 / h     &
            +   0.022474471128913d0 *   (f(i+3) - f(i-2))               / h     &
            -   0.0045027791193291d0 *   (f(i+4) - f(i-3))              / h     &
            +   5.3624403051412786d-4*   (f(i+5) - f(i-4))              / h ) * geom_coef(i)

        end if

        return

    end subroutine

    subroutine D1s_CTR_CLOSURE_SUP(f, d1f, n_max, h, nb_ghost, shifted)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(nb_ghost)

        implicit none

        integer :: i, n_max, nb_ghost, j
        real(kind=8) h, fn
        real(kind=8) d1f(:)
        real(kind=8) f(:)
        logical :: shifted

        real(kind=8) d1f_temp(n_max-CPT_MIN:n_max-1)

        if (nb_ghost.gt.5) then
            call exit(1)
        end if

        if (shifted) then

            if (treat_boundary_by_cpt_scheme) then

                !   ---f------d0f-----f------d0f-------f
                call D1s_CptCtr_O6Fp0(f(n_max-CPT_MIN-1:n_max-1), d1f_temp, CPT_MIN+1, h, .false., Dirichlet)

                d1f(n_max-1)=d1f_temp(n_max-1)

                d1f(n_max-2)=d1f_temp(n_max-2)
                if (nb_ghost.eq.2) return

            else

                ! First point-------------------------------------------
                i=n_max-1
                d1f(i)= (f(i)-f(i-1))    / h
                ! 2nd point----------------------------------------------
                i=n_max-2
                d1f(i)  =  9.0d0/8.d0 * (f(i)-f(i-1))  / h     &
                -   1.d0/24.d0    * (f(i+1)-f(i-2))     / h


                if (nb_ghost.eq.2) return

            endif

            ! 3rd point---------------------------------------------
            i=n_max-3
            d1f(i)  = 1.189101951200031d0  *   (f(i)-f(i-1))  / h       &
            -   0.073717642266682d0  *   (f(i+1)-f(i-2))        / h       &
            +   0.006410195120003d0 *   (f(i+2) - f(i-3))       / h

            if (nb_ghost.eq.3) return

            ! 4th point---------------------------------------------
            i=n_max-4
            d1f(i)  = 1.189101951200031d0  *   (f(i)-f(i-1))  / h     &
            -   0.073717642266682d0  *   (f(i+1)-f(i-2))        / h     &
            +   0.006410195120003d0 *   (f(i+2) - f(i-3))       / h

            if (nb_ghost.eq.4) return

        else

            if (treat_boundary_by_cpt_scheme) then

                !   ---f------d0f-----f------d0f-------f
                call D1s_CptCtr_O6Fp0(f(n_max-CPT_MIN:n_max), d1f_temp, CPT_MIN+1, h, .false., Dirichlet)

                d1f(n_max-1)=d1f_temp(n_max-1)

                d1f(n_max-2)=d1f_temp(n_max-2)
                if (nb_ghost.eq.2) return

            else

                ! First point-------------------------------------------
                i=n_max-1
                d1f(i)= (f(i+1)-f(i))    / h
                ! 2nd point----------------------------------------------
                i=n_max-2
                d1f(i)  = 9.0d0/8.d0 * (f(i+1)-f(i))  / h     &
                -   1.d0/24.d0    * (f(i+2)-f(i-1))     / h


                if (nb_ghost.eq.2) return

            endif

            ! 3rd point---------------------------------------------
            i=n_max-3
            d1f(i)  = 1.189101951200031d0  *   (f(i+1)-f(i))  / h         &
            -   0.073717642266682d0  *   (f(i+2)-f(i-1))        / h         &
            +   0.006410195120003d0 *   (f(i+3) - f(i-2))       / h

            if (nb_ghost.eq.3) return

            ! 4th point---------------------------------------------
            i=n_max-4
            d1f(i)  = 1.189101951200031d0  *   (f(i+1)-f(i))  / h     &
            -   0.073717642266682d0  *   (f(i+2)-f(i-1))        / h     &
            +   0.006410195120003d0 *   (f(i+3) - f(i-2))       / h

            if (nb_ghost.eq.4) return

            ! 5th point ----D0s_5pts.wxm-------wc=1.8---------------
            i=n_max-5
            d1f(i)  = 1.23519718551682d0  *   (f(i+1)-f(i))   / h         &
            -   0.10695865844795d0  *   (f(i+2)-f(i-1))         / h         &
            +   0.022474471128913d0 *   (f(i+3) - f(i-2))       / h         &
            -   0.0045027791193291d0 *   (f(i+4) - f(i-3))      / h         &
            +   5.3624403051412786d-4*   (f(i+5) - f(i-4))      / h

        end if

        return

    end subroutine

    subroutine D1s_CTR_CLOSURE_SUP_ACC(f, d1f, n_max, h, nb_ghost, shifted)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(nb_ghost)

        implicit none

        integer :: i, n_max, nb_ghost, j
        real(kind=8) h, fn
        real(kind=8) d1f(:)
        real(kind=8) f(:)
        logical :: shifted

        real(kind=8) d1f_temp(n_max-CPT_MIN:n_max-1)

        if (nb_ghost.gt.5) then
            call exit(1)
        end if

        if (shifted) then

            if (treat_boundary_by_cpt_scheme) then

                !   ---f------d0f-----f------d0f-------f
                call D1s_CptCtr_O6Fp0(f(n_max-CPT_MIN-1:n_max-1), d1f_temp, CPT_MIN+1, h, .false., Dirichlet)

                d1f(n_max-1)=d1f(n_max-1) + d1f_temp(n_max-1)

                d1f(n_max-2)=d1f(n_max-2) + d1f_temp(n_max-2)
                if (nb_ghost.eq.2) return

            else

                ! First point-------------------------------------------
                i=n_max-1
                d1f(i)= d1f(i) + (f(i)-f(i-1))    / h
                ! 2nd point----------------------------------------------
                i=n_max-2
                d1f(i)  = d1f(i) +  9.0d0/8.d0 * (f(i)-f(i-1))  / h     &
                -   1.d0/24.d0    * (f(i+1)-f(i-2))     / h


                if (nb_ghost.eq.2) return

            endif

            ! 3rd point---------------------------------------------
            i=n_max-3
            d1f(i)  = d1f(i) + 1.189101951200031d0  *   (f(i)-f(i-1))  / h       &
            -   0.073717642266682d0  *   (f(i+1)-f(i-2))        / h       &
            +   0.006410195120003d0 *   (f(i+2) - f(i-3))       / h

            if (nb_ghost.eq.3) return

            ! 4th point---------------------------------------------
            i=n_max-4
            d1f(i)  = d1f(i) + 1.189101951200031d0  *   (f(i)-f(i-1))  / h     &
            -   0.073717642266682d0  *   (f(i+1)-f(i-2))        / h     &
            +   0.006410195120003d0 *   (f(i+2) - f(i-3))       / h

            if (nb_ghost.eq.4) return

        else

            if (treat_boundary_by_cpt_scheme) then

                !   ---f------d0f-----f------d0f-------f
                call D1s_CptCtr_O6Fp0(f(n_max-CPT_MIN:n_max), d1f_temp, CPT_MIN+1, h, .false., Dirichlet)

                d1f(n_max-1)=d1f(n_max-1) + d1f_temp(n_max-1)

                d1f(n_max-2)=d1f(n_max-2) + d1f_temp(n_max-2)
                if (nb_ghost.eq.2) return

            else

                ! First point-------------------------------------------
                i=n_max-1
                d1f(i)= d1f(i) + (f(i+1)-f(i))    / h
                ! 2nd point----------------------------------------------
                i=n_max-2
                d1f(i)  = d1f(i) + 9.0d0/8.d0 * (f(i+1)-f(i))  / h     &
                -   1.d0/24.d0    * (f(i+2)-f(i-1))     / h


                if (nb_ghost.eq.2) return

            endif

            ! 3rd point---------------------------------------------
            i=n_max-3
            d1f(i)  = d1f(i) + 1.189101951200031d0  *   (f(i+1)-f(i))  / h         &
            -   0.073717642266682d0  *   (f(i+2)-f(i-1))        / h         &
            +   0.006410195120003d0 *   (f(i+3) - f(i-2))       / h

            if (nb_ghost.eq.3) return

            ! 4th point---------------------------------------------
            i=n_max-4
            d1f(i)  = d1f(i) + 1.189101951200031d0  *   (f(i+1)-f(i))  / h     &
            -   0.073717642266682d0  *   (f(i+2)-f(i-1))        / h     &
            +   0.006410195120003d0 *   (f(i+3) - f(i-2))       / h

            if (nb_ghost.eq.4) return

            ! 5th point ----D0s_5pts.wxm-------wc=1.8---------------
            i=n_max-5
            d1f(i)  = d1f(i) + 1.23519718551682d0  *   (f(i+1)-f(i))   / h         &
            -   0.10695865844795d0  *   (f(i+2)-f(i-1))         / h         &
            +   0.022474471128913d0 *   (f(i+3) - f(i-2))       / h         &
            -   0.0045027791193291d0 *   (f(i+4) - f(i-3))      / h         &
            +   5.3624403051412786d-4*   (f(i+5) - f(i-4))      / h

        end if

        return

    end subroutine

    subroutine D1s_CTR_CLOSURE_SUP_MULT(f, d1f, n_max, h, nb_ghost, shifted, geom_coef)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(nb_ghost)

        implicit none

        integer :: i, n_max, nb_ghost, j
        real(kind=8) h, fn
        real(kind=8) d1f(:), geom_coef(:)
        real(kind=8) f(:)
        logical :: shifted

        real(kind=8) d1f_temp(n_max-CPT_MIN:n_max-1)

        if (nb_ghost.gt.5) then
            call exit(1)
        end if

        if (shifted) then

            if (treat_boundary_by_cpt_scheme) then

                !   ---f------d0f-----f------d0f-------f
                call D1s_CptCtr_O6Fp0(f(n_max-CPT_MIN-1:n_max-1), d1f_temp, CPT_MIN+1, h, .false., Dirichlet)

                d1f(n_max-1)=d1f_temp(n_max-1) * geom_coef(n_max-1)

                d1f(n_max-2)=d1f_temp(n_max-2) * geom_coef(n_max-2)
                if (nb_ghost.eq.2) return

            else

                ! First point-------------------------------------------
                i=n_max-1
                d1f(i)= ( (f(i)-f(i-1))    / h ) * geom_coef(i)
                ! 2nd point----------------------------------------------
                i=n_max-2
                d1f(i)  = ( 9.0d0/8.d0 * (f(i)-f(i-1))    / h     &
                -   1.d0/24.d0    * (f(i+1)-f(i-2))                 / h ) * geom_coef(i)


                if (nb_ghost.eq.2) return

            endif

            ! 3rd point---------------------------------------------
            i=n_max-3
            d1f(i)  = ( 1.189101951200031d0  *   (f(i)-f(i-1)) / h       &
            -   0.073717642266682d0  *   (f(i+1)-f(i-2))                / h       &
            +   0.006410195120003d0 *   (f(i+2) - f(i-3))               / h ) * geom_coef(i)

            if (nb_ghost.eq.3) return

            ! 4th point---------------------------------------------
            i=n_max-4
            d1f(i)  = ( 1.189101951200031d0  *   (f(i)-f(i-1)) / h     &
            -   0.073717642266682d0  *   (f(i+1)-f(i-2))                / h     &
            +   0.006410195120003d0 *   (f(i+2) - f(i-3))               / h ) * geom_coef(i)

            if (nb_ghost.eq.4) return

        else

            if (treat_boundary_by_cpt_scheme) then

                !   ---f------d0f-----f------d0f-------f
                call D1s_CptCtr_O6Fp0(f(n_max-CPT_MIN:n_max), d1f_temp, CPT_MIN+1, h, .false., Dirichlet)

                d1f(n_max-1)=d1f_temp(n_max-1) * geom_coef(n_max-1)

                d1f(n_max-2)=d1f_temp(n_max-2) * geom_coef(n_max-2)
                if (nb_ghost.eq.2) return

            else

                ! First point-------------------------------------------
                i=n_max-1
                d1f(i)= ( (f(i+1)-f(i))    / h ) * geom_coef(i)
                ! 2nd point----------------------------------------------
                i=n_max-2
                d1f(i)  = ( 9.0d0/8.d0 * (f(i+1)-f(i)) / h     &
                -   1.d0/24.d0    * (f(i+2)-f(i-1))             / h ) * geom_coef(i)


                if (nb_ghost.eq.2) return

            endif

            ! 3rd point---------------------------------------------
            i=n_max-3
            d1f(i)  = ( 1.189101951200031d0  *   (f(i+1)-f(i)) / h         &
            -   0.073717642266682d0  *   (f(i+2)-f(i-1))                / h         &
            +   0.006410195120003d0 *   (f(i+3) - f(i-2))               / h ) * geom_coef(i)

            if (nb_ghost.eq.3) return

            ! 4th point---------------------------------------------
            i=n_max-4
            d1f(i)  = ( 1.189101951200031d0  *   (f(i+1)-f(i)) / h     &
            -   0.073717642266682d0  *   (f(i+2)-f(i-1))                / h     &
            +   0.006410195120003d0 *   (f(i+3) - f(i-2))               / h ) * geom_coef(i)

            if (nb_ghost.eq.4) return

            ! 5th point ----D0s_5pts.wxm-------wc=1.8---------------
            i=n_max-5
            d1f(i)  = ( 1.23519718551682d0  *   (f(i+1)-f(i))  / h         &
            -   0.10695865844795d0  *   (f(i+2)-f(i-1))                 / h         &
            +   0.022474471128913d0 *   (f(i+3) - f(i-2))               / h         &
            -   0.0045027791193291d0 *   (f(i+4) - f(i-3))              / h         &
            +   5.3624403051412786d-4*   (f(i+5) - f(i-4))              / h ) * geom_coef(i)

        end if

        return

    end subroutine

    subroutine D1s_CTR_CLOSURE_SUP_MULT_ACC(f, d1f, n_max, h, nb_ghost, shifted, geom_coef)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(nb_ghost)

        implicit none

        integer :: i, n_max, nb_ghost, j
        real(kind=8) h, fn
        real(kind=8) d1f(:), geom_coef(:)
        real(kind=8) f(:)
        logical :: shifted

        real(kind=8) d1f_temp(n_max-CPT_MIN:n_max-1)

        if (nb_ghost.gt.5) then
            call exit(1)
        end if

        if (shifted) then

            if (treat_boundary_by_cpt_scheme) then

                !   ---f------d0f-----f------d0f-------f
                call D1s_CptCtr_O6Fp0(f(n_max-CPT_MIN-1:n_max-1), d1f_temp, CPT_MIN+1, h, .false., Dirichlet)

                d1f(n_max-1)=d1f(n_max-1) + d1f_temp(n_max-1) * geom_coef(n_max-1)

                d1f(n_max-2)=d1f(n_max-2) + d1f_temp(n_max-2) * geom_coef(n_max-2)
                if (nb_ghost.eq.2) return

            else

                ! First point-------------------------------------------
                i=n_max-1
                d1f(i)= d1f(i) + ( (f(i)-f(i-1))    / h ) * geom_coef(i)
                ! 2nd point----------------------------------------------
                i=n_max-2
                d1f(i)  = d1f(i) +  ( 9.0d0/8.d0 * (f(i)-f(i-1))    / h     &
                -   1.d0/24.d0    * (f(i+1)-f(i-2))                 / h ) * geom_coef(i)


                if (nb_ghost.eq.2) return

            endif

            ! 3rd point---------------------------------------------
            i=n_max-3
            d1f(i)  = d1f(i) + ( 1.189101951200031d0  *   (f(i)-f(i-1)) / h       &
            -   0.073717642266682d0  *   (f(i+1)-f(i-2))                / h       &
            +   0.006410195120003d0 *   (f(i+2) - f(i-3))               / h ) * geom_coef(i)

            if (nb_ghost.eq.3) return

            ! 4th point---------------------------------------------
            i=n_max-4
            d1f(i)  = d1f(i) + ( 1.189101951200031d0  *   (f(i)-f(i-1)) / h     &
            -   0.073717642266682d0  *   (f(i+1)-f(i-2))                / h     &
            +   0.006410195120003d0 *   (f(i+2) - f(i-3))               / h ) * geom_coef(i)

            if (nb_ghost.eq.4) return

        else

            if (treat_boundary_by_cpt_scheme) then

                !   ---f------d0f-----f------d0f-------f
                call D1s_CptCtr_O6Fp0(f(n_max-CPT_MIN:n_max), d1f_temp, CPT_MIN+1, h, .false., Dirichlet)

                d1f(n_max-1)=d1f(n_max-1) + d1f_temp(n_max-1) * geom_coef(n_max-1)

                d1f(n_max-2)=d1f(n_max-2) + d1f_temp(n_max-2) * geom_coef(n_max-2)
                if (nb_ghost.eq.2) return

            else

                ! First point-------------------------------------------
                i=n_max-1
                d1f(i)= d1f(i) + ( (f(i+1)-f(i))    / h ) * geom_coef(i)
                ! 2nd point----------------------------------------------
                i=n_max-2
                d1f(i)  = d1f(i) + ( 9.0d0/8.d0 * (f(i+1)-f(i)) / h     &
                -   1.d0/24.d0    * (f(i+2)-f(i-1))             / h ) * geom_coef(i)


                if (nb_ghost.eq.2) return

            endif

            ! 3rd point---------------------------------------------
            i=n_max-3
            d1f(i)  = d1f(i) + ( 1.189101951200031d0  *   (f(i+1)-f(i)) / h         &
            -   0.073717642266682d0  *   (f(i+2)-f(i-1))                / h         &
            +   0.006410195120003d0 *   (f(i+3) - f(i-2))               / h ) * geom_coef(i)

            if (nb_ghost.eq.3) return

            ! 4th point---------------------------------------------
            i=n_max-4
            d1f(i)  = d1f(i) + ( 1.189101951200031d0  *   (f(i+1)-f(i)) / h     &
            -   0.073717642266682d0  *   (f(i+2)-f(i-1))                / h     &
            +   0.006410195120003d0 *   (f(i+3) - f(i-2))               / h ) * geom_coef(i)

            if (nb_ghost.eq.4) return

            ! 5th point ----D0s_5pts.wxm-------wc=1.8---------------
            i=n_max-5
            d1f(i)  = d1f(i) + ( 1.23519718551682d0  *   (f(i+1)-f(i))  / h         &
            -   0.10695865844795d0  *   (f(i+2)-f(i-1))                 / h         &
            +   0.022474471128913d0 *   (f(i+3) - f(i-2))               / h         &
            -   0.0045027791193291d0 *   (f(i+4) - f(i-3))              / h         &
            +   5.3624403051412786d-4*   (f(i+5) - f(i-4))              / h ) * geom_coef(i)

        end if

        return

    end subroutine

end module


