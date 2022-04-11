
module d0c_schemes

    use boundaries_types

end module

module d1c_schemes

    use boundaries_types
    use schemes_settings

    implicit none

contains

    subroutine D1c_ExpCtr_O2Fp0(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d1f(:)

        real(kind=8) A
        integer:: i


        A=0.5d0/h

        if (shifted) then

            do i=2,n_max-2
                d1f(i)= ( f(i+1) - f(i-1) ) * A
            enddo

            if (bc_type.eq.periodic) then
                d1f(1)= ( f(2) - f(n_max-1) ) * A
                d1f(n_max-1)= ( f(1) - f(n_max-2) ) * A
            end if

        else

            do i=2,n_max-2
                d1f(i)= ( f(i+1) - f(i-1) ) * A
            enddo

            if (bc_type.eq.periodic) then
                d1f(1)= ( f(2) - f(n_max-1) ) * A
                d1f(n_max-1)= ( f(1) - f(n_max-2) ) * A
            end if

            if (bc_type.eq.Dirichlet) then
                d1f(n_max-1)= ( f(n_max) - f(n_max-2) ) * A
            end if

        end if

        return

    end subroutine

    subroutine D1c_ExpCtr_O2Fp0_ACC(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(inout)   :: d1f(:)

        real(kind=8) A
        integer:: i


        A=0.5d0/h

        if (shifted) then

            do i=2,n_max-2
                d1f(i)= d1f(i) + ( f(i+1) - f(i-1) ) * A
            enddo

            if (bc_type.eq.periodic) then
                d1f(1)= d1f(1) + ( f(2) - f(n_max-1) ) * A
                d1f(n_max-1)= d1f(n_max-1) + ( f(1) - f(n_max-2) ) * A
            end if

        else

            do i=2,n_max-2
                d1f(i)= d1f(i) + ( f(i+1) - f(i-1) ) * A
            enddo

            if (bc_type.eq.periodic) then
                d1f(1)= d1f(1) + ( f(2) - f(n_max-1) ) * A
                d1f(n_max-1)= d1f(n_max-1) + ( f(1) - f(n_max-2) ) * A
            end if

            if (bc_type.eq.Dirichlet) then
                d1f(n_max-1)= d1f(n_max-1) + ( f(n_max) - f(n_max-2) ) * A
            end if

        end if

        return

    end subroutine

    subroutine D1c_ExpCtr_O2Fp0_MULT(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:), geom_coefs(:)
        real(kind=8), intent(out)     :: d1f(:)

        real(kind=8) A
        integer:: i


        A=0.5d0/h

        if (shifted) then

            do i=2,n_max-2
                d1f(i)= ( f(i+1) - f(i-1) ) * A * geom_coefs(i)
            enddo

            if (bc_type.eq.periodic) then
                d1f(1)= ( f(2) - f(n_max-1) ) * A * geom_coefs(1)
                d1f(n_max-1)= ( f(1) - f(n_max-2) ) * A * geom_coefs(n_max-1)
            end if

        else

            do i=2,n_max-2
                d1f(i)= ( f(i+1) - f(i-1) ) * A * geom_coefs(i)
            enddo

            if (bc_type.eq.periodic) then
                d1f(1)= ( f(2) - f(n_max-1) ) * A * geom_coefs(1)
                d1f(n_max-1)= ( f(1) - f(n_max-2) ) * A * geom_coefs(n_max-1)
            end if

            if (bc_type.eq.Dirichlet) then
                d1f(n_max-1)= ( f(n_max) - f(n_max-2) ) * A * geom_coefs(n_max-1)
            end if

        end if

        return

    end subroutine

    subroutine D1c_ExpCtr_O2Fp0_MULT_ACC(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:), geom_coefs(:)
        real(kind=8), intent(inout)   :: d1f(:)

        real(kind=8) A
        integer:: i


        A=0.5d0/h

        if (shifted) then

            do i=2,n_max-2
                d1f(i)= d1f(i) + ( f(i+1) - f(i-1) ) * A * geom_coefs(i)
            enddo

            if (bc_type.eq.periodic) then
                d1f(1)= d1f(1) + ( f(2) - f(n_max-1) ) * A * geom_coefs(1)
                d1f(n_max-1)= d1f(n_max-1) + ( f(1) - f(n_max-2) ) * A * geom_coefs(n_max-1)
            end if

        else

            do i=2,n_max-2
                d1f(i)= d1f(i) + ( f(i+1) - f(i-1) ) * A * geom_coefs(i)
            enddo

            if (bc_type.eq.periodic) then
                d1f(1)= d1f(1) + ( f(2) - f(n_max-1) ) * A * geom_coefs(1)
                d1f(n_max-1)= d1f(n_max-1) + ( f(1) - f(n_max-2) ) * A * geom_coefs(n_max-1)
            end if

            if (bc_type.eq.Dirichlet) then
                d1f(n_max-1)= d1f(n_max-1) + ( f(n_max) - f(n_max-2) ) * A * geom_coefs(n_max-1)
            end if

        end if

        return

    end subroutine

    subroutine D1_Tamm(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d1f(:)

        integer :: i
        real(kind=8) A(3)

        A(1) =  0.79926643d0 / h
        A(2) =  -0.18941314d0 / h
        A(3) =  0.02651995d0  / h

        do i=4,n_max-4
            d1f(i)=   A(3)*(f(i+3) - f(i-3))    &
            + A(2)*(f(i+2) - f(i-2))    &
            + A(1)*(f(i+1) - f(i-1))
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 3
                d1f(i)=   A(3)*(f(i+3) - f(mod(n_max+i-5,n_max-1)+1))     &
                + A(2)*(f(i+2) - f(mod(n_max+i-4,n_max-1)+1))     &
                + A(1)*(f(i+1) - f(mod(n_max+i-3,n_max-1)+1))
            enddo

            do i=n_max-3,n_max-1
                d1f(i)=   A(3)*(f(mod(i+2,n_max-1)+1) - f(i-3))       &
                + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-2))       &
                + A(1)*(f(mod(i,n_max-1)+1) - f(i-1))
            enddo

        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_INF(f, d1f, h, 3)
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_SUP(f, d1f, n_max, h, n_max-3, shifted)
        end if

        return

    end subroutine

    subroutine D1_Tamm_ACC(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(inout)   :: d1f(:)

        integer :: i
        real(kind=8) A(3)

        A(1) =  0.79926643d0 / h
        A(2) =  -0.18941314d0 / h
        A(3) =  0.02651995d0  / h

        do i=4,n_max-4
            d1f(i)=   d1f(i) + A(3)*(f(i+3) - f(i-3))    &
            + A(2)*(f(i+2) - f(i-2))    &
            + A(1)*(f(i+1) - f(i-1))
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 3
                d1f(i)=   d1f(i) + A(3)*(f(i+3) - f(mod(n_max+i-5,n_max-1)+1))     &
                + A(2)*(f(i+2) - f(mod(n_max+i-4,n_max-1)+1))     &
                + A(1)*(f(i+1) - f(mod(n_max+i-3,n_max-1)+1))
            enddo

            do i=n_max-3,n_max-1
                d1f(i)=   d1f(i) + A(3)*(f(mod(i+2,n_max-1)+1) - f(i-3))       &
                + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-2))       &
                + A(1)*(f(mod(i,n_max-1)+1) - f(i-1))
            enddo

        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_INF_ACC(f, d1f, h, 3)
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_SUP_ACC(f, d1f, n_max, h, n_max-3, shifted)
        end if

        return

    end subroutine

    subroutine D1_Tamm_MULT_ACC(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:), geom_coefs(:)
        real(kind=8), intent(inout)   :: d1f(:)

        integer :: i
        real(kind=8) A(3)

        A(1) =  0.79926643d0 / h
        A(2) =  -0.18941314d0 / h
        A(3) =  0.02651995d0  / h

        do i=4,n_max-4
            d1f(i)=   d1f(i) + ( A(3)*(f(i+3) - f(i-3))    &
            + A(2)*(f(i+2) - f(i-2))    &
            + A(1)*(f(i+1) - f(i-1)) ) * geom_coefs(i)
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 3
                d1f(i)=   d1f(i) + ( A(3)*(f(i+3) - f(mod(n_max+i-5,n_max-1)+1))     &
                + A(2)*(f(i+2) - f(mod(n_max+i-4,n_max-1)+1))     &
                + A(1)*(f(i+1) - f(mod(n_max+i-3,n_max-1)+1)) ) * geom_coefs(i)
            enddo

            do i=n_max-3,n_max-1
                d1f(i)=   d1f(i) + ( A(3)*(f(mod(i+2,n_max-1)+1) - f(i-3))       &
                + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-2))       &
                + A(1)*(f(mod(i,n_max-1)+1) - f(i-1)) ) * geom_coefs(i)
            enddo

        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_INF_MULT_ACC(f, d1f, h, 3, geom_coefs)
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_SUP_MULT_ACC(f, d1f, n_max, h, n_max-3, shifted, geom_coefs)
        end if

        return

    end subroutine

    subroutine D1_Tamm_7pts(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d1f(:)

        integer :: i
        real(kind=8) A(7)

        A(1) =  9.1942501110343045059277722885D-1    /h
        A(2) =  -3.5582959926835268755667642401D-1   /h
        A(3) =  1.5251501608406492469104928679D-1    /h
        A(4) =  -5.9463040829715772666828596899D-2   /h
        A(5) =  1.9010752709508298659849167988D-2   /h
        A(6) = -4.3808649297336481851137000907D-3   /h
        A(7) = 5.3896121868623384659692955878D-4   /h

        do i=8,n_max-8
            d1f(i)= A(7)*(f(i+7) - f(i-7))                  &
            + A(6)*(f(i+6) - f(i-6))     &
            + A(5)*(f(i+5) - f(i-5))     &
            + A(4)*(f(i+4) - f(i-4))     &
            + A(3)*(f(i+3) - f(i-3))     &
            + A(2)*(f(i+2) - f(i-2))     &
            + A(1)*(f(i+1) - f(i-1))
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 7
                d1f(i)=   A(7)*(f(i+7) - f(mod(n_max+i-9,n_max-1)+1))             &
                + A(6)*(f(i+6) - f(mod(n_max+i-8,n_max-1)+1))  &
                + A(5)*(f(i+5) - f(mod(n_max+i-7,n_max-1)+1))  &
                + A(4)*(f(i+4) - f(mod(n_max+i-6,n_max-1)+1))  &
                + A(3)*(f(i+3) - f(mod(n_max+i-5,n_max-1)+1))  &
                + A(2)*(f(i+2) - f(mod(n_max+i-4,n_max-1)+1))  &
                + A(1)*(f(i+1) - f(mod(n_max+i-3,n_max-1)+1))
            enddo

            do i=n_max-7, n_max-1
                d1f(i)= A(7)*(f(mod(i+6,n_max-1)+1) - f(i-7))       &
                + A(6)*(f(mod(i+5,n_max-1)+1) - f(i-6))             &
                + A(5)*(f(mod(i+4,n_max-1)+1) - f(i-5))             &
                + A(4)*(f(mod(i+3,n_max-1)+1) - f(i-4))             &
                + A(3)*(f(mod(i+2,n_max-1)+1) - f(i-3))             &
                + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-2))             &
                + A(1)*(f(mod(i,n_max-1)+1) - f(i-1))
            enddo
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_INF(f, d1f, h, 7)
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_SUP(f, d1f, n_max, h, n_max-7, shifted)
        end if

        return

    end subroutine

    subroutine D1_ExpCtr_O0Fp6(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d1f(:)

        integer :: i
        real(kind=8) A(6)

        A(1) =  0.91595792650492d0    /h
        A(2) =  -0.34922251636223d0   /h
        A(3) =  0.14398145036906d0    /h
        A(4) =  -0.051236991729043d0  /h
        A(5) =  0.013273181125903d0   /h
        A(6) = -0.0018126562894445d0  /h



        do i=7, n_max-7
            d1f(i)=   A(6)*(f(i+6) - f(i-6))    &
            + A(5)*(f(i+5) - f(i-5))    &
            + A(4)*(f(i+4) - f(i-4))    &
            + A(3)*(f(i+3) - f(i-3))    &
            + A(2)*(f(i+2) - f(i-2))    &
            + A(1)*(f(i+1) - f(i-1))
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 6
                d1f(i)= A(6)*(f(i+6) - f(mod(n_max+i-8,n_max-1)+1))     &
                + A(5)*(f(i+5) - f(mod(n_max+i-7,n_max-1)+1))   &
                + A(4)*(f(i+4) - f(mod(n_max+i-6,n_max-1)+1))   &
                + A(3)*(f(i+3) - f(mod(n_max+i-5,n_max-1)+1))   &
                + A(2)*(f(i+2) - f(mod(n_max+i-4,n_max-1)+1))   &
                + A(1)*(f(i+1) - f(mod(n_max+i-3,n_max-1)+1))
            enddo

            do i=n_max-6, n_max-1
                d1f(i)=   A(6)*(f(mod(i+5,n_max-1)+1) - f(i-6))     &
                + A(5)*(f(mod(i+4,n_max-1)+1) - f(i-5))     &
                + A(4)*(f(mod(i+3,n_max-1)+1) - f(i-4))     &
                + A(3)*(f(mod(i+2,n_max-1)+1) - f(i-3))     &
                + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-2))     &
                + A(1)*(f(mod(i,n_max-1)+1) - f(i-1))
            enddo

        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_INF(f, d1f, h, 6)
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_SUP(f, d1f, n_max, h, n_max-6, shifted)
        end if

        return

    end subroutine

    subroutine D1_ExpCtr_O0Fp6_ACC(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(inout)   :: d1f(:)

        integer :: i
        real(kind=8) A(6)

        A(1) =  0.91595792650492d0    /h
        A(2) =  -0.34922251636223d0   /h
        A(3) =  0.14398145036906d0    /h
        A(4) =  -0.051236991729043d0  /h
        A(5) =  0.013273181125903d0   /h
        A(6) = -0.0018126562894445d0  /h



        do i=7, n_max-7
            d1f(i)=   d1f(i) + A(6)*(f(i+6) - f(i-6))    &
            + A(5)*(f(i+5) - f(i-5))    &
            + A(4)*(f(i+4) - f(i-4))    &
            + A(3)*(f(i+3) - f(i-3))    &
            + A(2)*(f(i+2) - f(i-2))    &
            + A(1)*(f(i+1) - f(i-1))
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 6
                d1f(i)= d1f(i) + A(6)*(f(i+6) - f(mod(n_max+i-8,n_max-1)+1))     &
                + A(5)*(f(i+5) - f(mod(n_max+i-7,n_max-1)+1))   &
                + A(4)*(f(i+4) - f(mod(n_max+i-6,n_max-1)+1))   &
                + A(3)*(f(i+3) - f(mod(n_max+i-5,n_max-1)+1))   &
                + A(2)*(f(i+2) - f(mod(n_max+i-4,n_max-1)+1))   &
                + A(1)*(f(i+1) - f(mod(n_max+i-3,n_max-1)+1))
            enddo

            do i=n_max-6, n_max-1
                d1f(i)=   d1f(i) + A(6)*(f(mod(i+5,n_max-1)+1) - f(i-6))     &
                + A(5)*(f(mod(i+4,n_max-1)+1) - f(i-5))     &
                + A(4)*(f(mod(i+3,n_max-1)+1) - f(i-4))     &
                + A(3)*(f(mod(i+2,n_max-1)+1) - f(i-3))     &
                + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-2))     &
                + A(1)*(f(mod(i,n_max-1)+1) - f(i-1))
            enddo

        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_INF_ACC(f, d1f, h, 6)
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_SUP_ACC(f, d1f, n_max, h, n_max-6, shifted)
        end if

        return

    end subroutine

    subroutine D1_ExpCtr_O0Fp6_MULT(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:), geom_coefs(:)
        real(kind=8), intent(out)   :: d1f(:)

        integer :: i
        real(kind=8) A(6)

        A(1) =  0.91595792650492d0    /h
        A(2) =  -0.34922251636223d0   /h
        A(3) =  0.14398145036906d0    /h
        A(4) =  -0.051236991729043d0  /h
        A(5) =  0.013273181125903d0   /h
        A(6) = -0.0018126562894445d0  /h



        do i=7, n_max-7
            d1f(i)=   ( A(6)*(f(i+6) - f(i-6))    &
            + A(5)*(f(i+5) - f(i-5))    &
            + A(4)*(f(i+4) - f(i-4))    &
            + A(3)*(f(i+3) - f(i-3))    &
            + A(2)*(f(i+2) - f(i-2))    &
            + A(1)*(f(i+1) - f(i-1)) ) * geom_coefs(i)
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 6
                d1f(i)= ( A(6)*(f(i+6) - f(mod(n_max+i-8,n_max-1)+1))     &
                + A(5)*(f(i+5) - f(mod(n_max+i-7,n_max-1)+1))   &
                + A(4)*(f(i+4) - f(mod(n_max+i-6,n_max-1)+1))   &
                + A(3)*(f(i+3) - f(mod(n_max+i-5,n_max-1)+1))   &
                + A(2)*(f(i+2) - f(mod(n_max+i-4,n_max-1)+1))   &
                + A(1)*(f(i+1) - f(mod(n_max+i-3,n_max-1)+1)) ) * geom_coefs(i)
            enddo

            do i=n_max-6, n_max-1
                d1f(i)=   ( A(6)*(f(mod(i+5,n_max-1)+1) - f(i-6))     &
                + A(5)*(f(mod(i+4,n_max-1)+1) - f(i-5))     &
                + A(4)*(f(mod(i+3,n_max-1)+1) - f(i-4))     &
                + A(3)*(f(mod(i+2,n_max-1)+1) - f(i-3))     &
                + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-2))     &
                + A(1)*(f(mod(i,n_max-1)+1) - f(i-1)) ) * geom_coefs(i)
            enddo

        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_INF_MULT(f, d1f, h, 6, geom_coefs)
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_SUP_MULT(f, d1f, n_max, h, n_max-6, shifted, geom_coefs)
        end if

        return

    end subroutine

    subroutine D1_ExpCtr_O0Fp6_MULT_ACC(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:), geom_coefs(:)
        real(kind=8), intent(inout)   :: d1f(:)

        integer :: i
        real(kind=8) A(6)

        A(1) =  0.91595792650492d0    /h
        A(2) =  -0.34922251636223d0   /h
        A(3) =  0.14398145036906d0    /h
        A(4) =  -0.051236991729043d0  /h
        A(5) =  0.013273181125903d0   /h
        A(6) = -0.0018126562894445d0  /h



        do i=7, n_max-7
            d1f(i)=   d1f(i) + ( A(6)*(f(i+6) - f(i-6))    &
            + A(5)*(f(i+5) - f(i-5))    &
            + A(4)*(f(i+4) - f(i-4))    &
            + A(3)*(f(i+3) - f(i-3))    &
            + A(2)*(f(i+2) - f(i-2))    &
            + A(1)*(f(i+1) - f(i-1)) ) * geom_coefs(i)
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 6
                d1f(i)= d1f(i) + ( A(6)*(f(i+6) - f(mod(n_max+i-8,n_max-1)+1))     &
                + A(5)*(f(i+5) - f(mod(n_max+i-7,n_max-1)+1))   &
                + A(4)*(f(i+4) - f(mod(n_max+i-6,n_max-1)+1))   &
                + A(3)*(f(i+3) - f(mod(n_max+i-5,n_max-1)+1))   &
                + A(2)*(f(i+2) - f(mod(n_max+i-4,n_max-1)+1))   &
                + A(1)*(f(i+1) - f(mod(n_max+i-3,n_max-1)+1)) ) * geom_coefs(i)
            enddo

            do i=n_max-6, n_max-1
                d1f(i)=   d1f(i) + ( A(6)*(f(mod(i+5,n_max-1)+1) - f(i-6))     &
                + A(5)*(f(mod(i+4,n_max-1)+1) - f(i-5))     &
                + A(4)*(f(mod(i+3,n_max-1)+1) - f(i-4))     &
                + A(3)*(f(mod(i+2,n_max-1)+1) - f(i-3))     &
                + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-2))     &
                + A(1)*(f(mod(i,n_max-1)+1) - f(i-1)) ) * geom_coefs(i)
            enddo

        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_INF_MULT_ACC(f, d1f, h, 6, geom_coefs)
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_SUP_MULT_ACC(f, d1f, n_max, h, n_max-6, shifted, geom_coefs)
        end if

        return

    end subroutine

    subroutine D1_ExpCtr_O0Fp7(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d1f(:)

        integer :: i
        real(kind=8) A(7)

        A(1) =  0.95396219562045d0    /h
        A(2) =  -0.41234590494721d0   /h
        A(3) =  0.21233981563217d0    /h
        A(4) =  -0.10672135533957d0   /h
        A(5) =  0.046775295199319d0   /h
        A(6) = -0.015323662211638d0   /h
        A(7) = 0.0026664396358809d0   /h

        do i=8,n_max-8
            d1f(i)= A(7)*(f(i+7) - f(i-7))                  &
            + A(6)*(f(i+6) - f(i-6))     &
            + A(5)*(f(i+5) - f(i-5))     &
            + A(4)*(f(i+4) - f(i-4))     &
            + A(3)*(f(i+3) - f(i-3))     &
            + A(2)*(f(i+2) - f(i-2))     &
            + A(1)*(f(i+1) - f(i-1))
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 7
                d1f(i)=   A(7)*(f(i+7) - f(mod(n_max+i-9,n_max-1)+1))             &
                + A(6)*(f(i+6) - f(mod(n_max+i-8,n_max-1)+1))  &
                + A(5)*(f(i+5) - f(mod(n_max+i-7,n_max-1)+1))  &
                + A(4)*(f(i+4) - f(mod(n_max+i-6,n_max-1)+1))  &
                + A(3)*(f(i+3) - f(mod(n_max+i-5,n_max-1)+1))  &
                + A(2)*(f(i+2) - f(mod(n_max+i-4,n_max-1)+1))  &
                + A(1)*(f(i+1) - f(mod(n_max+i-3,n_max-1)+1))
            enddo

            do i=n_max-7, n_max-1
                d1f(i)= A(7)*(f(mod(i+6,n_max-1)+1) - f(i-7))       &
                + A(6)*(f(mod(i+5,n_max-1)+1) - f(i-6))             &
                + A(5)*(f(mod(i+4,n_max-1)+1) - f(i-5))             &
                + A(4)*(f(mod(i+3,n_max-1)+1) - f(i-4))             &
                + A(3)*(f(mod(i+2,n_max-1)+1) - f(i-3))             &
                + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-2))             &
                + A(1)*(f(mod(i,n_max-1)+1) - f(i-1))
            enddo
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_INF(f, d1f, h, 7)
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_SUP(f, d1f, n_max, h, n_max-7, shifted)
        end if

        return

    end subroutine

    subroutine D1_ExpCtr_O0Fp7_ACC(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(inout)   :: d1f(:)

        integer :: i
        real(kind=8) A(7)



        A(1) =  0.95396219562045d0    /h
        A(2) =  -0.41234590494721d0   /h
        A(3) =  0.21233981563217d0    /h
        A(4) =  -0.10672135533957d0   /h
        A(5) =  0.046775295199319d0   /h
        A(6) = -0.015323662211638d0   /h
        A(7) = 0.0026664396358809d0   /h

        do i=8,n_max-8
            d1f(i)= d1f(i) + A(7)*(f(i+7) - f(i-7))                  &
            + A(6)*(f(i+6) - f(i-6))     &
            + A(5)*(f(i+5) - f(i-5))     &
            + A(4)*(f(i+4) - f(i-4))     &
            + A(3)*(f(i+3) - f(i-3))     &
            + A(2)*(f(i+2) - f(i-2))     &
            + A(1)*(f(i+1) - f(i-1))
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 7
                d1f(i)=   d1f(i) + A(7)*(f(i+7) - f(mod(n_max+i-9,n_max-1)+1))             &
                + A(6)*(f(i+6) - f(mod(n_max+i-8,n_max-1)+1))  &
                + A(5)*(f(i+5) - f(mod(n_max+i-7,n_max-1)+1))  &
                + A(4)*(f(i+4) - f(mod(n_max+i-6,n_max-1)+1))  &
                + A(3)*(f(i+3) - f(mod(n_max+i-5,n_max-1)+1))  &
                + A(2)*(f(i+2) - f(mod(n_max+i-4,n_max-1)+1))  &
                + A(1)*(f(i+1) - f(mod(n_max+i-3,n_max-1)+1))
            enddo

            do i=n_max-7, n_max-1
                d1f(i)= d1f(i) + A(7)*(f(mod(i+6,n_max-1)+1) - f(i-7))       &
                + A(6)*(f(mod(i+5,n_max-1)+1) - f(i-6))             &
                + A(5)*(f(mod(i+4,n_max-1)+1) - f(i-5))             &
                + A(4)*(f(mod(i+3,n_max-1)+1) - f(i-4))             &
                + A(3)*(f(mod(i+2,n_max-1)+1) - f(i-3))             &
                + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-2))             &
                + A(1)*(f(mod(i,n_max-1)+1) - f(i-1))
            enddo
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_INF_ACC(f, d1f, h, 7)
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_SUP_ACC(f, d1f, n_max, h, n_max-7, shifted)
        end if

        return

    end subroutine

    subroutine D1_ExpCtr_O0Fp7_MULT(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:), geom_coefs(:)
        real(kind=8), intent(out)     :: d1f(:)

        integer :: i
        real(kind=8) A(7)

        A(1) =  0.95396219562045d0    /h
        A(2) =  -0.41234590494721d0   /h
        A(3) =  0.21233981563217d0    /h
        A(4) =  -0.10672135533957d0   /h
        A(5) =  0.046775295199319d0   /h
        A(6) = -0.015323662211638d0   /h
        A(7) = 0.0026664396358809d0   /h

        do i=8,n_max-8
            d1f(i)= ( A(7)*(f(i+7) - f(i-7))                  &
            + A(6)*(f(i+6) - f(i-6))     &
            + A(5)*(f(i+5) - f(i-5))     &
            + A(4)*(f(i+4) - f(i-4))     &
            + A(3)*(f(i+3) - f(i-3))     &
            + A(2)*(f(i+2) - f(i-2))     &
            + A(1)*(f(i+1) - f(i-1)) ) * geom_coefs(i)
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 7
                d1f(i)=   ( A(7)*(f(i+7) - f(mod(n_max+i-9,n_max-1)+1))             &
                + A(6)*(f(i+6) - f(mod(n_max+i-8,n_max-1)+1))  &
                + A(5)*(f(i+5) - f(mod(n_max+i-7,n_max-1)+1))  &
                + A(4)*(f(i+4) - f(mod(n_max+i-6,n_max-1)+1))  &
                + A(3)*(f(i+3) - f(mod(n_max+i-5,n_max-1)+1))  &
                + A(2)*(f(i+2) - f(mod(n_max+i-4,n_max-1)+1))  &
                + A(1)*(f(i+1) - f(mod(n_max+i-3,n_max-1)+1)) ) * geom_coefs(i)
            enddo

            do i=n_max-7, n_max-1
                d1f(i)= ( A(7)*(f(mod(i+6,n_max-1)+1) - f(i-7))       &
                + A(6)*(f(mod(i+5,n_max-1)+1) - f(i-6))             &
                + A(5)*(f(mod(i+4,n_max-1)+1) - f(i-5))             &
                + A(4)*(f(mod(i+3,n_max-1)+1) - f(i-4))             &
                + A(3)*(f(mod(i+2,n_max-1)+1) - f(i-3))             &
                + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-2))             &
                + A(1)*(f(mod(i,n_max-1)+1) - f(i-1)) ) * geom_coefs(i)
            enddo
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_INF_MULT(f, d1f, h, 7, geom_coefs)
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_SUP_MULT(f, d1f, n_max, h, n_max-7, shifted, geom_coefs)
        end if

        return

    end subroutine

    subroutine D1_ExpCtr_O0Fp7_MULT_ACC(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:), geom_coefs(:)
        real(kind=8), intent(inout)   :: d1f(:)

        integer :: i
        real(kind=8) A(7)

        A(1) =  0.95396219562045d0    /h
        A(2) =  -0.41234590494721d0   /h
        A(3) =  0.21233981563217d0    /h
        A(4) =  -0.10672135533957d0   /h
        A(5) =  0.046775295199319d0   /h
        A(6) = -0.015323662211638d0   /h
        A(7) = 0.0026664396358809d0   /h

        do i=8,n_max-8
            d1f(i)= d1f(i) + ( A(7)*(f(i+7) - f(i-7))                  &
            + A(6)*(f(i+6) - f(i-6))     &
            + A(5)*(f(i+5) - f(i-5))     &
            + A(4)*(f(i+4) - f(i-4))     &
            + A(3)*(f(i+3) - f(i-3))     &
            + A(2)*(f(i+2) - f(i-2))     &
            + A(1)*(f(i+1) - f(i-1)) ) * geom_coefs(i)
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 7
                d1f(i)=   d1f(i) + ( A(7)*(f(i+7) - f(mod(n_max+i-9,n_max-1)+1))             &
                + A(6)*(f(i+6) - f(mod(n_max+i-8,n_max-1)+1))  &
                + A(5)*(f(i+5) - f(mod(n_max+i-7,n_max-1)+1))  &
                + A(4)*(f(i+4) - f(mod(n_max+i-6,n_max-1)+1))  &
                + A(3)*(f(i+3) - f(mod(n_max+i-5,n_max-1)+1))  &
                + A(2)*(f(i+2) - f(mod(n_max+i-4,n_max-1)+1))  &
                + A(1)*(f(i+1) - f(mod(n_max+i-3,n_max-1)+1)) ) * geom_coefs(i)
            enddo

            do i=n_max-7, n_max-1
                d1f(i)= d1f(i) + ( A(7)*(f(mod(i+6,n_max-1)+1) - f(i-7))       &
                + A(6)*(f(mod(i+5,n_max-1)+1) - f(i-6))             &
                + A(5)*(f(mod(i+4,n_max-1)+1) - f(i-5))             &
                + A(4)*(f(mod(i+3,n_max-1)+1) - f(i-4))             &
                + A(3)*(f(mod(i+2,n_max-1)+1) - f(i-3))             &
                + A(2)*(f(mod(i+1,n_max-1)+1) - f(i-2))             &
                + A(1)*(f(mod(i,n_max-1)+1) - f(i-1)) ) * geom_coefs(i)
            enddo
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_INF_MULT_ACC(f, d1f, h, 7, geom_coefs)
        end if

        if (bc_type.eq.Dirichlet) then
            call D1c_CTR_CLOSURE_SUP_MULT_ACC(f, d1f, n_max, h, n_max-7, shifted, geom_coefs)
        end if

        return

    end subroutine

    subroutine D1_CptCtr_O6Fp0(f, d1f, n_max, h, shifted, bc_type)

        use Thomas_solver
        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d1f(:)

        integer :: i
        real(kind=8), dimension(:,:), save, allocatable :: A

        real(kind=8) rhs(n_max)

        real(kind=8), save    :: B(0:2),  B_lim(-1:2)

        integer, save :: last_nmax=0, last_bc_type1, last_bc_type2
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h



        update_matrix=(last_nmax/=n_max).or.(last_bc_type1/=bc_type).or.(last_bc_type2/=bc_type)    &
        .or.(last_shifted.neqv.shifted).or.(last_h/=h)

        ! MATRIX A, B definition -----------------------------------------------------------
        !if (update_matrix) then
        call fill_matrix
        !end if

        last_nmax=n_max
        last_bc_type1=bc_type
        last_bc_type2=bc_type
        last_shifted=shifted
        last_h=h

        ! SCHEME CALCULATION ---------------------------------------------------------------
        do i=3, n_max-3
            rhs(i)= B(1)*( f(i+1) - f(i-1) ) + B(2)*( f(i+2) - f(i-2) )
        enddo

        if (bc_type.eq.periodic) then

            rhs(1)= B(1)*( f(2) - f(n_max-1) ) + B(2)*( f(3) - f(n_max-2) )
            rhs(2)= B(1)*( f(3) - f(1) ) + B(2)*( f(4) - f(n_max-1) )

            rhs(n_max-1)= B(1)*( f(1) - f(n_max-2) ) + B(2)*( f(2) - f(n_max-3) )
            rhs(n_max-2)= B(1)*( f(n_max-1) - f(n_max-3) ) + B(2)*( f(1) - f(n_max-4) )

            call TS_Pr(A(-1,1:n_max-1), A(0,1:n_max-1), A(1,1:n_max-1), rhs(1:n_max-1), d1f(1:n_max-1), n_max-1)

        end if

        if ((bc_type.eq.Dirichlet)) then

            rhs(2)= B_lim(-1)*f(1) + B_lim(0) * f(2) + B_lim(1) * f(3) + B_lim(2) * f(4)

            if (shifted) then

                rhs(n_max-2)= -( B_lim(-1)*f(n_max-1) + B_lim(0)*f(n_max-2) + B_lim(1)*f(n_max-3)+ B_lim(2)*f(n_max-4) )

                call TS_NPr(rhs(2:n_max-2), d1f(2:n_max-2), n_max-3)

            else

                rhs(n_max-2)= B(1)*( f(n_max-1) - f(n_max-3) ) + B(2)*( f(n_max) - f(n_max-4) )

                rhs(n_max-1)= -(B_lim(-1)*f(n_max) + B_lim(0)*f(n_max-1) + B_lim(1)*f(n_max-2)+ B_lim(2)*f(n_max-3))

                call TS_NPr(rhs(2:n_max-1), d1f(2:n_max-1), n_max-2)

            end if

        end if

        return

    contains


        subroutine fill_matrix()

            if (allocated(A)) then
                deallocate(A)
            end if

            allocate(A(-1:1, n_max))

            B(1)=14.d0      / (18.d0*h)
            B(2)= 1.d0      / (36.d0*h)

            do i=3, n_max-3
                A(-1,i) = 1.d0/3.d0
                A(0,i) = 1.d0
                A(1,i) = 1.d0/3.d0
            enddo

            select case (bc_type)

                case (periodic)

                    A(-1,1:2) = 1.d0/3.d0
                    A(0,1:2) = 1.d0
                    A(1,1:2) = 1.d0/3.d0

                    A(-1,n_max-2:n_max-1) = 1.d0/3.d0
                    A(0,n_max-2: n_max-1) = 1.d0
                    A(1,n_max-2:n_max-1) = 1.d0/3.d0

                case (Dirichlet)

                    A(-1,2)=0.d0
                    A(0,2)=1.d0
                    A(1,2)=1.d0

                    B_lim(-1)   =   -1.d0/6.d0              / h
                    B_lim(0)    =   -3.d0/2.d0              / h
                    B_lim(1)    =   3.d0/2.d0               / h
                    B_lim(2)    =   1.d0/6.d0               / h

                    if (shifted) then

                        A(-1,n_max-2)=1.d0
                        A(0,n_max-2)=1.d0
                        A(1,n_max-2)=0.d0

                        call TS_init(A(-1,2:n_max-2), A(0,2:n_max-2), A(1,2:n_max-2), n_max-3)

                    else

                        A(-1,n_max-2) = 1.d0/3.d0
                        A(0,n_max-2) = 1.d0
                        A(1,n_max-2) = 1.d0/3.d0

                        A(-1,n_max-1)=1.d0
                        A(0,n_max-1)=1.d0
                        A(1,n_max-1)=0.d0

                        call TS_init(A(-1,2:n_max-1), A(0,2:n_max-1), A(1,2:n_max-1), n_max-2)

                    end if

                case default

            end select

        end subroutine fill_matrix


    end subroutine

    subroutine D1_CptCtr_O6Fp0_ACC(f, d1f, n_max, h, shifted, bc_type)

        implicit none
        logical, intent(in)             :: shifted
        integer, intent(in)             :: n_max, bc_type
        real(kind=8) , intent(in)             :: h,  f(:)
        real(kind=8), intent(inout)         :: d1f(:)

        real(kind=8), dimension(n_max) :: d1f_tmp
        integer :: i

        select case (bc_type)

            case (periodic)

                d1f_tmp(1:n_max-1)=d1f(1:n_max-1)
                call D1_CptCtr_O6Fp0(f, d1f_tmp, n_max, h, shifted, bc_type)
                d1f(1:n_max-1) = d1f(1:n_max-1) + d1f_tmp(1:n_max-1)

            case (Dirichlet)

                if (shifted) then

                    d1f_tmp(2:n_max-2)=d1f(2:n_max-2)
                    call D1_CptCtr_O6Fp0(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(2:n_max-2) = d1f(2:n_max-2) + d1f_tmp(2:n_max-2)

                else

                    d1f_tmp(2:n_max-1)=d1f(2:n_max-1)
                    call D1_CptCtr_O6Fp0(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(2:n_max-1) = d1f(2:n_max-1) + d1f_tmp(2:n_max-1)
                end if

            case default

        end select

        return

    end subroutine

    subroutine D1_CptCtr_O6Fp0_MULT(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        logical, intent(in)             :: shifted
        integer, intent(in)             :: n_max, bc_type
        real(kind=8) , intent(in)             :: h,  geom_coefs(:), f(:)
        real(kind=8), intent(out)             :: d1f(:)

        integer :: i

        call D1_CptCtr_O6Fp0(f, d1f, n_max, h, shifted, bc_type)

        select case (bc_type)

            case (periodic)

                do i = 1, n_max-1
                    d1f(i)=d1f(i)*geom_coefs(i)
                end do

            case (Dirichlet)
                if (shifted) then
                    d1f(2:n_max-2)=d1f(2:n_max-2)*geom_coefs(2:n_max-2)
                else
                    d1f(2:n_max-1)=d1f(2:n_max-1)*geom_coefs(2:n_max-1)
                end if

            case default

        end select

        return

    end subroutine

    subroutine D1_CptCtr_O6Fp0_MULT_ACC(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        logical, intent(in)             :: shifted
        integer, intent(in)             :: n_max, bc_type
        real(kind=8) , intent(in)             :: h,  geom_coefs(:), f(:)
        real(kind=8), intent(inout)             :: d1f(:)

        real(kind=8), dimension(n_max) :: d1f_tmp
        integer :: i

        select case (bc_type)

            case (periodic)

                d1f_tmp(1:n_max-1)=d1f(1:n_max-1)
                call D1_CptCtr_O6Fp0(f, d1f_tmp, n_max, h, shifted, bc_type)
                d1f(1:n_max-1) = d1f(1:n_max-1) + d1f_tmp(1:n_max-1)*geom_coefs(1:n_max-1)

            case (Dirichlet)

                if (shifted) then

                    d1f_tmp(2:n_max-2)=d1f(2:n_max-2)
                    call D1_CptCtr_O6Fp0(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(2:n_max-2) = d1f(2:n_max-2) + d1f_tmp(2:n_max-2)*geom_coefs(2:n_max-2)

                else

                    d1f_tmp(2:n_max-1)=d1f(2:n_max-1)
                    call D1_CptCtr_O6Fp0(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(2:n_max-1) = d1f(2:n_max-1) + d1f_tmp(2:n_max-1)*geom_coefs(2:n_max-1)
                end if

            case default

        end select

        return

    end subroutine

    subroutine D1c_KIM1996_CPT_O6_opt(f, d1f, n_max, h, shifted, bc_type)

        use Thomas_solver
        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d1f(:)

        integer :: i
        real(kind=8), dimension(:,:), save, allocatable :: A, B
        real(kind=8), dimension(:), save, allocatable   :: u, v
        real(kind=8), save                              :: vq
        real(kind=8) rhs(n_max)
        real(kind=8), save    :: Bs(-1:3), Be(-1:3)

        integer, save :: last_nmax=0, last_bc_type1, last_bc_type2
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h



        update_matrix=(last_nmax/=n_max).or.(last_bc_type1/=bc_type).or.(last_bc_type2/=bc_type)    &
        .or.(last_shifted.neqv.shifted).or.(last_h/=h)


        !write(*,*)'D1_OUCS3'

        ! MATRIX A, B definition -----------------------------------------------------------


        ! MATRIX A, B definition -----------------------------------------------------------
        if (update_matrix) then
            write(*,*)'updating'

            call fill_matrix
        end if

        if (n_max==0) then

            write(*,*)'update_matrix', update_matrix
            write(*,*)'nmax', last_nmax, n_max
            write(*,*)'bc_type', last_bc_type1, bc_type
            write(*,*)'bc_type2', last_bc_type2, bc_type
            write(*,*)'shifted', last_shifted, shifted
            write(*,*)'h', last_h, h
        end if


        last_nmax=n_max
        last_bc_type1=bc_type
        last_bc_type2=bc_type
        last_shifted=shifted
        last_h=h

        ! SCHEME CALCULATION ---------------------------------------------------------------
        do i=4, n_max-4
            rhs(i)=     B(i, 2)*( f(i+2) - f(i-2))     &
            +   B(i, 1)*( f(i+1) - f(i-1) )            &
            +   B(i, 3)*( f(i+3) - f(i-3) )
        enddo


        if (bc_type.eq.periodic) then

            rhs(1)=     B(1, 1)*( f(2) - f(n_max-1) )  &
            +   B(1, 2)*( f(3) - f(n_max-2) )          &
            +   B(1, 3)*( f(4) - f(n_max-3) )

            rhs(2)=     B(2, 1)*( f(3) - f(1) )  &
            +   B(2, 2)*( f(4) - f(n_max-1) )      &
            +   B(2, 3)*( f(5) - f(n_max-2) )

            rhs(3)=     B(3, 1)*( f(4) - f(2) )  &
            +   B(3, 2)*( f(5) - f(1) )      &
            +   B(3, 3)*( f(6) - f(n_max-1) )


            rhs(n_max-1)=   B(n_max-1, 1)*( f(1) - f(n_max-2) )  &
            +   B(n_max-1, 2)*( f(2) - f(n_max-3) )              &
            +   B(n_max-1, 3)*( f(3) - f(n_max-4) )


            rhs(n_max-2)=   B(n_max-2, 1)*( f(n_max-1) - f(n_max-3) )  &
            +   B(n_max-2, 2)*( f(1) - f(n_max-4) )                      &
            +   B(n_max-2, 3)*( f(2) - f(n_max-5) )


            rhs(n_max-3)=   B(n_max-3, 1)*( f(n_max-2) - f(n_max-4) )  &
            +   B(n_max-3, 2)*( f(n_max-1) - f(n_max-5) )                      &
            +   B(n_max-3, 3)*( f(1) - f(n_max-6) )

            call TS_Pr_opt(A(-1,1:n_max-1), A(0,1:n_max-1), A(1,1:n_max-1), rhs(1:n_max-1), d1f(1:n_max-1), u,v,vq,n_max-1)

        end if


        if ((bc_type.eq.Dirichlet)) then

            rhs(2)= Bs(-1)*f(1) + Bs(0) * f(2)          &
            + Bs(1) * f(3) + Bs(2) * f(4)

            rhs(3)=     B(3, 1)*( f(4) - f(2) )  &
            +   B(3, 2)*( f(5) - f(1) )

            if (shifted) then

                rhs(n_max-2)= -( + Be(-1)*f(n_max-1) + Be(0)*f(n_max-2)       &
                + Be(1)*f(n_max-3)+ Be(2)*f(n_max-4))


                rhs(n_max-3)=   B(n_max-3, 1)*( f(n_max-2) - f(n_max-4) )  &
                +   B(n_max-3, 2)*( f(n_max-1) - f(n_max-5) )

                call TS_NPr(rhs(2:n_max-2), d1f(2:n_max-2), n_max-3)

            else


                rhs(n_max-3)=   B(n_max-3, 1)*( f(n_max-2) - f(n_max-4) )  &
                +   B(n_max-3, 2)*( f(n_max-1) - f(n_max-5) )                      &
                +   B(n_max-3, 3)*( f(1) - f(n_max-6) )

                rhs(n_max-2)= B(n_max-2, 1)*( f(n_max-1) - f(n_max-3) )  &
                + B(n_max-2, 2)*( f(n_max)   - f(n_max-4) )

                rhs(n_max-1)= -( + Be(-1)*f(n_max) + Be(0)*f(n_max-1)     &
                + Be(1)*f(n_max-2)+ Be(2)*f(n_max-3))

                call TS_NPr(rhs(2:n_max-1), d1f(2:n_max-1), n_max-2)

            end if

        end if

        return


    contains

        subroutine fill_matrix()

            if (allocated(A)) then
                deallocate(A)
                deallocate(B)
                deallocate(u)
                deallocate(v)
            end if

            allocate(A(-1:1, n_max))
            allocate(u(n_max), v(n_max))
            allocate(B(n_max, 3))

            do i=4, n_max-4
                A(-1,i) = 0.408589269d0
                A(0,i) = 1.d0
                A(1,i) = 0.408589269d0

                B(i, 1)=1.568098212d0 / (2.d0*h)
                B(i, 2)= 0.271657107d0 / (4.d0*h)
                B(i, 3)= -0.022576781d0 / (6.d0*h)
            enddo

            select case (bc_type)

                case (periodic)

                    A(-1,1:3) = 0.408589269
                    A(0,1:3) = 1.d0
                    A(1,1:3) = 0.408589269

                    A(-1,n_max-3:n_max-1) = 0.408589269
                    A(0,n_max-3: n_max-1) = 1.d0
                    A(1,n_max-3:n_max-1) = 0.408589269

                    B(1:3, 1)=1.568098212 / (2.d0*h)
                    B(1:3, 2)= 0.271657107 / (4.d0*h)
                    B(1:3, 3)= -0.022576781 / (6.d0*h)

                    B(n_max-3:n_max, 1)=1.568098212 / (2.d0*h)
                    B(n_max-3:n_max, 2)= 0.271657107 / (4.d0*h)
                    B(n_max-3:n_max, 3)= -0.022576781 / (6.d0*h)

                    call TS_init_Pr(A(-1,1:n_max-1), A(0,1:n_max-1), A(1,1:n_max-1), u, v, vq,n_max-1)

                case (Dirichlet)

                    A(-1,2)=0.d0
                    A(0,2)=1.d0
                    A(1,2)=1.d0

                    A(-1,3)=0.33333333d0
                    A(0,3)=1.d0
                    A(1,3)=0.33333333d0

                    Bs(-1)=-1.d0/6.d0           / h
                    Bs(0)=-3.d0/2.d0            / h
                    Bs(1)=3.d0/2.d0             / h
                    Bs(2)=1.d0/6.d0             / h

                    B(3, 1)=14.0/ (18.d0*h)
                    B(3, 2)= 1.d0 / (36.d0*h)
                    B(3, 3)= 0.d0

                    ! End coefficients -----------------------

                    Be(-1)=-1.d0/6.d0           / h
                    Be(0)=-3.d0/2.d0            / h
                    Be(1)=3.d0/2.d0             / h
                    Be(2)=1.d0/6.d0             / h

                    if (shifted) then

                        A(-1,n_max-3) = 0.33333333d0
                        A(0,n_max-3) = 1.d0
                        A(1,n_max-3) = 0.33333333d0

                        A(-1,n_max-2)=1.d0
                        A(0,n_max-2)=1.d0
                        A(1,n_max-2)=0.d0

                        B(n_max-3, 1)=14.0/ (18.d0*h)
                        B(n_max-3, 2)= 1.d0 / (36.d0*h)
                        B(n_max-3, 3)= 0.d0


                        call TS_init(A(-1,2:n_max-2), A(0,2:n_max-2), A(1,2:n_max-2), n_max-3)

                    else

                        A(-1,n_max-3) = 0.408589269d0
                        A(0,n_max-3) = 1.d0
                        A(1,n_max-3) = 0.408589269d0

                        A(-1,n_max-2) = 0.33333333d0
                        A(0,n_max-2) = 1.d0
                        A(1,n_max-2) = 0.33333333d0

                        A(-1,n_max-1)=1.d0
                        A(0,n_max-1)=1.d0
                        A(1,n_max-1)=0.d0

                        B(n_max-3, 1)=1.568098212 / (2.d0*h)
                        B(n_max-3, 2)= 0.271657107 / (4.d0*h)
                        B(n_max-3, 3)= -0.022576781 / (6.d0*h)

                        B(n_max-2, 1)=14.0/ (18.d0*h)
                        B(n_max-2, 2)= 1.d0 / (36.d0*h)
                        B(n_max-2, 3)= 0.d0


                        call TS_init(A(-1,2:n_max-1), A(0,2:n_max-1), A(1,2:n_max-1), n_max-2)

                    end if


                case default

            end select

        end subroutine fill_matrix


    end subroutine

    subroutine D1c_KIM1996_CPT_O6_opt_ACC(f, d1f, n_max, h, shifted, bc_type)

        implicit none
        logical, intent(in)             :: shifted
        integer, intent(in)             :: n_max, bc_type
        real(kind=8) , intent(in)             :: h,  f(:)
        real(kind=8), intent(inout)         :: d1f(:)

        real(kind=8), dimension(n_max) :: d1f_tmp
        integer :: i

        select case (bc_type)

            case (periodic)

                d1f_tmp(1:n_max-1)=d1f(1:n_max-1)
                call D1c_KIM1996_CPT_O6_opt(f, d1f_tmp, n_max, h, shifted, bc_type)
                d1f(1:n_max-1) = d1f(1:n_max-1) + d1f_tmp(1:n_max-1)

            case (Dirichlet)

                if (shifted) then

                    d1f_tmp(2:n_max-2)=d1f(2:n_max-2)
                    call D1c_KIM1996_CPT_O6_opt(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(2:n_max-2) = d1f(2:n_max-2) + d1f_tmp(2:n_max-2)

                else

                    d1f_tmp(2:n_max-1)=d1f(2:n_max-1)
                    call D1c_KIM1996_CPT_O6_opt(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(2:n_max-1) = d1f(2:n_max-1) + d1f_tmp(2:n_max-1)
                end if

            case default

        end select

        return


    end subroutine

    subroutine D1c_KIM1996_CPT_O6_opt_MULT(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        logical, intent(in)             :: shifted
        integer, intent(in)             :: n_max, bc_type
        real(kind=8) , intent(in)             :: h,  geom_coefs(:), f(:)
        real(kind=8), intent(out)             :: d1f(:)

        integer :: i

        call D1c_KIM1996_CPT_O6_opt(f, d1f, n_max, h, shifted, bc_type)

        select case (bc_type)

            case (periodic)

                do i = 1, n_max-1
                    d1f(i)=d1f(i)*geom_coefs(i)
                end do

            case (Dirichlet)
                if (shifted) then
                    d1f(2:n_max-2)=d1f(2:n_max-2)*geom_coefs(2:n_max-2)
                else
                    d1f(2:n_max-1)=d1f(2:n_max-1)*geom_coefs(2:n_max-1)
                end if

            case default

        end select

        return

    end subroutine

    subroutine D1c_KIM1996_CPT_O6_opt_MULT_ACC(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        logical, intent(in)             :: shifted
        integer, intent(in)             :: n_max, bc_type
        real(kind=8) , intent(in)             :: h,  geom_coefs(:), f(:)
        real(kind=8), intent(inout)             :: d1f(:)

        real(kind=8), dimension(n_max) :: d1f_tmp
        integer :: i

        select case (bc_type)

            case (periodic)

                d1f_tmp(1:n_max-1)=d1f(1:n_max-1)
                call D1c_KIM1996_CPT_O6_opt(f, d1f_tmp, n_max, h, shifted, bc_type)
                d1f(1:n_max-1) = d1f(1:n_max-1) + d1f_tmp(1:n_max-1)*geom_coefs(1:n_max-1)

            case (Dirichlet)

                if (shifted) then

                    d1f_tmp(2:n_max-2)=d1f(2:n_max-2)
                    call D1c_KIM1996_CPT_O6_opt(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(2:n_max-2) = d1f(2:n_max-2) + d1f_tmp(2:n_max-2)*geom_coefs(2:n_max-2)

                else

                    d1f_tmp(2:n_max-1)=d1f(2:n_max-1)
                    call D1c_KIM1996_CPT_O6_opt(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(2:n_max-1) = d1f(2:n_max-1) + d1f_tmp(2:n_max-1)*geom_coefs(2:n_max-1)
                end if

            case default

        end select

        return

    end subroutine

    subroutine D1_OUCS3(f, d1f, n_max, h, shifted, bc_type)

        use Thomas_solver
        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d1f(:)

        integer :: i
        real(kind=8), dimension(:,:), save, allocatable :: A
        real(kind=8), dimension(:), save, allocatable :: u, v
        real(kind=8), save                              :: vq
        real(kind=8) rhs(n_max)
        real(kind=8), save    :: B(2), Bs(-1:3), Be(-1:3)

        integer, save :: last_nmax=0, last_bc_type1, last_bc_type2
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h



        update_matrix=(last_nmax/=n_max).or.(last_bc_type1/=bc_type).or.(last_bc_type2/=bc_type)    &
        .or.(last_shifted.neqv.shifted).or.(last_h/=h)


        !write(*,*)'D1_OUCS3'

        ! MATRIX A, B definition -----------------------------------------------------------
        if (update_matrix) then
            call fill_matrix
        end if
        call fill_matrix

        last_nmax=n_max
        last_bc_type1=bc_type
        last_bc_type2=bc_type
        last_shifted=shifted
        last_h=h

        ! SCHEME CALCULATION ---------------------------------------------------------------
        do i=3, n_max-3
            rhs(i)=     B(2)*( f(i+2) - f(i-2))   &
            +   B(1)*( f(i+1) - f(i-1) )
        enddo


        if (bc_type.eq.periodic) then

            rhs(1)=     B(1)*( f(2) - f(n_max-1) )  &
            +   B(2)*( f(3) - f(n_max-2) )

            rhs(2)=     B(1)*( f(3) - f(1) )  &
            +   B(2)*( f(4) - f(n_max-1) )


            rhs(n_max-1)=   B(1)*( f(1) - f(n_max-2) )  &
            +   B(2)*( f(2) - f(n_max-3) )


            rhs(n_max-2)=   B(1)*( f(n_max-1) - f(n_max-3) )  &
            +   B(2)*( f(1) - f(n_max-4) )

            call TS_Pr_opt(A(-1,1:n_max-1), A(0,1:n_max-1), A(1,1:n_max-1), rhs(1:n_max-1), d1f(1:n_max-1), u,v,vq,n_max-1)

        end if


        if ((bc_type.eq.Dirichlet)) then

            rhs(2)= Bs(-1)*f(1) + Bs(0) * f(2)          &
            + Bs(1) * f(3) + Bs(2) * f(4)        &
            + Bs(3) * f(5)

            if (shifted) then

                rhs(n_max-2)= -( + Be(-1)*f(n_max-1) + Be(0)*f(n_max-2)       &
                + Be(1)*f(n_max-3)+ Be(2)*f(n_max-4)  &
                + Be(3)*f(n_max-5) )

                call TS_NPr(rhs(2:n_max-2), d1f(2:n_max-2), n_max-3)

            else

                rhs(n_max-2)= B(1)*( f(n_max-1) - f(n_max-3) )  &
                + B(2)*( f(n_max)   - f(n_max-4) )

                rhs(n_max-1)= -( + Be(-1)*f(n_max) + Be(0)*f(n_max-1)     &
                + Be(1)*f(n_max-2)+ Be(2)*f(n_max-3)      &
                + Be(3)*f(n_max-4) )

                call TS_NPr(rhs(2:n_max-1), d1f(2:n_max-1), n_max-2)

            end if

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

            B(1)=1.57557379d0 / (2.d0*h)
            B(2)= 0.1832051925d0 / (4.d0*h)

            do i=3, n_max-3
                A(-1,i) = 0.3793894912d0
                A(0,i) = 1.d0
                A(1,i) = 0.3793894912d0
            enddo

            select case (bc_type)

                case (periodic)

                    A(-1,1:2) = 0.3793894912d0
                    A(0,1:2) = 1.d0
                    A(1,1:2) = 0.3793894912d0

                    A(-1,n_max-2:n_max-1) = 0.3793894912d0
                    A(0,n_max-2: n_max-1) = 1.d0
                    A(1,n_max-2:n_max-1) = 0.3793894912d0

                    call TS_init_Pr(A(-1,1:n_max-1), A(0,1:n_max-1), A(1,1:n_max-1), u, v, vq, n_max-1)

                case (Dirichlet)

                    A(-1,2)=0.d0
                    A(0,2)=1.d0
                    A(1,2)=0.d0

                    Bs(-1)=-0.35d0             / h
                    Bs(0)=-0.43333333333333d0  / h
                    Bs(1)=0.9d0                / h
                    Bs(2)=-0.1d0               / h
                    Bs(3)=-0.016666666666667d0 / h

                    ! End coefficients -----------------------

                    Be(-1)=-0.35d0             / h
                    Be(0)=-0.43333333333333d0  / h
                    Be(1)=0.9d0                / h
                    Be(2)=-0.1d0               / h
                    Be(3)=-0.016666666666667d0 / h

                    if (shifted) then

                        A(-1,n_max-2)=0.d0
                        A(0,n_max-2)=1.d0
                        A(1,n_max-2)=0.d0

                        call TS_init(A(-1,2:n_max-2), A(0,2:n_max-2), A(1,2:n_max-2), n_max-3)

                    else

                        A(-1,n_max-2) = 0.3793894912d0
                        A(0,n_max-2) = 1.d0
                        A(1,n_max-2) = 0.3793894912d0

                        A(-1,n_max-1)=0.d0
                        A(0,n_max-1)=1.d0
                        A(1,n_max-1)=0.d0

                        call TS_init(A(-1,2:n_max-1), A(0,2:n_max-1), A(1,2:n_max-1), n_max-2)

                    end if


                case default

            end select

        end subroutine fill_matrix


    end subroutine

    subroutine D1_OUCS3_ACC(f, d1f, n_max, h, shifted, bc_type)

        implicit none
        logical, intent(in)             :: shifted
        integer, intent(in)             :: n_max, bc_type
        real(kind=8) , intent(in)             :: h,  f(:)
        real(kind=8), intent(inout)         :: d1f(:)

        real(kind=8), dimension(n_max) :: d1f_tmp
        integer :: i

        select case (bc_type)

            case (periodic)

                d1f_tmp(1:n_max-1)=d1f(1:n_max-1)
                call D1_OUCS3(f, d1f_tmp, n_max, h, shifted, bc_type)
                d1f(1:n_max-1) = d1f(1:n_max-1) + d1f_tmp(1:n_max-1)

            case (Dirichlet)

                if (shifted) then

                    d1f_tmp(2:n_max-2)=d1f(2:n_max-2)
                    call D1_OUCS3(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(2:n_max-2) = d1f(2:n_max-2) + d1f_tmp(2:n_max-2)

                else

                    d1f_tmp(2:n_max-1)=d1f(2:n_max-1)
                    call D1_OUCS3(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(2:n_max-1) = d1f(2:n_max-1) + d1f_tmp(2:n_max-1)
                end if

            case default

        end select

        return

    end subroutine

    subroutine D1_OUCS3_MULT(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        logical, intent(in)             :: shifted
        integer, intent(in)             :: n_max, bc_type
        real(kind=8) , intent(in)             :: h,  geom_coefs(:), f(:)
        real(kind=8), intent(out)             :: d1f(:)

        integer :: i

        call D1_OUCS3(f, d1f, n_max, h, shifted, bc_type)

        select case (bc_type)

            case (periodic)

                do i = 1, n_max-1
                    d1f(i)=d1f(i)*geom_coefs(i)
                end do

            case (Dirichlet)
                if (shifted) then
                    d1f(2:n_max-2)=d1f(2:n_max-2)*geom_coefs(2:n_max-2)
                else
                    d1f(2:n_max-1)=d1f(2:n_max-1)*geom_coefs(2:n_max-1)
                end if

            case default

        end select

        return

    end subroutine

    subroutine D1_OUCS3_MULT_ACC(f, d1f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        logical, intent(in)             :: shifted
        integer, intent(in)             :: n_max, bc_type
        real(kind=8) , intent(in)             :: h,  geom_coefs(:), f(:)
        real(kind=8), intent(inout)             :: d1f(:)

        real(kind=8), dimension(n_max) :: d1f_tmp
        integer :: i

        select case (bc_type)

            case (periodic)

                d1f_tmp(1:n_max-1)=d1f(1:n_max-1)
                call D1_OUCS3(f, d1f_tmp, n_max, h, shifted, bc_type)
                d1f(1:n_max-1) = d1f(1:n_max-1) + d1f_tmp(1:n_max-1)*geom_coefs(1:n_max-1)

            case (Dirichlet)

                if (shifted) then

                    d1f_tmp(2:n_max-2)=d1f(2:n_max-2)
                    call D1_OUCS3(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(2:n_max-2) = d1f(2:n_max-2) + d1f_tmp(2:n_max-2)*geom_coefs(2:n_max-2)

                else

                    d1f_tmp(2:n_max-1)=d1f(2:n_max-1)
                    call D1_OUCS3(f, d1f_tmp, n_max, h, shifted, bc_type)
                    d1f(2:n_max-1) = d1f(2:n_max-1) + d1f_tmp(2:n_max-1)*geom_coefs(2:n_max-1)
                end if

            case default

        end select

        return

    end subroutine

    subroutine D1c_CTR_CLOSURE_INF(f, d1f, h, imax)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(imax)

        implicit none

        integer :: i, n_max, imax
        real(kind=8) h
        real(kind=8) d1f(:)
        real(kind=8) f(:)

        real(kind=8) d1f_temp(CPT_MIN)


        if (treat_boundary_by_cpt_scheme) then

            call D1_OUCS3(  &
            f(1:CPT_MIN), d1f_temp(1:CPT_MIN), CPT_MIN, h, .true.,      &
            Dirichlet)

            d1f(2)=d1f_temp(2)
            if (imax.eq.2) return

            d1f(3)=d1f_temp(3)
            if (imax.eq.3) return

        else

            ! i=2---------------------------------------------------
            i=2
            d1f(i)= (f(i+1)-f(i-1)) /   (2.d0*h)

            if (i.eq.imax) return
            ! i=3---------------------------------------------------
            i=3
            d1f(i)    =   2.d0*(f(i+1)-f(i-1))  / (3.d0*h)      &
            -   (f(i+2) - f(i-2))         / (12.d0*h)

            if (i.eq.imax) return

        endif

        ! i=4---------------------------------------------------
        i=4
        d1f(i)  =   0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h

        if (i.eq.imax) return

        ! i=5---------------------------------------------------
        i=5
        d1f(i)  =   0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h

        if (i.eq.imax) return

        ! i=6---SEE D1_EXP_5pts.wxm-----------------------------
        i=6
        d1f(i)  =   0.88149995153015d0*(f(i+1)-f(i-1))      / h     &
        -   0.29786830893263d0*(f(i+2)-f(i-2))              / h     &
        +   0.098184076271814d0* (f(i+3) - f(i-3))          / h     &
        -   0.02388959221513d0* (f(i+4) - f(i-4))           / h     &
        +   0.0030486998975861d0* (f(i+5) - f(i-5))         / h

        if (i.eq.imax) return

        ! i=7---SEE D1_EXP_6pts.wxm-----------------------------
        i=7
        d1f(i)  =   0.91595792650492d0  * (f(i+1)-f(i-1))   / h     &
        -   0.34922251636223d0  * (f(i+2)-f(i-2))           / h     &
        +   0.14398145036906d0  * (f(i+3) - f(i-3))         / h     &
        -   0.051236991729043d0 * (f(i+4) - f(i-4))         / h     &
        +   0.013273181125903d0 * (f(i+5) - f(i-5))         / h     &
        -   0.0018126562894445d0* (f(i+6) - f(i-6))         / h

        return

    end subroutine

    subroutine D1c_CTR_CLOSURE_INF_ACC(f, d1f, h, imax)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(nb_of_pts)

        implicit none

        integer :: i, n_max, imax
        real(kind=8) h
        real(kind=8) d1f(:)
        real(kind=8) f(:)

        real(kind=8) d1f_temp(CPT_MIN)


        if (treat_boundary_by_cpt_scheme) then

            d1f_temp(1:CPT_MIN)=d1f(1:CPT_MIN)
            call D1_OUCS3_ACC(f(1:CPT_MIN), d1f_temp(1:CPT_MIN), CPT_MIN, h, .true., Dirichlet)

            d1f(2)=d1f_temp(2)
            if (imax.eq.2) return

            d1f(3)=d1f_temp(3)
            if (imax.eq.3) return

        else

            ! i=2---------------------------------------------------
            i=2
            d1f(i)= d1f(i) + (f(i+1)-f(i-1)) /   (2.d0*h)

            if (i.eq.imax) return
            ! i=3---------------------------------------------------
            i=3
            d1f(i)    =   d1f(i) + 2.d0*(f(i+1)-f(i-1))  / (3.d0*h)      &
            -   (f(i+2) - f(i-2))         / (12.d0*h)

            if (i.eq.imax) return

        end if

        ! i=4---------------------------------------------------
        i=4
        d1f(i)  =   d1f(i) + 0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h

        if (i.eq.imax) return

        ! i=5---------------------------------------------------
        i=5
        d1f(i)  =   d1f(i) + 0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h

        if (i.eq.imax) return

        ! i=6---SEE D1_EXP_5pts.wxm-----------------------------
        i=6
        d1f(i)  =   d1f(i) + 0.88149995153015d0*(f(i+1)-f(i-1))      / h     &
        -   0.29786830893263d0*(f(i+2)-f(i-2))              / h     &
        +   0.098184076271814d0* (f(i+3) - f(i-3))          / h     &
        -   0.02388959221513d0* (f(i+4) - f(i-4))           / h     &
        +   0.0030486998975861d0* (f(i+5) - f(i-5))         / h

        if (i.eq.imax) return

        ! i=7---SEE D1_EXP_6pts.wxm-----------------------------
        i=7
        d1f(i)  =   d1f(i) + 0.91595792650492d0  * (f(i+1)-f(i-1))   / h     &
        -   0.34922251636223d0  * (f(i+2)-f(i-2))           / h     &
        +   0.14398145036906d0  * (f(i+3) - f(i-3))         / h     &
        -   0.051236991729043d0 * (f(i+4) - f(i-4))         / h     &
        +   0.013273181125903d0 * (f(i+5) - f(i-5))         / h     &
        -   0.0018126562894445d0* (f(i+6) - f(i-6))         / h


        return

    end subroutine

    subroutine D1c_CTR_CLOSURE_INF_MULT(f, d1f, h, imax, geom_coefs)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(nb_of_pts)

        implicit none

        integer :: i, n_max, imax
        real(kind=8) h
        real(kind=8) d1f(:), geom_coefs(:)
        real(kind=8) f(:)

        real(kind=8) d1f_temp(CPT_MIN)


        if (treat_boundary_by_cpt_scheme) then

            call D1_OUCS3_MULT( f(1:CPT_MIN), d1f_temp(1:CPT_MIN), CPT_MIN, h, .true.,   &
            Dirichlet, geom_coefs)

            d1f(2)=d1f_temp(2)
            if (imax.eq.2) return

            d1f(3)=d1f_temp(3)
            if (imax.eq.3) return

        else

            ! i=2---------------------------------------------------
            i=2
            d1f(i)= ( (f(i+1)-f(i-1)) /   (2.d0*h) ) * geom_coefs(i)

            if (i.eq.imax) return
            ! i=3---------------------------------------------------
            i=3
            d1f(i)    =   ( 2.d0*(f(i+1)-f(i-1))  / (3.d0*h)      &
            -   (f(i+2) - f(i-2))         / (12.d0*h) ) * geom_coefs(i)

            if (i.eq.imax) return

        endif

        ! i=4---------------------------------------------------
        i=4
        d1f(i)  =   ( 0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h ) * geom_coefs(i)

        if (i.eq.imax) return

        ! i=5---------------------------------------------------
        i=5
        d1f(i)  =  ( 0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h ) * geom_coefs(i)

        if (i.eq.imax) return

        ! i=6---SEE D1_EXP_5pts.wxm-----------------------------
        i=6
        d1f(i)  =   ( 0.88149995153015d0*(f(i+1)-f(i-1))      / h     &
        -   0.29786830893263d0*(f(i+2)-f(i-2))              / h     &
        +   0.098184076271814d0* (f(i+3) - f(i-3))          / h     &
        -   0.02388959221513d0* (f(i+4) - f(i-4))           / h     &
        +   0.0030486998975861d0* (f(i+5) - f(i-5))         / h ) * geom_coefs(i)

        if (i.eq.imax) return

        ! i=7---SEE D1_EXP_6pts.wxm-----------------------------
        i=7
        d1f(i)  =   ( 0.91595792650492d0  * (f(i+1)-f(i-1))   / h     &
        -   0.34922251636223d0  * (f(i+2)-f(i-2))           / h     &
        +   0.14398145036906d0  * (f(i+3) - f(i-3))         / h     &
        -   0.051236991729043d0 * (f(i+4) - f(i-4))         / h     &
        +   0.013273181125903d0 * (f(i+5) - f(i-5))         / h     &
        -   0.0018126562894445d0* (f(i+6) - f(i-6))         / h ) * geom_coefs(i)


        return

    end subroutine

    subroutine D1c_CTR_CLOSURE_INF_MULT_ACC(f, d1f, h, imax, geom_coefs)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(nb_of_pts)

        implicit none

        integer :: i, n_max, imax
        real(kind=8) h
        real(kind=8) d1f(:), geom_coefs(:), d1f_temp(CPT_MIN)
        real(kind=8) f(:)

        if (treat_boundary_by_cpt_scheme) then

            d1f_temp(1:CPT_MIN)=d1f(1:CPT_MIN)
            call D1_OUCS3_MULT_ACC(f(1:CPT_MIN), d1f_temp(1:CPT_MIN), CPT_MIN, h, .true.,    &
            Dirichlet, geom_coefs)

            d1f(2)=d1f_temp(2)
            if (imax.eq.2) return

            d1f(3)=d1f_temp(3)
            if (imax.eq.3) return

        else

            ! i=2---------------------------------------------------
            i=2
            d1f(i)= d1f(i) + ( (f(i+1)-f(i-1)) /   (2.d0*h) ) * geom_coefs(i)

            if (i.eq.imax) return
            ! i=3---------------------------------------------------
            i=3
            d1f(i)    =   d1f(i) + ( 2.d0*(f(i+1)-f(i-1))  / (3.d0*h)      &
            -   (f(i+2) - f(i-2))         / (12.d0*h) ) * geom_coefs(i)

            if (i.eq.imax) return

        endif

        ! i=4---------------------------------------------------
        i=4
        d1f(i)  =   d1f(i) +( 0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h ) * geom_coefs(i)

        if (i.eq.imax) return

        ! i=5---------------------------------------------------
        i=5
        d1f(i)  =   d1f(i) + ( 0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h ) * geom_coefs(i)

        if (i.eq.imax) return

        ! i=6---SEE D1_EXP_5pts.wxm-----------------------------
        i=6
        d1f(i)  =   d1f(i) + ( 0.88149995153015d0*(f(i+1)-f(i-1))      / h     &
        -   0.29786830893263d0*(f(i+2)-f(i-2))              / h     &
        +   0.098184076271814d0* (f(i+3) - f(i-3))          / h     &
        -   0.02388959221513d0* (f(i+4) - f(i-4))           / h     &
        +   0.0030486998975861d0* (f(i+5) - f(i-5))         / h ) * geom_coefs(i)

        if (i.eq.imax) return

        ! i=7---SEE D1_EXP_6pts.wxm-----------------------------
        i=7
        d1f(i)  =   d1f(i) + ( 0.91595792650492d0  * (f(i+1)-f(i-1))   / h     &
        -   0.34922251636223d0  * (f(i+2)-f(i-2))           / h     &
        +   0.14398145036906d0  * (f(i+3) - f(i-3))         / h     &
        -   0.051236991729043d0 * (f(i+4) - f(i-4))         / h     &
        +   0.013273181125903d0 * (f(i+5) - f(i-5))         / h     &
        -   0.0018126562894445d0* (f(i+6) - f(i-6))         / h ) * geom_coefs(i)


        return

    end subroutine

    subroutine D1c_CTR_CLOSURE_SUP(f, d1f, n_max, h, imin, shifted)

        implicit none

        integer :: i, n_max, imin, j
        real(kind=8) h, fn
        real(kind=8) d1f(:)
        real(kind=8) f(:)
        logical :: shifted

        real(kind=8) d1f_temp(n_max-CPT_MIN:n_max-1)

        if (shifted) then
            i=n_max-2
        else
            i=n_max-1
        end if

        if (treat_boundary_by_cpt_scheme) then

            call D1_OUCS3(f(n_max-CPT_MIN:n_max), d1f_temp, CPT_MIN+1, h, shifted, Dirichlet)

            d1f(i)=d1f_temp(i)
            if (i.eq.imin) return

            i=i-1
            d1f(i)=d1f_temp(i)
            if (i.eq.imin) return

        else

            ! First point-------------------------------------------
            d1f(i)= (f(i+1)-f(i-1)) /   (2.d0*h)

            if (i.eq.imin) return
            ! 2nd point----------------------------------------------
            i=i-1
            d1f(i)    =   2.d0*(f(i+1)-f(i-1))  / (3.d0*h)      &
            -   (f(i+2) - f(i-2))         / (12.d0*h)

            if (i.eq.imin) return

        endif

        ! 3rd point---------------------------------------------
        i=i-1
        d1f(i)  =   0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h

        if (i.eq.imin) return

        ! 4th point---------------------------------------------
        i=i-1
        d1f(i)  =   0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h

        if (i.eq.imin)return

        ! 5th point ----D1_EXP_5pts.wxm------------------------
        i=i-1
        d1f(i)  =   0.88149995153015d0*(f(i+1)-f(i-1))      / h     &
        -   0.29786830893263d0*(f(i+2)-f(i-2))              / h     &
        +   0.098184076271814d0* (f(i+3) - f(i-3))          / h     &
        -   0.02388959221513d0* (f(i+4) - f(i-4))           / h     &
        +   0.0030486998975861d0* (f(i+5) - f(i-5))         / h

        if (i.eq.imin) return

        ! 6th point ----D1_EXP_6pts.wxm------------------------
        i=i-1
        d1f(i)  =   0.91595792650492d0  * (f(i+1)-f(i-1))   / h     &
        -   0.34922251636223d0  * (f(i+2)-f(i-2))           / h     &
        +   0.14398145036906d0  * (f(i+3) - f(i-3))         / h     &
        -   0.051236991729043d0 * (f(i+4) - f(i-4))         / h     &
        +   0.013273181125903d0 * (f(i+5) - f(i-5))         / h     &
        -   0.0018126562894445d0* (f(i+6) - f(i-6))         / h

        if (i.eq.imin) return

        ! 7th point ---------------------------------------------
        i=i-1
        d1f(i)  =   0.95396219562045d0  * (f(i+1)-f(i-1))   / h     &
        -   0.41234590494721d0  * (f(i+2)-f(i-2))           / h     &
        +   0.21233981563217d0  * (f(i+3) - f(i-3))         / h     &
        -   0.10672135533957d0 * (f(i+4) - f(i-4))          / h     &
        +   0.046775295199319d0 * (f(i+5) - f(i-5))         / h     &
        -   0.015323662211638d0* (f(i+6) - f(i-6))          / h     &
        +   0.0026664396358809d0* (f(i+7) - f(i-7))          / h


        return

    end subroutine

    subroutine D1c_CTR_CLOSURE_SUP_ACC(f, d1f, n_max, h, imin, shifted)

        implicit none

        integer :: i, n_max, imin, j
        real(kind=8) h, fn
        real(kind=8) d1f(:)
        real(kind=8) f(:)
        logical :: shifted

        real(kind=8) d1f_temp(n_max-CPT_MIN:n_max-1)

        if (shifted) then
            i=n_max-2
        else
            i=n_max-1
        end if

        if (treat_boundary_by_cpt_scheme) then
            d1f_temp=d1f(n_max-CPT_MIN:n_max-1)
            call D1_OUCS3_ACC(f(n_max-CPT_MIN:n_max), d1f_temp, CPT_MIN+1, h, shifted, Dirichlet)

            d1f(i)=d1f_temp(i)
            if (i.eq.imin) return

            i=i-1
            d1f(i)=d1f_temp(i)
            if (i.eq.imin) return

        else

            ! First point-------------------------------------------
            d1f(i)= d1f(i) + (f(i+1)-f(i-1)) /   (2.d0*h)

            if (i.eq.imin) return
            ! 2nd point----------------------------------------------
            i=i-1
            d1f(i)    =   d1f(i) + 2.d0*(f(i+1)-f(i-1))  / (3.d0*h)      &
            -   (f(i+2) - f(i-2))         / (12.d0*h)

            if (i.eq.imin) return

        endif

        ! 3rd point---------------------------------------------
        i=i-1
        d1f(i)  =   d1f(i) + 0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h

        if (i.eq.imin) return

        ! 4th point---------------------------------------------
        i=i-1
        d1f(i)  =   d1f(i) + 0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h

        if (i.eq.imin) return

        ! 5th point ----D1_EXP_5pts.wxm------------------------
        i=i-1
        d1f(i)  =   d1f(i) + 0.88149995153015d0*(f(i+1)-f(i-1))      / h     &
        -   0.29786830893263d0*(f(i+2)-f(i-2))              / h     &
        +   0.098184076271814d0* (f(i+3) - f(i-3))          / h     &
        -   0.02388959221513d0* (f(i+4) - f(i-4))           / h     &
        +   0.0030486998975861d0* (f(i+5) - f(i-5))         / h

        if (i.eq.imin) return

        ! 6th point ----D1_EXP_6pts.wxm------------------------
        i=i-1
        d1f(i)  =   d1f(i) + 0.91595792650492d0  * (f(i+1)-f(i-1))   / h     &
        -   0.34922251636223d0  * (f(i+2)-f(i-2))           / h     &
        +   0.14398145036906d0  * (f(i+3) - f(i-3))         / h     &
        -   0.051236991729043d0 * (f(i+4) - f(i-4))         / h     &
        +   0.013273181125903d0 * (f(i+5) - f(i-5))         / h     &
        -   0.0018126562894445d0* (f(i+6) - f(i-6))         / h

        if (i.eq.imin) return

        ! 7th point ----D1_EXP_6pts.wxm------------------------
        i=i-1
        d1f(i)  =   d1f(i) + 0.95396219562045d0  * (f(i+1)-f(i-1))   / h     &
        -   0.41234590494721d0  * (f(i+2)-f(i-2))           / h     &
        +   0.21233981563217d0  * (f(i+3) - f(i-3))         / h     &
        -   0.10672135533957d0 * (f(i+4) - f(i-4))          / h     &
        +   0.046775295199319d0 * (f(i+5) - f(i-5))         / h     &
        -   0.015323662211638d0* (f(i+6) - f(i-6))          / h     &
        +   0.0026664396358809d0* (f(i+7) - f(i-7))          / h


        return

    end subroutine

    subroutine D1c_CTR_CLOSURE_SUP_MULT(f, d1f, n_max, h, imin, shifted, geom_coefs)

        implicit none

        integer :: i, n_max, imin, j
        real(kind=8) h, fn
        real(kind=8) d1f(:), geom_coefs(:)
        real(kind=8) f(:)
        logical :: shifted

        real(kind=8) d1f_temp(n_max-CPT_MIN:n_max-1)

        if (shifted) then
            i=n_max-2
        else
            i=n_max-1
        end if

        if (treat_boundary_by_cpt_scheme) then

            call D1_OUCS3_MULT(f(n_max-CPT_MIN:n_max), d1f_temp, CPT_MIN+1, h, shifted,     &
            Dirichlet, geom_coefs(n_max-CPT_MIN:n_max))

            d1f(i)=d1f_temp(i)
            if (i.eq.imin) return

            i=i-1
            d1f(i)=d1f_temp(i)
            if (i.eq.imin) return

        else

            ! First point-------------------------------------------
            d1f(i)= ( (f(i+1)-f(i-1)) /   (2.d0*h) ) * geom_coefs(i)

            if (i.eq.imin) return
            ! 2nd point----------------------------------------------
            i=i-1
            d1f(i)    =   ( 2.d0*(f(i+1)-f(i-1))  / (3.d0*h)      &
            -   (f(i+2) - f(i-2))         / (12.d0*h) ) * geom_coefs(i)

            if (i.eq.imin) return

        endif

        ! 3rd point---------------------------------------------
        i=i-1
        d1f(i)  =   ( 0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h ) * geom_coefs(i)

        if (i.eq.imin) return

        ! 4th point---------------------------------------------
        i=i-1
        d1f(i)  =   ( 0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h ) * geom_coefs(i)

        if (i.eq.imin) return

        ! 5th point ----D1_EXP_5pts.wxm------------------------
        i=i-1
        d1f(i)  =   ( 0.88149995153015d0*(f(i+1)-f(i-1))      / h     &
        -   0.29786830893263d0*(f(i+2)-f(i-2))              / h     &
        +   0.098184076271814d0* (f(i+3) - f(i-3))          / h     &
        -   0.02388959221513d0* (f(i+4) - f(i-4))           / h     &
        +   0.0030486998975861d0* (f(i+5) - f(i-5))         / h ) * geom_coefs(i)

        if (i.eq.imin) return

        ! 6th point ----D1_EXP_6pts.wxm------------------------
        i=i-1
        d1f(i)  =   ( 0.91595792650492d0  * (f(i+1)-f(i-1))   / h     &
        -   0.34922251636223d0  * (f(i+2)-f(i-2))           / h     &
        +   0.14398145036906d0  * (f(i+3) - f(i-3))         / h     &
        -   0.051236991729043d0 * (f(i+4) - f(i-4))         / h     &
        +   0.013273181125903d0 * (f(i+5) - f(i-5))         / h     &
        -   0.0018126562894445d0* (f(i+6) - f(i-6))         / h ) * geom_coefs(i)

        if (i.eq.imin) return

        ! 7th point ----D1_EXP_6pts.wxm------------------------
        i=i-1
        d1f(i)  =   ( 0.95396219562045d0  * (f(i+1)-f(i-1))   / h     &
        -   0.41234590494721d0  * (f(i+2)-f(i-2))           / h     &
        +   0.21233981563217d0  * (f(i+3) - f(i-3))         / h     &
        -   0.10672135533957d0 * (f(i+4) - f(i-4))          / h     &
        +   0.046775295199319d0 * (f(i+5) - f(i-5))         / h     &
        -   0.015323662211638d0* (f(i+6) - f(i-6))          / h     &
        +   0.0026664396358809d0* (f(i+7) - f(i-7))          / h ) * geom_coefs(i)


        return

    end subroutine

    subroutine D1c_CTR_CLOSURE_SUP_MULT_ACC(f, d1f, n_max, h, imin, shifted, geom_coefs)

        implicit none

        integer :: i, n_max, imin, j
        real(kind=8) h, fn
        real(kind=8) d1f(:), geom_coefs(:)
        real(kind=8) f(:)
        logical :: shifted

        real(kind=8) d1f_temp(n_max-CPT_MIN:n_max-1)

        if (shifted) then
            i=n_max-2
        else
            i=n_max-1
        end if

        if (treat_boundary_by_cpt_scheme) then
            d1f_temp=d1f(n_max-CPT_MIN:n_max-1)
            call D1_OUCS3_MULT_ACC(f(n_max-CPT_MIN:n_max), d1f_temp, CPT_MIN+1, h, shifted,     &
            Dirichlet, geom_coefs(n_max-CPT_MIN:n_max))

            d1f(i)=d1f_temp(i)
            if (i.eq.imin) return

            i=i-1
            d1f(i)=d1f_temp(i)
            if (i.eq.imin) return

        else

            ! First point-------------------------------------------
            d1f(i)= d1f(i) + ( (f(i+1)-f(i-1)) /   (2.d0*h) ) * geom_coefs(i)

            if (i.eq.imin) return
            ! 2nd point----------------------------------------------
            i=i-1
            d1f(i)    =   d1f(i) + ( 2.d0*(f(i+1)-f(i-1))  / (3.d0*h)      &
            -   (f(i+2) - f(i-2))         / (12.d0*h) ) * geom_coefs(i)

            if (i.eq.imin) return

        endif

        ! 3rd point---------------------------------------------
        i=i-1
        d1f(i)  =   d1f(i) + ( 0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h ) * geom_coefs(i)

        if (i.eq.imin) return

        ! 4th point---------------------------------------------
        i=i-1
        d1f(i)  =   d1f(i) + ( 0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h ) * geom_coefs(i)

        if (i.eq.imin) return

        ! 5th point ----D1_EXP_5pts.wxm------------------------
        i=i-1
        d1f(i)  =   d1f(i) + ( 0.88149995153015d0*(f(i+1)-f(i-1))      / h     &
        -   0.29786830893263d0*(f(i+2)-f(i-2))              / h     &
        +   0.098184076271814d0* (f(i+3) - f(i-3))          / h     &
        -   0.02388959221513d0* (f(i+4) - f(i-4))           / h     &
        +   0.0030486998975861d0* (f(i+5) - f(i-5))         / h ) * geom_coefs(i)

        if (i.eq.imin) return

        ! 6th point ----D1_EXP_6pts.wxm------------------------
        i=i-1
        d1f(i)  =   d1f(i) + ( 0.91595792650492d0  * (f(i+1)-f(i-1))   / h     &
        -   0.34922251636223d0  * (f(i+2)-f(i-2))           / h     &
        +   0.14398145036906d0  * (f(i+3) - f(i-3))         / h     &
        -   0.051236991729043d0 * (f(i+4) - f(i-4))         / h     &
        +   0.013273181125903d0 * (f(i+5) - f(i-5))         / h     &
        -   0.0018126562894445d0* (f(i+6) - f(i-6))         / h ) * geom_coefs(i)

        if (i.eq.imin) return

        ! 7th point ----D1_EXP_6pts.wxm------------------------
        i=i-1
        d1f(i)  =   d1f(i) + ( 0.95396219562045d0  * (f(i+1)-f(i-1))   / h     &
        -   0.41234590494721d0  * (f(i+2)-f(i-2))           / h     &
        +   0.21233981563217d0  * (f(i+3) - f(i-3))         / h     &
        -   0.10672135533957d0 * (f(i+4) - f(i-4))          / h     &
        +   0.046775295199319d0 * (f(i+5) - f(i-5))         / h     &
        -   0.015323662211638d0* (f(i+6) - f(i-6))          / h     &
        +   0.0026664396358809d0* (f(i+7) - f(i-7))          / h ) * geom_coefs(i)


        return

    end subroutine

end module

module d2c_schemes

    use boundaries_types
    use schemes_settings

    implicit none

contains

    subroutine D2c_ExpCtr_O2Fp0(f, d2f, n_max, h, shifted, bc_type)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h
        real(kind=8), intent(in)      :: f(:)
        real(kind=8), intent(out)     :: d2f(:)

        integer :: i
        real(kind=8) A, inv_h2

        A=0.5d0/h
        inv_h2=1.d0/h**2

        do i=2,n_max-2
            d2f(i)= ( f(i-1) - f(i)*2.d0 + f(i+1) ) *inv_h2
        enddo

        if (bc_type.eq.periodic) then

            d2f(1)= ( f(n_max-1) - f(1)*2.d0 + f(2) ) *inv_h2
            d2f(n_max-1)=( f(n_max-2) - f(n_max-1)*2.d0 + f(1) ) *inv_h2

        end if

        if (bc_type.eq.Dirichlet) then

            if (shifted.eqv.(.false.)) then
                d2f(n_max-1)= ( f(n_max-2) - f(n_max-1)*2.d0 + f(n_max) ) *inv_h2
            end if
        end if

        return

    end subroutine

    subroutine D2c_ExpCtr_O2Fp0_ACC(f, d2f, n_max, h, shifted, bc_type)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h
        real(kind=8), intent(in)      ::f(:)
        real(kind=8), intent(inout)   :: d2f(:)

        integer :: i
        real(kind=8) A, inv_h2

        A=0.5d0/h
        inv_h2=1.d0/h**2

        do i=2,n_max-2
            d2f(i)= d2f(i) + ( f(i-1) - f(i)*2.d0 + f(i+1) ) *inv_h2
        enddo

        if (bc_type.eq.periodic) then

            d2f(1)= d2f(1) + ( f(n_max-1) - f(1)*2.d0 + f(2) ) *inv_h2
            d2f(n_max-1)= d2f(n_max-1) + ( f(n_max-2) - f(n_max-1)*2.d0 + f(1) ) *inv_h2

        end if

        if (bc_type.eq.Dirichlet) then

            if (shifted.eqv.(.false.)) then
                d2f(n_max-1)= d2f(n_max-1) + ( f(n_max-2) - f(n_max-1)*2.d0 + f(n_max) ) *inv_h2
            end if
        end if

        return

    end subroutine

    subroutine D2c_ExpCtr_O2Fp0_MULT(f, d2f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h,  geom_coefs(:)
        real(kind=8), intent(in)      ::f(:)
        real(kind=8), intent(out)   :: d2f(:)

        integer :: i
        real(kind=8) A, inv_h2

        A=0.5d0/h
        inv_h2=1.d0/h**2

        do i=2,n_max-2
            d2f(i)= ( f(i-1) - f(i)*2.d0 + f(i+1) ) *inv_h2 * geom_coefs(i)
        enddo

        if (bc_type.eq.periodic) then

            d2f(1)= ( f(n_max-1) - f(1)*2.d0 + f(2) ) *inv_h2 * geom_coefs(1)
            d2f(n_max-1)= ( f(n_max-2) - f(n_max-1)*2.d0 + f(1) ) *inv_h2 * geom_coefs(n_max-1)

        end if

        if (bc_type.eq.Dirichlet) then

            if (shifted.eqv.(.false.)) then
                d2f(n_max-1)= ( f(n_max-2) - f(n_max-1)*2.d0 + f(n_max) ) *inv_h2 * geom_coefs(n_max-1)
            end if
        end if

        return

    end subroutine

    subroutine D2c_ExpCtr_O2Fp0_MULT_ACC(f, d2f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h,  geom_coefs(:)
        real(kind=8), intent(in)      ::f(:)
        real(kind=8), intent(inout)   :: d2f(:)

        integer :: i
        real(kind=8) A, inv_h2

        A=0.5d0/h
        inv_h2=1.d0/h**2

        do i=2,n_max-2
            d2f(i)= d2f(i) + ( f(i-1) - f(i)*2.d0 + f(i+1) ) *inv_h2 * geom_coefs(i)
        enddo

        if (bc_type.eq.periodic) then

            d2f(1)= d2f(1) + ( f(n_max-1) - f(1)*2.d0 + f(2) ) *inv_h2 * geom_coefs(1)
            d2f(n_max-1)= d2f(n_max-1) + ( f(n_max-2) - f(n_max-1)*2.d0 + f(1) ) *inv_h2 * geom_coefs(n_max-1)

        end if

        if (bc_type.eq.Dirichlet) then

            if (shifted.eqv.(.false.)) then
                d2f(n_max-1)= d2f(n_max-1) + ( f(n_max-2) - f(n_max-1)*2.d0 + f(n_max) ) *inv_h2 * geom_coefs(n_max-1)
            end if
        end if

        return

    end subroutine

    subroutine D2c_Tamm_twice(f, d2f, n_max, h, shifted, bc_type)


        use d1s_schemes

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h
        real(kind=8), intent(in)      :: f(:)
        real(kind=8), intent(out)     :: d2f(:)

        real(kind=8):: d1f(n_max)

        call D1s_ExpCtr_O0Fp5(f, d1f, n_max, h, shifted, bc_type)

        if (shifted) then

            if (bc_type.eq.periodic) then
                call D1s_Tamm(d1f, d2f, n_max, h, .false., periodic)
            else
                call D1s_Tamm(d1f(2:n_max-1), d2f(2:n_max-2), n_max-2, h, .false., Dirichlet)
            end if

        else

            call D1s_Tamm(d1f, d2f, n_max, h, .true., bc_type)

        end if

        return

    end subroutine

    subroutine D2c_Tamm_twice_ACC(f, d2f, n_max, h, shifted, bc_type)


        use d1s_schemes

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h
        real(kind=8), intent(in)      ::f(:)
        real(kind=8), intent(inout)   :: d2f(:)

        real(kind=8):: d1f(n_max)

        call D1s_ExpCtr_O0Fp5(f, d1f, n_max, h, shifted, bc_type)

        if (shifted) then

            if (bc_type.eq.periodic) then
                call D1s_Tamm_ACC(d1f, d2f, n_max, h, .false., periodic)
            else
                call D1s_Tamm_ACC(d1f(2:n_max-1), d2f(2:n_max-2), n_max-2, h, .false., Dirichlet)
            end if

        else

            call D1s_Tamm_ACC(d1f, d2f, n_max, h, .true., bc_type)

        end if

        return

    end subroutine

    subroutine D2c_Tamm_twice_MULT(f, d2f, n_max, h, shifted, bc_type, geom_coefs)


        use d1s_schemes

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h,  geom_coefs(:)
        real(kind=8), intent(in)      ::f(:)
        real(kind=8), intent(out)   :: d2f(:)

        real(kind=8):: d1f(n_max)

        call D1s_ExpCtr_O0Fp5(f, d1f, n_max, h, shifted, bc_type)

        if (shifted) then

            if (bc_type.eq.periodic) then
                call D1s_Tamm_MULT(d1f, d2f, n_max, h, .false., periodic, geom_coefs)
            else
                call D1s_Tamm_MULT(d1f(2:n_max-1), d2f(2:n_max-2), n_max-2, h, .false., Dirichlet, geom_coefs(2:n_max-2))
            end if

        else

            call D1s_Tamm_MULT(d1f, d2f, n_max, h, .true., bc_type, geom_coefs)

        end if

        return

    end subroutine

    subroutine D2c_Tamm_twice_MULT_ACC(f, d2f, n_max, h, shifted, bc_type, geom_coefs)


        use d1s_schemes

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h,  geom_coefs(:)
        real(kind=8), intent(in)      ::f(:)
        real(kind=8), intent(inout)   :: d2f(:)

        real(kind=8):: d1f(n_max)

        call D1s_ExpCtr_O0Fp5(f, d1f, n_max, h, shifted, bc_type)

        if (shifted) then

            if (bc_type.eq.periodic) then
                call D1s_Tamm_MULT_ACC(d1f, d2f, n_max, h, .false., periodic, geom_coefs)
            else
                call D1s_Tamm_MULT_ACC(d1f(2:n_max-1), d2f(2:n_max-2), n_max-2, h, .false., Dirichlet, geom_coefs(2:n_max-2))
            end if

        else

            call D1s_Tamm_MULT_ACC(d1f, d2f, n_max, h, .true., bc_type, geom_coefs)

        end if

        return

    end subroutine

    subroutine D2_ExpCtr_O0Fp6(f, d2f, n_max, h, shifted, bc_type)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h
        real(kind=8), intent(in)      :: f(:)
        real(kind=8), intent(out)     :: d2f(:)

        integer :: i
        real(kind=8) A(0:6)

        A(1) =  1.813617457417053d0         / h**2
        A(2) =  -0.33529477396055d0         / h**2
        A(3) =  0.087416628525326d0         / h**2
        A(4) =  -0.021599958070134d0        / h**2
        A(5) =  0.0040422194201875d0        / h**2
        A(6) = -4.0678692104879698d-4       / h**2
        A(0) = -2.d0*(A(1) + A(2) + A(3) + A(4) + A(5) + A(6))

        do i=7, n_max-7
            d2f(i)=   A(6)*(f(i+6) + f(i-6))    &
            + A(5)*(f(i+5) + f(i-5))    &
            + A(4)*(f(i+4) + f(i-4))    &
            + A(3)*(f(i+3) + f(i-3))    &
            + A(2)*(f(i+2) + f(i-2))    &
            + A(1)*(f(i+1) + f(i-1))    &
            + A(0)*f(i)
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 6
                d2f(i)=   A(6)*(f(i+6) + f(mod(n_max+i-8,n_max-1)+1))     &
                + A(5)*(f(i+5) + f(mod(n_max+i-7,n_max-1)+1))       &
                + A(4)*(f(i+4) + f(mod(n_max+i-6,n_max-1)+1))       &
                + A(3)*(f(i+3) + f(mod(n_max+i-5,n_max-1)+1))       &
                + A(2)*(f(i+2) + f(mod(n_max+i-4,n_max-1)+1))       &
                + A(1)*(f(i+1) + f(mod(n_max+i-3,n_max-1)+1))       &
                + A(0)*f(i)
            enddo

            do i=n_max-6, n_max-1
                d2f(i)=   A(6)*(f(mod(i+5,n_max-1)+1) + f(i-6))       &
                + A(5)*(f(mod(i+4,n_max-1)+1) + f(i-5))         &
                + A(4)*(f(mod(i+3,n_max-1)+1) + f(i-4))         &
                + A(3)*(f(mod(i+2,n_max-1)+1) + f(i-3))         &
                + A(2)*(f(mod(i+1,n_max-1)+1) + f(i-2))         &
                + A(1)*(f(mod(i,n_max-1)+1) + f(i-1))           &
                + A(0)*f(i)
            enddo

        end if

        if (bc_type.eq.Dirichlet) then
            call D2c_CTR_CLOSURE_INF(f, d2f, h, 6)
        end if

        if (bc_type.eq.Dirichlet) then
            call D2c_CTR_CLOSURE_SUP(f, d2f, n_max, h, n_max-6, shifted)
        end if

        return

    end subroutine

    subroutine D2_ExpCtr_O0Fp6_ACC(f, d2f, n_max, h, shifted, bc_type)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h
        real(kind=8), intent(in)      ::f(:)
        real(kind=8), intent(inout)   :: d2f(:)

        integer :: i
        real(kind=8) A(0:6)

        A(1) =  1.813617457417053d0         / h**2
        A(2) =  -0.33529477396055d0         / h**2
        A(3) =  0.087416628525326d0         / h**2
        A(4) =  -0.021599958070134d0        / h**2
        A(5) =  0.0040422194201875d0        / h**2
        A(6) = -4.0678692104879698d-4       / h**2
        A(0) = -2.d0*(A(1) + A(2) + A(3) + A(4) + A(5) + A(6))

        do i=7, n_max-7
            d2f(i)=   d2f(i) + A(6)*(f(i+6) + f(i-6))    &
            + A(5)*(f(i+5) + f(i-5))    &
            + A(4)*(f(i+4) + f(i-4))    &
            + A(3)*(f(i+3) + f(i-3))    &
            + A(2)*(f(i+2) + f(i-2))    &
            + A(1)*(f(i+1) + f(i-1))    &
            + A(0)*f(i)
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 6
                d2f(i)=   d2f(i) + A(6)*(f(i+6) + f(mod(n_max+i-8,n_max-1)+1))     &
                + A(5)*(f(i+5) + f(mod(n_max+i-7,n_max-1)+1))       &
                + A(4)*(f(i+4) + f(mod(n_max+i-6,n_max-1)+1))       &
                + A(3)*(f(i+3) + f(mod(n_max+i-5,n_max-1)+1))       &
                + A(2)*(f(i+2) + f(mod(n_max+i-4,n_max-1)+1))       &
                + A(1)*(f(i+1) + f(mod(n_max+i-3,n_max-1)+1))       &
                + A(0)*f(i)
            enddo

            do i=n_max-6, n_max-1
                d2f(i)=   d2f(i) + A(6)*(f(mod(i+5,n_max-1)+1) + f(i-6))       &
                + A(5)*(f(mod(i+4,n_max-1)+1) + f(i-5))         &
                + A(4)*(f(mod(i+3,n_max-1)+1) + f(i-4))         &
                + A(3)*(f(mod(i+2,n_max-1)+1) + f(i-3))         &
                + A(2)*(f(mod(i+1,n_max-1)+1) + f(i-2))         &
                + A(1)*(f(mod(i,n_max-1)+1) + f(i-1))           &
                + A(0)*f(i)
            enddo

        end if

        if (bc_type.eq.Dirichlet) then
            call D2c_CTR_CLOSURE_INF_ACC(f, d2f, h, 6)
        end if

        if (bc_type.eq.Dirichlet) then
            call D2c_CTR_CLOSURE_SUP_ACC(f, d2f, n_max, h, n_max-6, shifted)
        end if

        return

    end subroutine

    subroutine D2_ExpCtr_O0Fp6_MULT(f, d2f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h,  geom_coefs(:)
        real(kind=8), intent(in)      ::f(:)
        real(kind=8), intent(out)   :: d2f(:)

        integer :: i
        real(kind=8) A(0:6)

        A(1) =  1.813617457417053d0         / h**2
        A(2) =  -0.33529477396055d0         / h**2
        A(3) =  0.087416628525326d0         / h**2
        A(4) =  -0.021599958070134d0        / h**2
        A(5) =  0.0040422194201875d0        / h**2
        A(6) = -4.0678692104879698d-4       / h**2
        A(0) = -2.d0*(A(1) + A(2) + A(3) + A(4) + A(5) + A(6))

        do i=7, n_max-7
            d2f(i)=   ( A(6)*(f(i+6) + f(i-6))    &
            + A(5)*(f(i+5) + f(i-5))    &
            + A(4)*(f(i+4) + f(i-4))    &
            + A(3)*(f(i+3) + f(i-3))    &
            + A(2)*(f(i+2) + f(i-2))    &
            + A(1)*(f(i+1) + f(i-1))    &
            + A(0)*f(i) ) * geom_coefs(i)
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 6
                d2f(i)=   ( A(6)*(f(i+6) + f(mod(n_max+i-8,n_max-1)+1))     &
                + A(5)*(f(i+5) + f(mod(n_max+i-7,n_max-1)+1))       &
                + A(4)*(f(i+4) + f(mod(n_max+i-6,n_max-1)+1))       &
                + A(3)*(f(i+3) + f(mod(n_max+i-5,n_max-1)+1))       &
                + A(2)*(f(i+2) + f(mod(n_max+i-4,n_max-1)+1))       &
                + A(1)*(f(i+1) + f(mod(n_max+i-3,n_max-1)+1))       &
                + A(0)*f(i) ) * geom_coefs(i)
            enddo

            do i=n_max-6, n_max-1
                d2f(i)=   ( A(6)*(f(mod(i+5,n_max-1)+1) + f(i-6))       &
                + A(5)*(f(mod(i+4,n_max-1)+1) + f(i-5))         &
                + A(4)*(f(mod(i+3,n_max-1)+1) + f(i-4))         &
                + A(3)*(f(mod(i+2,n_max-1)+1) + f(i-3))         &
                + A(2)*(f(mod(i+1,n_max-1)+1) + f(i-2))         &
                + A(1)*(f(mod(i,n_max-1)+1) + f(i-1))           &
                + A(0)*f(i) ) * geom_coefs(i)
            enddo

        end if

        if (bc_type.eq.Dirichlet) then
            call D2c_CTR_CLOSURE_INF_MULT(f, d2f, h, 6, geom_coefs)
        end if

        if (bc_type.eq.Dirichlet) then
            call D2c_CTR_CLOSURE_SUP_MULT(f, d2f, n_max, h, n_max-6, shifted, geom_coefs)
        end if

        return

    end subroutine

    subroutine D2_ExpCtr_O0Fp6_MULT_ACC(f, d2f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h,  geom_coefs(:)
        real(kind=8), intent(in)      ::f(:)
        real(kind=8), intent(inout)   :: d2f(:)

        integer :: i
        real(kind=8) A(0:6)

        A(1) =  1.813617457417053d0         / h**2
        A(2) =  -0.33529477396055d0         / h**2
        A(3) =  0.087416628525326d0         / h**2
        A(4) =  -0.021599958070134d0        / h**2
        A(5) =  0.0040422194201875d0        / h**2
        A(6) = -4.0678692104879698d-4       / h**2
        A(0) = -2.d0*(A(1) + A(2) + A(3) + A(4) + A(5) + A(6))

        do i=7, n_max-7
            d2f(i)=   d2f(i) + ( A(6)*(f(i+6) + f(i-6))    &
            + A(5)*(f(i+5) + f(i-5))    &
            + A(4)*(f(i+4) + f(i-4))    &
            + A(3)*(f(i+3) + f(i-3))    &
            + A(2)*(f(i+2) + f(i-2))    &
            + A(1)*(f(i+1) + f(i-1))    &
            + A(0)*f(i) ) * geom_coefs(i)
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 6
                d2f(i)=   d2f(i) + ( A(6)*(f(i+6) + f(mod(n_max+i-8,n_max-1)+1))     &
                + A(5)*(f(i+5) + f(mod(n_max+i-7,n_max-1)+1))       &
                + A(4)*(f(i+4) + f(mod(n_max+i-6,n_max-1)+1))       &
                + A(3)*(f(i+3) + f(mod(n_max+i-5,n_max-1)+1))       &
                + A(2)*(f(i+2) + f(mod(n_max+i-4,n_max-1)+1))       &
                + A(1)*(f(i+1) + f(mod(n_max+i-3,n_max-1)+1))       &
                + A(0)*f(i) ) * geom_coefs(i)
            enddo

            do i=n_max-6, n_max-1
                d2f(i)=   d2f(i) + ( A(6)*(f(mod(i+5,n_max-1)+1) + f(i-6))       &
                + A(5)*(f(mod(i+4,n_max-1)+1) + f(i-5))         &
                + A(4)*(f(mod(i+3,n_max-1)+1) + f(i-4))         &
                + A(3)*(f(mod(i+2,n_max-1)+1) + f(i-3))         &
                + A(2)*(f(mod(i+1,n_max-1)+1) + f(i-2))         &
                + A(1)*(f(mod(i,n_max-1)+1) + f(i-1))           &
                + A(0)*f(i) ) * geom_coefs(i)
            enddo

        end if

        if (bc_type.eq.Dirichlet) then
            call D2c_CTR_CLOSURE_INF_MULT_ACC(f, d2f, h, 6, geom_coefs)
        end if

        if (bc_type.eq.Dirichlet) then
            call D2c_CTR_CLOSURE_SUP_MULT_ACC(f, d2f, n_max, h, n_max-6, shifted, geom_coefs)
        end if

        return

    end subroutine

    subroutine D2_ExpCtr_O0Fp6_hump(f, d2f, n_max, h, shifted, bc_type)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h
        real(kind=8), intent(in)      :: f(:)
        real(kind=8), intent(out)     :: d2f(:)

        integer :: i
        real(kind=8) A(0:6)

        A(1) =  2.178831413991164d0     / h**2
        A(2) =  -0.62233954789826d0     / h**2
        A(3) =  0.27324369530023d0      / h**2
        A(4) =  -0.11423388546519d0     / h**2
        A(5) =  0.035021195634051d0     / h**2
        A(6) = -0.0054574850439034d0    / h**2
        A(0) = -2.d0*(A(1) + A(2) + A(3) + A(4) + A(5) + A(6))

        do i=7, n_max-7
            d2f(i)=   A(6)*(f(i+6) + f(i-6))    &
            + A(5)*(f(i+5) + f(i-5))    &
            + A(4)*(f(i+4) + f(i-4))    &
            + A(3)*(f(i+3) + f(i-3))    &
            + A(2)*(f(i+2) + f(i-2))    &
            + A(1)*(f(i+1) + f(i-1))    &
            + A(0)*f(i)
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 6
                d2f(i)=   A(6)*(f(i+6) + f(mod(n_max+i-8,n_max-1)+1))     &
                + A(5)*(f(i+5) + f(mod(n_max+i-7,n_max-1)+1))       &
                + A(4)*(f(i+4) + f(mod(n_max+i-6,n_max-1)+1))       &
                + A(3)*(f(i+3) + f(mod(n_max+i-5,n_max-1)+1))       &
                + A(2)*(f(i+2) + f(mod(n_max+i-4,n_max-1)+1))       &
                + A(1)*(f(i+1) + f(mod(n_max+i-3,n_max-1)+1))       &
                + A(0)*f(i)
            enddo

            do i=n_max-6, n_max-1
                d2f(i)=   A(6)*(f(mod(i+5,n_max-1)+1) + f(i-6))       &
                + A(5)*(f(mod(i+4,n_max-1)+1) + f(i-5))         &
                + A(4)*(f(mod(i+3,n_max-1)+1) + f(i-4))         &
                + A(3)*(f(mod(i+2,n_max-1)+1) + f(i-3))         &
                + A(2)*(f(mod(i+1,n_max-1)+1) + f(i-2))         &
                + A(1)*(f(mod(i,n_max-1)+1) + f(i-1))           &
                + A(0)*f(i)
            enddo

        end if

        if (bc_type.eq.Dirichlet) then
            call D2c_CTR_CLOSURE_INF(f, d2f, h, 6)
        end if

        if (bc_type.eq.Dirichlet) then
            call D2c_CTR_CLOSURE_SUP(f, d2f, n_max, h, n_max-6, shifted)
        end if

        return

    end subroutine

    subroutine D2_ExpCtr_O0Fp6_hump_ACC(f, d2f, n_max, h, shifted, bc_type)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h
        real(kind=8), intent(in)      ::f(:)
        real(kind=8), intent(inout)   :: d2f(:)

        integer :: i
        real(kind=8) A(0:6)

        A(1) =  2.178831413991164d0     / h**2
        A(2) =  -0.62233954789826d0     / h**2
        A(3) =  0.27324369530023d0      / h**2
        A(4) =  -0.11423388546519d0     / h**2
        A(5) =  0.035021195634051d0     / h**2
        A(6) = -0.0054574850439034d0    / h**2
        A(0) = -2.d0*(A(1) + A(2) + A(3) + A(4) + A(5) + A(6))

        do i=7, n_max-7
            d2f(i)=   d2f(i) + A(6)*(f(i+6) + f(i-6))    &
            + A(5)*(f(i+5) + f(i-5))    &
            + A(4)*(f(i+4) + f(i-4))    &
            + A(3)*(f(i+3) + f(i-3))    &
            + A(2)*(f(i+2) + f(i-2))    &
            + A(1)*(f(i+1) + f(i-1))    &
            + A(0)*f(i)
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 6
                d2f(i)=   d2f(i) + A(6)*(f(i+6) + f(mod(n_max+i-8,n_max-1)+1))     &
                + A(5)*(f(i+5) + f(mod(n_max+i-7,n_max-1)+1))       &
                + A(4)*(f(i+4) + f(mod(n_max+i-6,n_max-1)+1))       &
                + A(3)*(f(i+3) + f(mod(n_max+i-5,n_max-1)+1))       &
                + A(2)*(f(i+2) + f(mod(n_max+i-4,n_max-1)+1))       &
                + A(1)*(f(i+1) + f(mod(n_max+i-3,n_max-1)+1))       &
                + A(0)*f(i)
            enddo

            do i=n_max-6, n_max-1
                d2f(i)=   d2f(i) + A(6)*(f(mod(i+5,n_max-1)+1) + f(i-6))       &
                + A(5)*(f(mod(i+4,n_max-1)+1) + f(i-5))         &
                + A(4)*(f(mod(i+3,n_max-1)+1) + f(i-4))         &
                + A(3)*(f(mod(i+2,n_max-1)+1) + f(i-3))         &
                + A(2)*(f(mod(i+1,n_max-1)+1) + f(i-2))         &
                + A(1)*(f(mod(i,n_max-1)+1) + f(i-1))           &
                + A(0)*f(i)
            enddo

        end if

        if (bc_type.eq.Dirichlet) then
            call D2c_CTR_CLOSURE_INF_ACC(f, d2f, h, 6)
        end if

        if (bc_type.eq.Dirichlet) then
            call D2c_CTR_CLOSURE_SUP_ACC(f, d2f, n_max, h, n_max-6, shifted)
        end if

        return

    end subroutine

    subroutine D2_ExpCtr_O0Fp6_hump_MULT(f, d2f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h,  geom_coefs(:)
        real(kind=8), intent(in)      ::f(:)
        real(kind=8), intent(out)   :: d2f(:)

        integer :: i
        real(kind=8) A(0:6)

        A(1) =  2.178831413991164d0     / h**2
        A(2) =  -0.62233954789826d0     / h**2
        A(3) =  0.27324369530023d0      / h**2
        A(4) =  -0.11423388546519d0     / h**2
        A(5) =  0.035021195634051d0     / h**2
        A(6) = -0.0054574850439034d0    / h**2
        A(0) = -2.d0*(A(1) + A(2) + A(3) + A(4) + A(5) + A(6))

        do i=7, n_max-7
            d2f(i)=   ( A(6)*(f(i+6) + f(i-6))    &
            + A(5)*(f(i+5) + f(i-5))    &
            + A(4)*(f(i+4) + f(i-4))    &
            + A(3)*(f(i+3) + f(i-3))    &
            + A(2)*(f(i+2) + f(i-2))    &
            + A(1)*(f(i+1) + f(i-1))    &
            + A(0)*f(i) ) * geom_coefs(i)
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 6
                d2f(i)=   ( A(6)*(f(i+6) + f(mod(n_max+i-8,n_max-1)+1))     &
                + A(5)*(f(i+5) + f(mod(n_max+i-7,n_max-1)+1))       &
                + A(4)*(f(i+4) + f(mod(n_max+i-6,n_max-1)+1))       &
                + A(3)*(f(i+3) + f(mod(n_max+i-5,n_max-1)+1))       &
                + A(2)*(f(i+2) + f(mod(n_max+i-4,n_max-1)+1))       &
                + A(1)*(f(i+1) + f(mod(n_max+i-3,n_max-1)+1))       &
                + A(0)*f(i) ) * geom_coefs(i)
            enddo

            do i=n_max-6, n_max-1
                d2f(i)=   ( A(6)*(f(mod(i+5,n_max-1)+1) + f(i-6))       &
                + A(5)*(f(mod(i+4,n_max-1)+1) + f(i-5))         &
                + A(4)*(f(mod(i+3,n_max-1)+1) + f(i-4))         &
                + A(3)*(f(mod(i+2,n_max-1)+1) + f(i-3))         &
                + A(2)*(f(mod(i+1,n_max-1)+1) + f(i-2))         &
                + A(1)*(f(mod(i,n_max-1)+1) + f(i-1))           &
                + A(0)*f(i) ) * geom_coefs(i)
            enddo

        end if

        if (bc_type.eq.Dirichlet) then
            call D2c_CTR_CLOSURE_INF_MULT(f, d2f, h, 6, geom_coefs)
        end if

        if (bc_type.eq.Dirichlet) then
            call D2c_CTR_CLOSURE_SUP_MULT(f, d2f, n_max, h, n_max-6, shifted, geom_coefs)
        end if

        return

    end subroutine

    subroutine D2_ExpCtr_O0Fp6_hump_MULT_ACC(f, d2f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h,  geom_coefs(:)
        real(kind=8), intent(in)      ::f(:)
        real(kind=8), intent(inout)   :: d2f(:)

        integer :: i
        real(kind=8) A(0:6)

        A(1) =  2.178831413991164d0     / h**2
        A(2) =  -0.62233954789826d0     / h**2
        A(3) =  0.27324369530023d0      / h**2
        A(4) =  -0.11423388546519d0     / h**2
        A(5) =  0.035021195634051d0     / h**2
        A(6) = -0.0054574850439034d0    / h**2
        A(0) = -2.d0*(A(1) + A(2) + A(3) + A(4) + A(5) + A(6))

        do i=7, n_max-7
            d2f(i)=   d2f(i) + ( A(6)*(f(i+6) + f(i-6))    &
            + A(5)*(f(i+5) + f(i-5))    &
            + A(4)*(f(i+4) + f(i-4))    &
            + A(3)*(f(i+3) + f(i-3))    &
            + A(2)*(f(i+2) + f(i-2))    &
            + A(1)*(f(i+1) + f(i-1))    &
            + A(0)*f(i) ) * geom_coefs(i)
        enddo

        if (bc_type.eq.periodic) then

            do i=1, 6
                d2f(i)=   d2f(i) + ( A(6)*(f(i+6) + f(mod(n_max+i-8,n_max-1)+1))     &
                + A(5)*(f(i+5) + f(mod(n_max+i-7,n_max-1)+1))       &
                + A(4)*(f(i+4) + f(mod(n_max+i-6,n_max-1)+1))       &
                + A(3)*(f(i+3) + f(mod(n_max+i-5,n_max-1)+1))       &
                + A(2)*(f(i+2) + f(mod(n_max+i-4,n_max-1)+1))       &
                + A(1)*(f(i+1) + f(mod(n_max+i-3,n_max-1)+1))       &
                + A(0)*f(i) ) * geom_coefs(i)
            enddo

            do i=n_max-6, n_max-1
                d2f(i)=   d2f(i) + ( A(6)*(f(mod(i+5,n_max-1)+1) + f(i-6))       &
                + A(5)*(f(mod(i+4,n_max-1)+1) + f(i-5))         &
                + A(4)*(f(mod(i+3,n_max-1)+1) + f(i-4))         &
                + A(3)*(f(mod(i+2,n_max-1)+1) + f(i-3))         &
                + A(2)*(f(mod(i+1,n_max-1)+1) + f(i-2))         &
                + A(1)*(f(mod(i,n_max-1)+1) + f(i-1))           &
                + A(0)*f(i) ) * geom_coefs(i)
            enddo

        end if

        if (bc_type.eq.Dirichlet) then
            call D2c_CTR_CLOSURE_INF_MULT_ACC(f, d2f, h, 6, geom_coefs)
        end if

        if (bc_type.eq.Dirichlet) then
            call D2c_CTR_CLOSURE_SUP_MULT_ACC(f, d2f, n_max, h, n_max-6, shifted, geom_coefs)
        end if

        return

    end subroutine

    subroutine D2_CptCtr_O6Fp0(f, d2f, n_max, h, shifted, bc_type)

        use Thomas_solver
        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d2f(:)

        integer :: i
        real(kind=8), dimension(:,:), save, allocatable :: A

        real(kind=8) rhs(n_max)

        real(kind=8), save    :: B(0:2),  B_lim(-1:2)

        integer, save :: last_nmax=0, last_bc_type1, last_bc_type2
        logical, save   :: last_shifted, update_matrix
        real(kind=8), save   :: last_h



        update_matrix=(last_nmax/=n_max).or.(last_bc_type1/=bc_type).or.(last_bc_type2/=bc_type)    &
        .or.(last_shifted.neqv.shifted).or.(last_h/=h)

        ! MATRIX A, B definition -----------------------------------------------------------
        !if (update_matrix) then
        call fill_matrix
        !end if

        last_nmax=n_max
        last_bc_type1=bc_type
        last_bc_type2=bc_type
        last_shifted=shifted
        last_h=h

        ! SCHEME CALCULATION ---------------------------------------------------------------
        do i=3, n_max-3
            rhs(i)= B(0)*f(i) + B(1)*( f(i+1) + f(i-1) ) + B(2)*( f(i+2) + f(i-2) )
        enddo

        if (bc_type.eq.periodic) then

            rhs(1)= B(0)*f(1) + B(1)*( f(2) + f(n_max-1) ) + B(2)*( f(3) + f(n_max-2) )
            rhs(2)= B(0)*f(2) + B(1)*( f(3) + f(1) ) + B(2)*( f(4) + f(n_max-1) )

            rhs(n_max-1)= B(0)*f(n_max-1) + B(1)*( f(1) + f(n_max-2) ) + B(2)*( f(2) + f(n_max-3) )
            rhs(n_max-2)= B(0)*f(n_max-2) + B(1)*( f(n_max-1) + f(n_max-3) ) + B(2)*( f(1) + f(n_max-4) )

            call TS_Pr(A(-1,1:n_max-1), A(0,1:n_max-1), A(1,1:n_max-1), rhs(1:n_max-1), d2f(1:n_max-1), n_max-1)

        end if

        if ((bc_type.eq.Dirichlet)) then

            rhs(2)= B_lim(-1)*f(1) + B_lim(0) * f(2) + B_lim(1) * f(3) + B_lim(2) * f(4)

            if (shifted) then

                rhs(n_max-2)= ( B_lim(-1)*f(n_max-1) + B_lim(0)*f(n_max-2) + B_lim(1)*f(n_max-3)+ B_lim(2)*f(n_max-4) )

                call TS_NPr(rhs(2:n_max-2), d2f(2:n_max-2), n_max-3)

            else

                rhs(n_max-2)= B(0)*f(n_max-2) + B(1)*( f(n_max-1) + f(n_max-3) ) + B(2)*( f(n_max) + f(n_max-4) )

                rhs(n_max-1)= B_lim(-1)*f(n_max) + B_lim(0)*f(n_max-1) + B_lim(1)*f(n_max-2)+ B_lim(2)*f(n_max-3)

                call TS_NPr(rhs(2:n_max-1), d2f(2:n_max-1), n_max-2)

            end if

        end if

        return

    contains


        subroutine fill_matrix()

            if (allocated(A)) then
                deallocate(A)
            end if

            allocate(A(-1:1, n_max))

            B(0)= -102.d0    / (44.d0*h**2)
            B(1)=12.d0      / (11.d0*h**2)
            B(2)= 3.d0      / (44.d0*h**2)

            do i=3, n_max-3
                A(-1,i) = 2.d0/11.d0
                A(0,i) = 1.d0
                A(1,i) = 2.d0/11.d0
            enddo

            select case (bc_type)

                case (periodic)

                    A(-1,1:2) = 2.d0/11.d0
                    A(0,1:2) = 1.d0
                    A(1,1:2) = 2.d0/11.d0

                    A(-1,n_max-2:n_max-1) = 2.d0/11.d0
                    A(0,n_max-2: n_max-1) = 1.d0
                    A(1,n_max-2:n_max-1) = 2.d0/11.d0

                case (Dirichlet)

                    A(-1,2)=0.d0
                    A(0,2)=1.d0
                    A(1,2)=-1.d0

                    B_lim(-1)   =   1.d0    / h**2
                    B_lim(0)    =   -3.d0   / h**2
                    B_lim(1)    =   3.d0    / h**2
                    B_lim(2)    =   -1.d0   / h**2

                    if (shifted) then

                        A(-1,n_max-2)=-1.d0
                        A(0,n_max-2)=1.d0
                        A(1,n_max-2)=0.d0

                        call TS_init(A(-1,2:n_max-2), A(0,2:n_max-2), A(1,2:n_max-2), n_max-3)

                    else

                        A(-1,n_max-2) = 2.d0/11.d0
                        A(0,n_max-2) = 1.d0
                        A(1,n_max-2) = 2.d0/11.d0

                        A(-1,n_max-1)=-1.d0
                        A(0,n_max-1)=1.d0
                        A(1,n_max-1)=0.d0

                        call TS_init(A(-1,2:n_max-1), A(0,2:n_max-1), A(1,2:n_max-1), n_max-2)

                    end if

                case default

            end select

        end subroutine fill_matrix


    end subroutine

    subroutine D2_CptCtr_O6Fp0_ACC(f, d2f, n_max, h, shifted, bc_type)

        implicit none
        logical, intent(in) :: shifted
        integer, intent(in) :: n_max, bc_type
        real(kind=8) , intent(in) :: h,  f(:)
        real(kind=8), intent(inout) :: d2f(:)

        real(kind=8), dimension(n_max) :: d2f_tmp
        integer :: i

        select case (bc_type)

            case (periodic)

                d2f_tmp(1:n_max-1)=d2f(1:n_max-1)
                call D2_CptCtr_O6Fp0(f, d2f_tmp, n_max, h, shifted, bc_type)
                d2f(1:n_max-1) = d2f(1:n_max-1) + d2f_tmp(1:n_max-1)

            case (Dirichlet)

                if (shifted) then

                    d2f_tmp(2:n_max-2)=d2f(2:n_max-2)
                    call D2_CptCtr_O6Fp0(f, d2f_tmp, n_max, h, shifted, bc_type)
                    d2f(2:n_max-2) = d2f(2:n_max-2) + d2f_tmp(2:n_max-2)

                else

                    d2f_tmp(2:n_max-1)=d2f(2:n_max-1)
                    call D2_CptCtr_O6Fp0(f, d2f_tmp, n_max, h, shifted, bc_type)
                    d2f(2:n_max-1) = d2f(2:n_max-1) + d2f_tmp(2:n_max-1)
                end if

            case default

        end select

        return

    end subroutine

    subroutine D2_CptCtr_O6Fp0_MULT(f, d2f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        logical, intent(in) :: shifted
        integer, intent(in) :: n_max, bc_type
        real(kind=8) , intent(in) :: h,  geom_coefs(:), f(:)
        real(kind=8), intent(out) :: d2f(:)

        integer :: i

        call D2_CptCtr_O6Fp0(f, d2f, n_max, h, shifted, bc_type)

        select case (bc_type)

            case (periodic)

                do i = 1, n_max-1
                    d2f(i)=d2f(i)*geom_coefs(i)
                end do

            case (Dirichlet)
                if (shifted) then
                    d2f(2:n_max-2)=d2f(2:n_max-2)*geom_coefs(2:n_max-2)
                else
                    d2f(2:n_max-1)=d2f(2:n_max-1)*geom_coefs(2:n_max-1)
                end if

            case default

        end select

        return

    end subroutine

    subroutine D2_CptCtr_O6Fp0_MULT_ACC(f, d2f, n_max, h, shifted, bc_type, geom_coefs)

        implicit none
        logical, intent(in) :: shifted
        integer, intent(in) :: n_max, bc_type
        real(kind=8) , intent(in) :: h,  geom_coefs(:), f(:)
        real(kind=8), intent(inout) :: d2f(:)

        real(kind=8), dimension(n_max) :: d2f_tmp
        integer :: i

        select case (bc_type)

            case (periodic)

                d2f_tmp(1:n_max-1)=d2f(1:n_max-1)
                call D2_CptCtr_O6Fp0(f, d2f_tmp, n_max, h, shifted, bc_type)
                d2f(1:n_max-1) = d2f(1:n_max-1) + d2f_tmp(1:n_max-1)*geom_coefs(1:n_max-1)

            case (Dirichlet)

                if (shifted) then

                    d2f_tmp(2:n_max-2)=d2f(2:n_max-2)
                    call D2_CptCtr_O6Fp0(f, d2f_tmp, n_max, h, shifted, bc_type)
                    d2f(2:n_max-2) = d2f(2:n_max-2) + d2f_tmp(2:n_max-2)*geom_coefs(2:n_max-2)

                else

                    d2f_tmp(2:n_max-1)=d2f(2:n_max-1)
                    call D2_CptCtr_O6Fp0(f, d2f_tmp, n_max, h, shifted, bc_type)
                    d2f(2:n_max-1) = d2f(2:n_max-1) + d2f_tmp(2:n_max-1)*geom_coefs(2:n_max-1)
                end if

            case default

        end select

        return

    end subroutine

    subroutine D2c_CTR_CLOSURE_INF(f, d2f, h, imax)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(imax)

        implicit none

        integer :: i, n_max, imax
        real(kind=8) h
        real(kind=8) d2f(:)
        real(kind=8) f(:)

        real(kind=8) d2f_temp(CPT_MIN)


        if (treat_boundary_by_cpt_scheme) then

            call D2_CptCtr_O6Fp0(  &
            f(1:CPT_MIN), d2f_temp(1:CPT_MIN), CPT_MIN, h, .true., Dirichlet)

            d2f(2)=d2f_temp(2)
            if (imax.eq.2) return

            d2f(3)=d2f_temp(3)
            if (imax.eq.3) return

        else

            ! i=2---------------------------------------------------
            i=2
            d2f(i)  =   ( -2.0d0*f(i) + (f(i+1)+f(i-1)) ) / h**2

            if (i.eq.imax) return
            ! i=3---------------------------------------------------
            i=3
            d2f(i)  =   -5.0d0  * f(i)  / (2.d0*h**2)   &
            +   4.d0  * (f(i+1)+f(i-1)) / (3.d0*h**2)   &
            -   1.d0 * (f(i+2) + f(i-2))    / (12.d0*h**2)

            if (i.eq.imax) return

        endif

        ! i=4---------------------------------------------------
        i=4
        d2f(i)  =   - 2.814728882213931d0 * f(i)    / h**2      &
        +   1.569379994993781d0*(f(i+1)+f(i-1))     / h**2      &
        -   0.177751998d0*(f(i+2)+f(i-2))           / h**2      &
        +   0.0157364441d0* (f(i+3) + f(i-3))           / h**2

        if (i.eq.imax) return

        ! i=5---------------------------------------------------
        i=5
        d2f(i)  =   - 2.814728882213931d0 * f(i)    / h**2      &
        +   1.569379994993781d0*(f(i+1)+f(i-1))     / h**2      &
        -   0.177751998d0*(f(i+2)+f(i-2))           / h**2      &
        +   0.0157364441d0* (f(i+3) + f(i-3))           / h**2

        if (i.eq.imax) return

        ! i=6---SEE D1_EXP_5pts.wxm-----------------------------
        i=6
        d2f(i)  =   (                                               &
        -   3.07898593599001d0      * f(i)                  &
        +   1.797719245977798d0     * (f(i+1)+f(i-1))       &
        -   0.32149582867581d0      * (f(i+2)+f(i-2))       &
        +   0.07716473796078d0      * (f(i+3) + f(i-3))     &
        -   0.015684734192127d0     * (f(i+4) + f(i-4))     &
        +   0.0017895469243599d0    * (f(i+5) + f(i-5))     &
        ) /h**2

        if (i.eq.imax) return

        ! i=7---SEE D1_EXP_6pts.wxm-----------------------------
        i=7
        d2f(i)  =   (                                       &
        -   3.177647869182281d0     * f(i)                  &
        +   1.889601043267912d0     * (f(i+1)+f(i-1))       &
        -   0.39502242251325d0      * (f(i+2)+f(i-2))       &
        +   0.12612555410575d0      * (f(i+3) + f(i-3))     &
        -   0.04095084311657d0      * (f(i+4) + f(i-4))     &
        +   0.0105430010251d0       * (f(i+5) + f(i-5))     &
        -   0.0014723981778037d0    * (f(i+6) + f(i-6))     &
        ) / h**2


        return

    end subroutine

    subroutine D2c_CTR_CLOSURE_INF_ACC(f, d2f, h, imax)

        implicit none

        integer :: i, n_max, imax
        real(kind=8) h
        real(kind=8) d2f(:)
        real(kind=8) f(:)

        real(kind=8) d2f_temp(CPT_MIN)


        if (treat_boundary_by_cpt_scheme) then

            d2f_temp(1:CPT_MIN)=d2f(1:CPT_MIN)
            call D2_CptCtr_O6Fp0_ACC(  &
            f(1:CPT_MIN), d2f_temp(1:CPT_MIN), CPT_MIN, h, .true., Dirichlet)

            d2f(2)=d2f_temp(2)
            if (imax.eq.2) return

            d2f(3)=d2f_temp(3)
            if (imax.eq.3) return

        else

            ! i=2---------------------------------------------------
            i=2
            d2f(i)  =  d2f(i) + ( -2.0d0*f(i) + (f(i+1)+f(i-1)) ) / h**2

            if (i.eq.imax) return
            ! i=3---------------------------------------------------
            i=3
            d2f(i)  =   d2f(i) + (-5.0d0  * f(i)  / (2.d0*h**2)   &
            +   4.d0  * (f(i+1)+f(i-1)) / (3.d0*h**2)   &
            -   1.d0 * (f(i+2) + f(i-2))    / (12.d0*h**2) )

            if (i.eq.imax) return

        endif

        ! i=4---------------------------------------------------
        i=4
        d2f(i)  =   d2f(i) + (- 2.814728882213931d0 * f(i)    / h**2      &
        +   1.569379994993781d0*(f(i+1)+f(i-1))     / h**2      &
        -   0.177751998d0*(f(i+2)+f(i-2))           / h**2      &
        +   0.0157364441d0* (f(i+3) + f(i-3))           / h**2)

        if (i.eq.imax) return

        ! i=5---------------------------------------------------
        i=5
        d2f(i)  =   d2f(i) + (- 2.814728882213931d0 * f(i)    / h**2      &
        +   1.569379994993781d0*(f(i+1)+f(i-1))     / h**2      &
        -   0.177751998d0*(f(i+2)+f(i-2))           / h**2      &
        +   0.0157364441d0* (f(i+3) + f(i-3))           / h**2)

        if (i.eq.imax) return

        ! i=6---SEE D1_EXP_5pts.wxm-----------------------------
        i=6
        d2f(i)  =   d2f(i) + (                                               &
        -   3.07898593599001d0      * f(i)                  &
        +   1.797719245977798d0     * (f(i+1)+f(i-1))       &
        -   0.32149582867581d0      * (f(i+2)+f(i-2))       &
        +   0.07716473796078d0      * (f(i+3) + f(i-3))     &
        -   0.015684734192127d0     * (f(i+4) + f(i-4))     &
        +   0.0017895469243599d0    * (f(i+5) + f(i-5))     &
        ) /h**2

        if (i.eq.imax) return

        ! i=7---SEE D1_EXP_6pts.wxm-----------------------------
        i=7
        d2f(i)  =   d2f(i) + (                                       &
        -   3.177647869182281d0     * f(i)                  &
        +   1.889601043267912d0     * (f(i+1)+f(i-1))       &
        -   0.39502242251325d0      * (f(i+2)+f(i-2))       &
        +   0.12612555410575d0      * (f(i+3) + f(i-3))     &
        -   0.04095084311657d0      * (f(i+4) + f(i-4))     &
        +   0.0105430010251d0       * (f(i+5) + f(i-5))     &
        -   0.0014723981778037d0    * (f(i+6) + f(i-6))     &
        ) / h**2


        return

    end subroutine

    subroutine D2c_CTR_CLOSURE_INF_MULT(f, d2f, h, imax, geom_coefs)

        implicit none

        integer :: i, n_max, imax
        real(kind=8) h
        real(kind=8) d2f(:), geom_coefs(:)
        real(kind=8) f(:)

        real(kind=8) d2f_temp(CPT_MIN)


        if (treat_boundary_by_cpt_scheme) then

            call D2_CptCtr_O6Fp0_MULT(  &
            f(1:CPT_MIN), d2f_temp(1:CPT_MIN), CPT_MIN, h, .true., Dirichlet, geom_coefs)

            d2f(2)=d2f_temp(2)
            if (imax.eq.2) return

            d2f(3)=d2f_temp(3)
            if (imax.eq.3) return

        else

            ! i=2---------------------------------------------------
            i=2
            d2f(i)  =  ( ( -2.0d0*f(i) + (f(i+1)+f(i-1)) ) / h**2) * geom_coefs(i)

            if (i.eq.imax) return
            ! i=3---------------------------------------------------
            i=3
            d2f(i)  =   (-5.0d0  * f(i)  / (2.d0*h**2)   &
            +   4.d0  * (f(i+1)+f(i-1)) / (3.d0*h**2)   &
            -   1.d0 * (f(i+2) + f(i-2))    / (12.d0*h**2) ) * geom_coefs(i)

            if (i.eq.imax) return

        endif

        ! i=4---------------------------------------------------
        i=4
        d2f(i)  =   (- 2.814728882213931d0 * f(i)    / h**2      &
        +   1.569379994993781d0*(f(i+1)+f(i-1))     / h**2      &
        -   0.177751998d0*(f(i+2)+f(i-2))           / h**2      &
        +   0.0157364441d0* (f(i+3) + f(i-3))           / h**2) * geom_coefs(i)

        if (i.eq.imax) return

        ! i=5---------------------------------------------------
        i=5
        d2f(i)  =   (- 2.814728882213931d0 * f(i)    / h**2      &
        +   1.569379994993781d0*(f(i+1)+f(i-1))     / h**2      &
        -   0.177751998d0*(f(i+2)+f(i-2))           / h**2      &
        +   0.0157364441d0* (f(i+3) + f(i-3))           / h**2) * geom_coefs(i)

        if (i.eq.imax) return

        ! i=6---SEE D1_EXP_5pts.wxm-----------------------------
        i=6
        d2f(i)  =   (                                               &
        -   3.07898593599001d0      * f(i)                  &
        +   1.797719245977798d0     * (f(i+1)+f(i-1))       &
        -   0.32149582867581d0      * (f(i+2)+f(i-2))       &
        +   0.07716473796078d0      * (f(i+3) + f(i-3))     &
        -   0.015684734192127d0     * (f(i+4) + f(i-4))     &
        +   0.0017895469243599d0    * (f(i+5) + f(i-5))     &
        )  * geom_coefs(i) /h**2

        if (i.eq.imax) return

        ! i=7---SEE D1_EXP_6pts.wxm-----------------------------
        i=7
        d2f(i)  =   (                                       &
        -   3.177647869182281d0     * f(i)                  &
        +   1.889601043267912d0     * (f(i+1)+f(i-1))       &
        -   0.39502242251325d0      * (f(i+2)+f(i-2))       &
        +   0.12612555410575d0      * (f(i+3) + f(i-3))     &
        -   0.04095084311657d0      * (f(i+4) + f(i-4))     &
        +   0.0105430010251d0       * (f(i+5) + f(i-5))     &
        -   0.0014723981778037d0    * (f(i+6) + f(i-6))     &
        )  * geom_coefs(i) / h**2


        return

    end subroutine

    subroutine D2c_CTR_CLOSURE_INF_MULT_ACC(f, d2f, h, imax, geom_coefs)

        implicit none

        integer :: i, n_max, imax
        real(kind=8) h
        real(kind=8) d2f(:), geom_coefs(:)
        real(kind=8) f(:)

        real(kind=8) d2f_temp(CPT_MIN)


        if (treat_boundary_by_cpt_scheme) then

            d2f_temp(1:CPT_MIN)=d2f(1:CPT_MIN)
            call D2_CptCtr_O6Fp0_MULT_ACC(  &
            f(1:CPT_MIN), d2f_temp(1:CPT_MIN), CPT_MIN, h, .true., Dirichlet, geom_coefs)

            d2f(2)=d2f_temp(2)
            if (imax.eq.2) return

            d2f(3)=d2f_temp(3)
            if (imax.eq.3) return

        else

            ! i=2---------------------------------------------------
            i=2
            d2f(i)  =  d2f(i) + ( ( -2.0d0*f(i) + (f(i+1)+f(i-1)) ) / h**2) * geom_coefs(i)

            if (i.eq.imax) return
            ! i=3---------------------------------------------------
            i=3
            d2f(i)  =   d2f(i) + (-5.0d0  * f(i)  / (2.d0*h**2)   &
            +   4.d0  * (f(i+1)+f(i-1)) / (3.d0*h**2)   &
            -   1.d0 * (f(i+2) + f(i-2))    / (12.d0*h**2) ) * geom_coefs(i)

            if (i.eq.imax) return

        endif

        ! i=4---------------------------------------------------
        i=4
        d2f(i)  =   d2f(i) + (- 2.814728882213931d0 * f(i)    / h**2      &
        +   1.569379994993781d0*(f(i+1)+f(i-1))     / h**2      &
        -   0.177751998d0*(f(i+2)+f(i-2))           / h**2      &
        +   0.0157364441d0* (f(i+3) + f(i-3))           / h**2) * geom_coefs(i)

        if (i.eq.imax) return

        ! i=5---------------------------------------------------
        i=5
        d2f(i)  =   d2f(i) + (- 2.814728882213931d0 * f(i)    / h**2      &
        +   1.569379994993781d0*(f(i+1)+f(i-1))     / h**2      &
        -   0.177751998d0*(f(i+2)+f(i-2))           / h**2      &
        +   0.0157364441d0* (f(i+3) + f(i-3))           / h**2) * geom_coefs(i)

        if (i.eq.imax) return

        ! i=6---SEE D1_EXP_5pts.wxm-----------------------------
        i=6
        d2f(i)  =   d2f(i) + (                                               &
        -   3.07898593599001d0      * f(i)                  &
        +   1.797719245977798d0     * (f(i+1)+f(i-1))       &
        -   0.32149582867581d0      * (f(i+2)+f(i-2))       &
        +   0.07716473796078d0      * (f(i+3) + f(i-3))     &
        -   0.015684734192127d0     * (f(i+4) + f(i-4))     &
        +   0.0017895469243599d0    * (f(i+5) + f(i-5))     &
        )  * geom_coefs(i) /h**2

        if (i.eq.imax) return

        ! i=7---SEE D1_EXP_6pts.wxm-----------------------------
        i=7
        d2f(i)  =   d2f(i) + (                                       &
        -   3.177647869182281d0     * f(i)                  &
        +   1.889601043267912d0     * (f(i+1)+f(i-1))       &
        -   0.39502242251325d0      * (f(i+2)+f(i-2))       &
        +   0.12612555410575d0      * (f(i+3) + f(i-3))     &
        -   0.04095084311657d0      * (f(i+4) + f(i-4))     &
        +   0.0105430010251d0       * (f(i+5) + f(i-5))     &
        -   0.0014723981778037d0    * (f(i+6) + f(i-6))     &
        )  * geom_coefs(i) / h**2


        return

    end subroutine

    subroutine D2c_CTR_CLOSURE_SUP(f, d2f, n_max, h, imin, shifted)

        implicit none

        integer :: i, n_max, imin, j
        real(kind=8) h, fn
        real(kind=8) d2f(:)
        real(kind=8) f(:)
        logical :: shifted

        real(kind=8) d2f_temp(n_max-CPT_MIN:n_max-1)

        if (shifted) then
            i=n_max-2
        else
            i=n_max-1
        end if

        if (treat_boundary_by_cpt_scheme) then
            call D2_CptCtr_O6Fp0(f(n_max-CPT_MIN:n_max), d2f_temp, CPT_MIN+1, h, shifted, Dirichlet)

            d2f(i)=d2f_temp(i)
            if (i.eq.imin) return

            i=i-1
            d2f(i)=d2f_temp(i)
            if (i.eq.imin) return

        else

            ! First point-------------------------------------------
            d2f(i)  =   ( -2.0d0*f(i) + (f(i+1)+f(i-1)) ) / h**2

            if (i.eq.imin) return
            ! 2nd point----------------------------------------------
            i=i-1
            d2f(i)  =   -5.0d0  * f(i)  / (2.d0*h**2)   &
            +   4.d0  * (f(i+1)+f(i-1)) / (3.d0*h**2)   &
            -   1.d0 * (f(i+2) + f(i-2))    / (12.d0*h**2)

            if (i.eq.imin) return

        endif

        ! 3rd point---------------------------------------------
        i=i-1
        d2f(i)  =   - 2.814728882213931d0 * f(i)    / h**2      &
        +   1.569379994993781d0*(f(i+1)+f(i-1))     / h**2      &
        -   0.177751998d0*(f(i+2)+f(i-2))           / h**2      &
        +   0.0157364441d0* (f(i+3) + f(i-3))           / h**2

        if (i.eq.imin) return

        ! 4th point---------------------------------------------
        i=i-1
        d2f(i)  =   - 2.814728882213931d0 * f(i)    / h**2      &
        +   1.569379994993781d0*(f(i+1)+f(i-1))     / h**2      &
        -   0.177751998d0*(f(i+2)+f(i-2))           / h**2      &
        +   0.0157364441d0* (f(i+3) + f(i-3))           / h**2

        if (i.eq.imin) return

        ! 5th point --SEE D2_EXP_6pts.wxm---wc=2.5---hump=1.0---
        i=i-1
        d2f(i)  =   (                                               &
        -   3.07898593599001d0      * f(i)                  &
        +   1.797719245977798d0     * (f(i+1)+f(i-1))       &
        -   0.32149582867581d0      * (f(i+2)+f(i-2))       &
        +   0.07716473796078d0      * (f(i+3) + f(i-3))     &
        -   0.015684734192127d0     * (f(i+4) + f(i-4))     &
        +   0.0017895469243599d0    * (f(i+5) + f(i-5))     &
        ) /h**2

        if (i.eq.imin) return

        ! 6th point --SEE D2_EXP_6pts.wxm---wc=2.5---hump=1.1---
        i=i-1

        d2f(i)  =   (                                           &
        -   3.4901307730361837d0     * f(i)                     &
        +   2.178831413991164d0     * (f(i+1)+f(i-1))           &
        -   0.62233954789826d0      * (f(i+2)+f(i-2))           &
        +   0.27324369530023d0      * (f(i+3) + f(i-3))         &
        -   0.11423388546519d0      * (f(i+4) + f(i-4))         &
        +   0.035021195634051d0       * (f(i+5) + f(i-5))       &
        -   0.0054574850439034d0    * (f(i+6) + f(i-6))         &
        ) / h**2

        if (i.eq.imin) return


        return

    end subroutine

    subroutine D2c_CTR_CLOSURE_SUP_ACC(f, d2f, n_max, h, imin, shifted)

        implicit none

        integer :: i, n_max, imin, j
        real(kind=8) h, fn
        real(kind=8) d2f(:)
        real(kind=8) f(:)
        logical :: shifted

        real(kind=8) d2f_temp(n_max-CPT_MIN:n_max-1)

        if (shifted) then
            i=n_max-2
        else
            i=n_max-1
        end if

        if (treat_boundary_by_cpt_scheme) then
            d2f_temp=d2f(n_max-CPT_MIN:n_max-1)
            call D2_CptCtr_O6Fp0_ACC(f(n_max-CPT_MIN:n_max), d2f_temp, CPT_MIN+1, h, shifted, Dirichlet)

            d2f(i)=d2f_temp(i)
            if (i.eq.imin) return

            i=i-1
            d2f(i)=d2f_temp(i)
            if (i.eq.imin) return

        else

            ! First point-------------------------------------------
            d2f(i)  =   d2f(i) + ( -2.0d0*f(i) + (f(i+1)+f(i-1)) ) / h**2

            if (i.eq.imin) return
            ! 2nd point----------------------------------------------
            i=i-1
            d2f(i)  =   d2f(i) + ( -5.0d0  * f(i)  / (2.d0*h**2)   &
            +   4.d0  * (f(i+1)+f(i-1)) / (3.d0*h**2)   &
            -   1.d0 * (f(i+2) + f(i-2))    / (12.d0*h**2) )

            if (i.eq.imin) return

        endif

        ! 3rd point---------------------------------------------
        i=i-1
        d2f(i)  =   d2f(i) + (- 2.814728882213931d0 * f(i)    / h**2      &
        +   1.569379994993781d0*(f(i+1)+f(i-1))     / h**2      &
        -   0.177751998d0*(f(i+2)+f(i-2))           / h**2      &
        +   0.0157364441d0* (f(i+3) + f(i-3))           / h**2 )

        if (i.eq.imin) return

        ! 4th point---------------------------------------------
        i=i-1
        d2f(i)  =   d2f(i) + ( - 2.814728882213931d0 * f(i)    / h**2      &
        +   1.569379994993781d0*(f(i+1)+f(i-1))     / h**2      &
        -   0.177751998d0*(f(i+2)+f(i-2))           / h**2      &
        +   0.0157364441d0* (f(i+3) + f(i-3))           / h**2 )

        if (i.eq.imin) return

        ! 5th point --SEE D2_EXP_6pts.wxm---wc=2.5---hump=1.0---
        i=i-1
        d2f(i)  =   d2f(i) + (                                               &
        -   3.07898593599001d0      * f(i)                  &
        +   1.797719245977798d0     * (f(i+1)+f(i-1))       &
        -   0.32149582867581d0      * (f(i+2)+f(i-2))       &
        +   0.07716473796078d0      * (f(i+3) + f(i-3))     &
        -   0.015684734192127d0     * (f(i+4) + f(i-4))     &
        +   0.0017895469243599d0    * (f(i+5) + f(i-5))     &
        ) /h**2

        if (i.eq.imin) return

        ! 6th point --SEE D2_EXP_6pts.wxm---wc=2.5---hump=1.1---
        i=i-1

        d2f(i)  =   d2f(i) + (                                           &
        -   3.4901307730361837d0     * f(i)                     &
        +   2.178831413991164d0     * (f(i+1)+f(i-1))           &
        -   0.62233954789826d0      * (f(i+2)+f(i-2))           &
        +   0.27324369530023d0      * (f(i+3) + f(i-3))         &
        -   0.11423388546519d0      * (f(i+4) + f(i-4))         &
        +   0.035021195634051d0       * (f(i+5) + f(i-5))       &
        -   0.0054574850439034d0    * (f(i+6) + f(i-6))         &
        ) / h**2

        if (i.eq.imin) return


        return

    end subroutine

    subroutine D2c_CTR_CLOSURE_SUP_MULT(f, d2f, n_max, h, imin, shifted, geom_coefs)

        implicit none

        integer :: i, n_max, imin, j
        real(kind=8) h, fn
        real(kind=8) d2f(:), geom_coefs(:)
        real(kind=8) f(:)
        logical :: shifted

        real(kind=8) d2f_temp(n_max-CPT_MIN:n_max-1)

        if (shifted) then
            i=n_max-2
        else
            i=n_max-1
        end if

        if (treat_boundary_by_cpt_scheme) then
            call D2_CptCtr_O6Fp0_MULT(f(n_max-CPT_MIN:n_max), d2f_temp, CPT_MIN+1, h, shifted,     &
            Dirichlet, geom_coefs(n_max-CPT_MIN:n_max))

            d2f(i)=d2f_temp(i)
            if (i.eq.imin) return

            i=i-1
            d2f(i)=d2f_temp(i)
            if (i.eq.imin) return

        else

            ! First point-------------------------------------------
            d2f(i)  =   ( -2.0d0*f(i) + (f(i+1)+f(i-1)) ) * geom_coefs(i) / h**2

            if (i.eq.imin) return
            ! 2nd point----------------------------------------------
            i=i-1
            d2f(i)  =   ( -5.0d0  * f(i)  / (2.d0*h**2)   &
            +   4.d0  * (f(i+1)+f(i-1)) / (3.d0*h**2)   &
            -   1.d0 * (f(i+2) + f(i-2))    / (12.d0*h**2) ) * geom_coefs(i)

            if (i.eq.imin) return

        endif

        ! 3rd point---------------------------------------------
        i=i-1
        d2f(i)  =   (- 2.814728882213931d0 * f(i)    / h**2      &
        +   1.569379994993781d0*(f(i+1)+f(i-1))     / h**2      &
        -   0.177751998d0*(f(i+2)+f(i-2))           / h**2      &
        +   0.0157364441d0* (f(i+3) + f(i-3))           / h**2 ) * geom_coefs(i)

        if (i.eq.imin) return

        ! 4th point---------------------------------------------
        i=i-1
        d2f(i)  =   ( - 2.814728882213931d0 * f(i)    / h**2      &
        +   1.569379994993781d0*(f(i+1)+f(i-1))     / h**2      &
        -   0.177751998d0*(f(i+2)+f(i-2))           / h**2      &
        +   0.0157364441d0* (f(i+3) + f(i-3))           / h**2 ) * geom_coefs(i)

        if (i.eq.imin) return

        ! 5th point --SEE D2_EXP_6pts.wxm---wc=2.5---hump=1.0---
        i=i-1
        d2f(i)  =   (                                               &
        -   3.07898593599001d0      * f(i)                  &
        +   1.797719245977798d0     * (f(i+1)+f(i-1))       &
        -   0.32149582867581d0      * (f(i+2)+f(i-2))       &
        +   0.07716473796078d0      * (f(i+3) + f(i-3))     &
        -   0.015684734192127d0     * (f(i+4) + f(i-4))     &
        +   0.0017895469243599d0    * (f(i+5) + f(i-5))     &
        )  * geom_coefs(i) /h**2

        if (i.eq.imin) return

        ! 6th point --SEE D2_EXP_6pts.wxm---wc=2.5---hump=1.1---
        i=i-1

        d2f(i)  =   (                                           &
        -   3.4901307730361837d0     * f(i)                     &
        +   2.178831413991164d0     * (f(i+1)+f(i-1))           &
        -   0.62233954789826d0      * (f(i+2)+f(i-2))           &
        +   0.27324369530023d0      * (f(i+3) + f(i-3))         &
        -   0.11423388546519d0      * (f(i+4) + f(i-4))         &
        +   0.035021195634051d0       * (f(i+5) + f(i-5))       &
        -   0.0054574850439034d0    * (f(i+6) + f(i-6))         &
        )  * geom_coefs(i) / h**2

        if (i.eq.imin) return

        return

    end subroutine

    subroutine D2c_CTR_CLOSURE_SUP_MULT_ACC(f, d2f, n_max, h, imin, shifted, geom_coefs)

        implicit none

        integer :: i, n_max, imin, j
        real(kind=8) h, fn
        real(kind=8) d2f(:), geom_coefs(:)
        real(kind=8) f(:)
        logical :: shifted

        real(kind=8) d2f_temp(n_max-CPT_MIN:n_max-1)

        if (shifted) then
            i=n_max-2
        else
            i=n_max-1
        end if

        if (treat_boundary_by_cpt_scheme) then

            d2f_temp=d2f(n_max-CPT_MIN:n_max-1)
            call D2_CptCtr_O6Fp0_MULT_ACC(f(n_max-CPT_MIN:n_max), d2f_temp, CPT_MIN+1, h, shifted,     &
            Dirichlet, geom_coefs(n_max-CPT_MIN:n_max))

            d2f(i)=d2f_temp(i)
            if (i.eq.imin) return

            i=i-1
            d2f(i)=d2f_temp(i)
            if (i.eq.imin) return

        else

            ! First point-------------------------------------------
            d2f(i)  =   d2f(i) + ( -2.0d0*f(i) + (f(i+1)+f(i-1)) ) * geom_coefs(i) / h**2

            if (i.eq.imin) return
            ! 2nd point----------------------------------------------
            i=i-1
            d2f(i)  =   d2f(i) + ( -5.0d0  * f(i)  / (2.d0*h**2)   &
            +   4.d0  * (f(i+1)+f(i-1)) / (3.d0*h**2)   &
            -   1.d0 * (f(i+2) + f(i-2))    / (12.d0*h**2) ) * geom_coefs(i)

            if (i.eq.imin) return

        endif

        ! 3rd point---------------------------------------------
        i=i-1
        d2f(i)  =   d2f(i) + (- 2.814728882213931d0 * f(i)    / h**2      &
        +   1.569379994993781d0*(f(i+1)+f(i-1))     / h**2      &
        -   0.177751998d0*(f(i+2)+f(i-2))           / h**2      &
        +   0.0157364441d0* (f(i+3) + f(i-3))           / h**2 ) * geom_coefs(i)

        if (i.eq.imin) return

        ! 4th point---------------------------------------------
        i=i-1
        d2f(i)  =   d2f(i) + ( - 2.814728882213931d0 * f(i)    / h**2      &
        +   1.569379994993781d0*(f(i+1)+f(i-1))     / h**2      &
        -   0.177751998d0*(f(i+2)+f(i-2))           / h**2      &
        +   0.0157364441d0* (f(i+3) + f(i-3))           / h**2 ) * geom_coefs(i)

        if (i.eq.imin) return

        ! 5th point --SEE D2_EXP_6pts.wxm---wc=2.5---hump=1.0---
        i=i-1
        d2f(i)  =   d2f(i) + (                                               &
        -   3.07898593599001d0      * f(i)                  &
        +   1.797719245977798d0     * (f(i+1)+f(i-1))       &
        -   0.32149582867581d0      * (f(i+2)+f(i-2))       &
        +   0.07716473796078d0      * (f(i+3) + f(i-3))     &
        -   0.015684734192127d0     * (f(i+4) + f(i-4))     &
        +   0.0017895469243599d0    * (f(i+5) + f(i-5))     &
        )  * geom_coefs(i) /h**2

        if (i.eq.imin) return

        ! 6th point --SEE D2_EXP_6pts.wxm---wc=2.5---hump=1.1---
        i=i-1

        d2f(i)  =   d2f(i) + (                                           &
        -   3.4901307730361837d0     * f(i)                     &
        +   2.178831413991164d0     * (f(i+1)+f(i-1))           &
        -   0.62233954789826d0      * (f(i+2)+f(i-2))           &
        +   0.27324369530023d0      * (f(i+3) + f(i-3))         &
        -   0.11423388546519d0      * (f(i+4) + f(i-4))         &
        +   0.035021195634051d0       * (f(i+5) + f(i-5))       &
        -   0.0054574850439034d0    * (f(i+6) + f(i-6))         &
        )  * geom_coefs(i) / h**2

        return

    end subroutine

end module d2c_schemes

! This module is used for Poisson solving
module d2c_from_d1s
    use boundaries_types
    use d0s_schemes
    use d1s_schemes_ibm

contains

    ! only for periodic case
    subroutine apply_d1s_twice(d1s_scheme, f, d2f, n_max, h, shifted)

        implicit none
        integer :: i, n_max
        procedure(D0s_ExpCtr_O2Fp0) :: d1s_scheme
        logical :: shifted, d1f_shifted
        real(kind=8) SoV(5)

        real(kind=8) h
        real(kind=8) d2f(n_max), d2f0(n_max+1), d1f(n_max)
        real(kind=8) f(n_max)


        call d1s_scheme(f, d1f, n_max, h, shifted, periodic)

        if (shifted) then
            d1f_shifted=.false.
        else
            d1f_shifted=.true.
        end if

        call d1s_scheme(d1f, d2f, n_max, h, d1f_shifted, periodic)

        return

    end subroutine

    real(kind=8) function K2eq_of_d1s_twice(d1s_scheme, N, k0, L, shifted)

        implicit none
        real(kind=8) :: dx, PI, L, scale_factor
        integer :: i, N, k, k0
        procedure(D0s_ExpCtr_O2Fp0) :: d1s_scheme
        logical :: shifted
        real(kind=8)  x(N),f1(N), f2(N), df1_calc(N), df2_calc(N)

        PI=3.141592653589793238462643d0
        k=k0-1

        dx=(2.d0*PI/(N-1))
        do i=1, N
            x(i)=1.d0+(i-1)*dx
        enddo

        do i=1, N
            f1(i)=dcos(k*x(i))
            f2(i)=dsin(k*x(i))
        enddo

        call apply_d1s_twice(d1s_scheme, f1, df1_calc, N, dx, shifted)
        call apply_d1s_twice(d1s_scheme, f2, df2_calc, N, dx, shifted)

        i=N/2
        scale_factor=(2.d0/L)**2
        K2eq_of_d1s_twice =  -(  f1(i)*df1_calc(i) + f2(i)*df2_calc(i) )*scale_factor


        return
    end function

    ! only for periodic case
    subroutine apply_d1s_twice_ibm(d1s_scheme, f, d2f, n_max, h, shifted, n_re_fc, n_objects_fc, n_re_cc, n_objects_cc)

        implicit none
        integer :: i, n_max
        procedure(D1s_ExpCtr_O2Fp0_ibm) :: d1s_scheme
        logical :: shifted, d1f_shifted
        real(kind=8) SoV(5)

        real(kind=8) h
        real(kind=8) d2f(n_max), d2f0(n_max+1), d1f(n_max)
        real(kind=8) f(n_max)
        integer :: n_objects_fc, n_objects_cc
        integer, dimension(:) :: n_re_fc, n_re_cc

        call d1s_scheme(f, d1f, n_max, h, shifted, n_re_fc, n_objects_fc)

        if (shifted) then
            d1f_shifted=.false.
        else
            d1f_shifted=.true.
        end if

        call d1s_scheme(d1f, d2f, n_max, h, d1f_shifted, n_re_cc, n_objects_cc)

        return

    end subroutine

    real(kind=8) function K2eq_of_d1s_twice_ibm(d1s_scheme, N, k0, L, shifted, n_re_fc, n_objects_fc, n_re_cc, n_objects_cc)

        implicit none
        real(kind=8) :: dx, PI, L, scale_factor
        integer :: i, N, k, k0
        procedure(D1s_ExpCtr_O2Fp0_ibm) :: d1s_scheme
        logical :: shifted
        real(kind=8)  x(N),f1(N), f2(N), df1_calc(N), df2_calc(N)
        integer :: n_objects_fc, n_objects_cc
        integer, dimension(:) :: n_re_fc, n_re_cc

        PI=3.141592653589793238462643d0
        k=k0-1

        dx=(2.d0*PI/(N-1))
        do i=1, N
            x(i)=1.d0+(i-1)*dx
        enddo

        do i=1, N
            f1(i)=dcos(k*x(i))
            f2(i)=dsin(k*x(i))
        enddo

        call apply_d1s_twice_ibm(d1s_scheme, f1, df1_calc, N, dx, shifted, n_re_fc, n_objects_fc, n_re_cc, n_objects_cc)
        call apply_d1s_twice_ibm(d1s_scheme, f2, df2_calc, N, dx, shifted, n_re_fc, n_objects_fc, n_re_cc, n_objects_cc)

        i=N/2
        scale_factor=(2.d0/L)**2
        K2eq_of_d1s_twice_ibm =  -(  f1(i)*df1_calc(i) + f2(i)*df2_calc(i) )*scale_factor


        return
    end function

    subroutine set_d2c_coeffs(d1sA_coeffs, d1sB_coeffs, d2c_coeffs, A_size, B_size, geom_A_coeffs, geom_B_coeff)

        implicit none
        integer A_size, B_size, i, j
        real(kind=8), dimension(-A_size:A_size)   :: d1sA_coeffs
        real(kind=8), dimension(-B_size:B_size)   :: geom_A_coeffs
        real(kind=8), dimension(-B_size:B_size)   :: d1sB_coeffs
        real(kind=8), dimension(-(A_size+B_size)+1:(A_size+B_size)-1) :: d2c_coeffs
        real(kind=8) ::geom_B_coeff

        d2c_coeffs=0.d0

        do i = -B_size, -1
            do j = -A_size, -1
                d2c_coeffs(i+j+1)= d2c_coeffs(i+j+1)+d1sB_coeffs(i)*d1sA_coeffs(j) * geom_A_coeffs(i)
            end do
            do j = 1, A_size
                d2c_coeffs(i+j)= d2c_coeffs(i+j)+d1sB_coeffs(i)*d1sA_coeffs(j) * geom_A_coeffs(i)
            end do
        end do

        do i = 1, B_size
            do j = -A_size, -1
                d2c_coeffs(i+j)= d2c_coeffs(i+j)+d1sB_coeffs(i)*d1sA_coeffs(j) * geom_A_coeffs(i)
            end do
            do j = 1, A_size
                d2c_coeffs(i+j-1)= d2c_coeffs(i+j-1)+d1sB_coeffs(i)*d1sA_coeffs(j) * geom_A_coeffs(i)
            end do
        end do

        d2c_coeffs=d2c_coeffs*geom_B_coeff

    end subroutine

    subroutine set_d2c_coeffs3(d1sA_coeffs, d1sB_coeffs, d2c_coeffs, A_size, A_size_max, B_size, geom_A_coeffs, geom_B_coeff)

        implicit none
        integer B_size, A_size_max, i, j
        integer, dimension(-B_size:B_size)  ::   A_size
        real(kind=8), dimension(-B_size:B_size,-A_size_max:A_size_max) :: d1sA_coeffs
        real(kind=8), dimension(-B_size:B_size)   :: geom_A_coeffs
        real(kind=8), dimension(-B_size:B_size)   :: d1sB_coeffs
        real(kind=8), dimension(-(A_size(-B_size)+B_size)+1:(A_size(B_size)+B_size)-1)    :: d2c_coeffs
        real(kind=8) ::geom_B_coeff

        d2c_coeffs=0.d0

        do i = -B_size, -1
            do j = -A_size(i), -1
                d2c_coeffs(i+j+1)= d2c_coeffs(i+j+1)+d1sB_coeffs(i) *d1sA_coeffs(i,j)* geom_A_coeffs(i)
            end do
            do j = 1, A_size(i)
                d2c_coeffs(i+j)= d2c_coeffs(i+j)+d1sB_coeffs(i)*d1sA_coeffs(i,j) * geom_A_coeffs(i)
            end do
        end do

        do i = 1, B_size
            do j = -A_size(i), -1
                d2c_coeffs(i+j)= d2c_coeffs(i+j)+d1sB_coeffs(i)*d1sA_coeffs(i,j) * geom_A_coeffs(i)
            end do
            do j = 1, A_size(i)
                d2c_coeffs(i+j-1)= d2c_coeffs(i+j-1)+d1sB_coeffs(i) * d1sA_coeffs(i,j) * geom_A_coeffs(i)
            end do
        end do

        d2c_coeffs=d2c_coeffs*geom_B_coeff

    end subroutine

end module d2c_from_d1s


module optimized_schemes

    use boundaries_types
    use schemes_settings
    implicit none
    real(kind=8), allocatable, dimension(:) :: m_loc, bp_loc, upper_loc

    real(kind=8), dimension(:,:), allocatable :: A
    real(kind=8)    :: B(2), Bs(-1:3), Be(-1:3)

    integer :: last_nmax=0, last_bc_type1, last_bc_type2
    logical   :: last_shifted, update_matrix
    real(kind=8)   :: last_h

contains

    ! Schéma le plus simple: d1f = B*f ou B est une matrice tridiagonale
    subroutine O2(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d1f(:)

        real(kind=8) B
        integer:: i

        B=0.5d0/h

        do i=2,n_max-2
            d1f(i)= ( f(i+1) - f(i-1) ) * B
        enddo

        d1f(n_max-1)= ( f(n_max) - f(n_max-2) ) * B

        return

    end subroutine

    ! Schéma compact: d1f est donné par [A]*d1f=[B]*f où A est une matrice tridiagonale

    subroutine CPT(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d1f(:)

        integer :: i
        real(kind=8)  :: temp

        ! Calcul des matrices necessaires à l'algorithme -----------------------------------------------------------
        ! S'execute uniquement lorsque le nombre de points sur lequel on applique le schéma varie d'une exécution à l'autre...
        if (last_nmax/=n_max) then
            call fill_matrix_opt(n_max, bc_type, shifted, h)
            write(*,*)'REMPLISSAGE DES MATRICES'
        end if

        last_nmax=n_max

        ! SCHEME CALCULATION ---------------------------------------------------------------
        ! [A]*d1f = [B]*f : d1f est obtenu par application d'une méthode LU constituée de deux étapes (forward et backward)


        d1f(2)= Bs(-1)*f(1) + Bs(0) * f(2)          &
        + Bs(1) * f(3) + Bs(2) * f(4)        &
        + Bs(3) * f(5)


        ! forward step
        temp=d1f(2)
        do i=3, n_max-2
            d1f(i)=     B(2)*( f(i+2) - f(i-2))   &
            +   B(1)*( f(i+1) - f(i-1) ) - m_loc(i-1)*temp ! ici remplacer d1f par f multiplie la vitesse par 2! (0.32 ->0.14)
            temp=d1f(i)
        enddo

        d1f(n_max-1)= (-( + Be(-1)*f(n_max) + Be(0)*f(n_max-1)     &
        + Be(1)*f(n_max-2)+ Be(2)*f(n_max-3)      &
        + Be(3)*f(n_max-4) ) - m_loc(n_max-2)*temp)*bp_loc(n_max-2)



        ! backward step
        temp=d1f(n_max-1)
        do i = n_max-2, 2, -1
            d1f(i) = (d1f(i) - upper_loc(i-1)*temp)*bp_loc(i-1)
            temp=d1f(i)
        enddo


        return

    end subroutine

    subroutine CPT_DGTTRF(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h,f(:)
        real(kind=8), intent(out)     :: d1f(:)
        real*8, dimension(:), save, allocatable    :: DL, D, DU
        real(kind=8),dimension (:), save, allocatable    :: ipiv, du2
        integer     :: info

        integer :: i
        real(kind=8)  :: temp


        !write(*,*)'D1_OUCS3'

        if (last_nmax/=n_max) then
            if (allocated(D)) deallocate(DL, D, DU, DU2, ipiv)
            allocate(DL(n_max-2), D(n_max-2), DU(n_max-2), DU2(n_max-2), ipiv(n_max-2))

            ! Definition de la matrice B -----------------------------------------------------------
            ! interieur de la matrice
            B(1)=1.57557379d0 / (2.d0*h)
            B(2)= 0.1832051925d0 / (4.d0*h)

            ! 1 ere ligne
            Bs(-1)=-0.35d0             / h
            Bs(0)=-0.43333333333333d0  / h
            Bs(1)=0.9d0                / h
            Bs(2)=-0.1d0               / h
            Bs(3)=-0.016666666666667d0 / h

            ! dernière ligne -----------------------
            Be(-1) =   -0.27333333333333d0 / h
            Be(0) =    -0.74d0             / h
            Be(1) =    1.36d0              / h
            Be(2) =    -0.40666666666667d0 / h
            Be(3) =    0.06d0              / h

            ! Definition de la matrice A -----------------------------------------------------------
                ! DL est la diagonale inférieure de A
                ! D la diagonale
                ! et DU la diagonale supérieure
            DL=0.3793894912d0
            D=1.d0
            DU=0.3793894912d0

            DU(1)=0.d0
            DL(n_max-3)=0.d0


            call DGTTRF( n_max-2, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), du2(1:n_max-4), ipiv, info )

        end if

        last_nmax=n_max



        ! Calcul du terme droit d=B*f ---------------------------------------------------------------

        d1f(2)= Bs(-1)*f(1) + Bs(0) * f(2)          &
        + Bs(1) * f(3) + Bs(2) * f(4)        &
        + Bs(3) * f(5)

        do i=3, n_max-2
            d1f(i)=     B(2)*( f(i+2) - f(i-2))   &
            +   B(1)*( f(i+1) - f(i-1) )
        enddo

        d1f(n_max-1)= -( + Be(-1)*f(n_max) + Be(0)*f(n_max-1)     &
        + Be(1)*f(n_max-2)+ Be(2)*f(n_max-3)      &
        + Be(3)*f(n_max-4) )

        ! Resolution du systeme [A] d1f = d
        call DGTTS2( 0, n_max-2, 1, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), du2(1:n_max-4), IPIV, d1f(2:n_max-1), n_max-2 )
            !call DGTTRS( 'N', n_max-2, 1, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), du2(1:n_max-4), ipiv, d1f(2:n_max-1), n_max-2, info )
            !call Lapack_solver(DL, D, DU, d1f(2:n_max-1), d1f(2:n_max-1), n_max-2)


        return

    end subroutine

#ifdef IFORT
    subroutine CPT_DDTTRFB(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h,f(:)
        real(kind=8), intent(out)     :: d1f(:)
        real*8, dimension(:), save, allocatable    :: DL, D, DU
        real(kind=8),dimension (:), save, allocatable    :: ipiv, du2
        integer     :: info

        integer :: i
        real(kind=8)  :: temp


        !write(*,*)'D1_OUCS3'

        if (last_nmax/=n_max) then
            if (allocated(D)) deallocate(DL, D, DU, DU2, ipiv)
            allocate(DL(n_max-2), D(n_max-2), DU(n_max-2), DU2(n_max-2), ipiv(n_max-2))

            ! Definition de la matrice B -----------------------------------------------------------
            ! interieur de la matrice
            B(1)=1.57557379d0 / (2.d0*h)
            B(2)= 0.1832051925d0 / (4.d0*h)

            ! 1 ere ligne
            Bs(-1)=-0.35d0             / h
            Bs(0)=-0.43333333333333d0  / h
            Bs(1)=0.9d0                / h
            Bs(2)=-0.1d0               / h
            Bs(3)=-0.016666666666667d0 / h

            ! dernière ligne -----------------------
            Be(-1) =   -0.27333333333333d0 / h
            Be(0) =    -0.74d0             / h
            Be(1) =    1.36d0              / h
            Be(2) =    -0.40666666666667d0 / h
            Be(3) =    0.06d0              / h

            ! Definition de la matrice A -----------------------------------------------------------
                ! DL est la diagonale inférieure de A
                ! D la diagonale
                ! et DU la diagonale supérieure
            DL=0.3793894912d0
            D=1.d0
            DU=0.3793894912d0

            DU(1)=0.d0
            DL(n_max-3)=0.d0

            call ddttrfb( n_max-2, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), info)
            !call DGTTRF( n_max-2, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), du2(1:n_max-4), ipiv, info )

        end if

        last_nmax=n_max



        ! Calcul du terme droit d=B*f ---------------------------------------------------------------

        d1f(2)= Bs(-1)*f(1) + Bs(0) * f(2)          &
        + Bs(1) * f(3) + Bs(2) * f(4)        &
        + Bs(3) * f(5)

        do i=3, n_max-2
            d1f(i)=     B(2)*( f(i+2) - f(i-2))   &
            +   B(1)*( f(i+1) - f(i-1) )
        enddo

        d1f(n_max-1)= -( + Be(-1)*f(n_max) + Be(0)*f(n_max-1)     &
        + Be(1)*f(n_max-2)+ Be(2)*f(n_max-3)      &
        + Be(3)*f(n_max-4) )

        ! Resolution du systeme [A] d1f = d
        !call DGTTS2( 0, n_max-2, 1, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), du2(1:n_max-4), IPIV, d1f(2:n_max-1), n_max-2 )
        call ddttrsb( 'N', n_max-2, 1, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), d1f(2:n_max-1), n_max-2, info )
            !call DGTTRS( 'N', n_max-2, 1, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), du2(1:n_max-4), ipiv, d1f(2:n_max-1), n_max-2, info )
            !call Lapack_solver(DL, D, DU, d1f(2:n_max-1), d1f(2:n_max-1), n_max-2)


        return

    end subroutine

#endif

    subroutine DRP7(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d1f(:)

        integer :: i
        real(kind=8), allocatable, save :: A(:)
        integer, save   :: last_nmax=0

        if (last_nmax/=n_max) then

            if (allocated(A)) deallocate(A)
            allocate(A(7))
            A(1) =  0.95396219562045d0    /h
            A(2) =  -0.41234590494721d0   /h
            A(3) =  0.21233981563217d0    /h
            A(4) =  -0.10672135533957d0   /h
            A(5) =  0.046775295199319d0   /h
            A(6) = -0.015323662211638d0   /h
            A(7) = 0.0026664396358809d0   /h
        end if

        last_nmax=n_max

        do i=8,n_max-8
            d1f(i)= A(7)*(f(i+7) - f(i-7))                  &
            + A(6)*(f(i+6) - f(i-6))     &
            + A(5)*(f(i+5) - f(i-5))     &
            + A(4)*(f(i+4) - f(i-4))     &
            + A(3)*(f(i+3) - f(i-3))     &
            + A(2)*(f(i+2) - f(i-2))     &
            + A(1)*(f(i+1) - f(i-1))
        enddo


        call D1c_CTR_CLOSURE_INF(f, d1f, h, 7)
        call D1c_CTR_CLOSURE_SUP(f, d1f, n_max, h, n_max-7, shifted)

        return

    end subroutine



    subroutine DRP7_blas(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d1f(:)

        integer :: j
        integer, save   :: last_nmax=0

        ! BLAS data
        integer, parameter :: ml=7, mu=7
        real*8, dimension(:,:), allocatable, save       :: A, matrix

        integer, save     :: N, M
        real*8  :: alpha=1.d0, beta=0.d0
        real*8  :: a1, a2, a3, a4, a5, a6, a7
        real*8  :: b11, b21, b22, b31, b32, b33, b41, b42, b43, b44, b51, b52, b53, b54, b55
        real*8  :: b61, b62, b63, b64, b65, b66


        !call example1
        !call example2
        !call example3
        !call example4
        !call example5
        !call example6
        !call exit




        if (last_nmax/=n_max) then

            if (allocated(A)) deallocate(A, matrix)
            allocate(A(ml+mu+1, n_max))
            allocate(matrix(n_max, n_max))

            N=n_max
            M=n_max

            b11=0.5d0 /h; b21=2.d0/(3.d0) /h; b22=-1.d0/(12.d0) /h
            b31=0.79926643d0 /h; b32=-0.18941314d0 /h; b33=0.02651995d0 /h
            b51=0.88149995153015d0 /h; b52=-0.29786830893263d0 /h; b53= 0.098184076271814d0 /h
            b54=-0.02388959221513d0 /h; b55=0.0030486998975861d0 /h
            b61= 0.91595792650492d0/h; b62=- 0.34922251636223d0/h; b63= 0.14398145036906d0/h
            b64=- 0.051236991729043d0/h; b65= 0.013273181125903d0/h; b66=- 0.0018126562894445d0/h;


            a1= 0.95396219562045d0/h; a2=- 0.41234590494721d0/h; a3= 0.21233981563217d0/h
            a4=- 0.10672135533957d0/h; a5= 0.046775295199319d0/h; a6=- 0.015323662211638d0/h; a7= 0.0026664396358809d0/h;

            matrix=0.d0

            j=2
            matrix(j,j+1)=b11
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!

            j=3
            matrix(j,j+1)=b21
            matrix(j,j+2)=b22
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!

            j=4
            matrix(j,j+1)=b31
            matrix(j,j+2)=b32
            matrix(j,j+3)=b33
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!

            j=5
            matrix(j,j+1)=b31
            matrix(j,j+2)=b32
            matrix(j,j+3)=b33
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!

            j=6
            matrix(j,j+1)=b51
            matrix(j,j+2)=b52
            matrix(j,j+3)=b53
            matrix(j,j+4)=b54
            matrix(j,j+5)=b55
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-4)=-matrix(j,j+4) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-5)=-matrix(j,j+5) !!!!!!!!!!!!!!!!!!!

            j=7
            matrix(j,j+1)=b61
            matrix(j,j+2)=b62
            matrix(j,j+3)=b63
            matrix(j,j+4)=b64
            matrix(j,j+5)=b65
            matrix(j,j+6)=b66
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-4)=-matrix(j,j+4) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-5)=-matrix(j,j+5) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-6)=-matrix(j,j+6) !!!!!!!!!!!!!!!!!!!


            do j = 8, M-7
                matrix(j,j+1)=a1
                matrix(j,j+2)=a2
                matrix(j,j+3)=a3
                matrix(j,j+4)=a4
                matrix(j,j+5)=a5
                matrix(j,j+6)=a6
                matrix(j,j+7)=a7
                matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
                matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
                matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!
                matrix(j,j-4)=-matrix(j,j+4) !!!!!!!!!!!!!!!!!!!
                matrix(j,j-5)=-matrix(j,j+5) !!!!!!!!!!!!!!!!!!!
                matrix(j,j-6)=-matrix(j,j+6) !!!!!!!!!!!!!!!!!!!
                matrix(j,j-7)=-matrix(j,j+7) !!!!!!!!!!!!!!!!!!!

            end do

            j=M-6
            matrix(j,j+1)=b61
            matrix(j,j+2)=b62
            matrix(j,j+3)=b63
            matrix(j,j+4)=b64
            matrix(j,j+5)=b65
            matrix(j,j+6)=b66
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-4)=-matrix(j,j+4) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-5)=-matrix(j,j+5) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-6)=-matrix(j,j+6) !!!!!!!!!!!!!!!!!!!

            j=M-5
            matrix(j,j+1)=b51
            matrix(j,j+2)=b52
            matrix(j,j+3)=b53
            matrix(j,j+4)=b54
            matrix(j,j+5)=b55
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-4)=-matrix(j,j+4) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-5)=-matrix(j,j+5) !!!!!!!!!!!!!!!!!!!

            j=M-4
            matrix(j,j+1)=b31
            matrix(j,j+2)=b32
            matrix(j,j+3)=b33
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!

            j=M-3
            matrix(j,j+1)=b31
            matrix(j,j+2)=b32
            matrix(j,j+3)=b33
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!

            j=M-2
            matrix(j,j+1)=b21
            matrix(j,j+2)=b22
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!

            j=M-1
            matrix(j,j+1)=b11
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!

            j=M
            matrix(j,j)=0.d0

            call create_bs_matrix(A, matrix, ml, mu, N, M)

        end if

        last_nmax=n_max

        !f=1.d0

        !CALL DGBMV( 'N' , M , N , ML , MU ,  alpha , A , ml+mu+1 , f , 1  , beta , d1f , 1)
        CALL DGBMV( 'N' , N , N , ML , MU ,  1.d0 , A , ml+mu+1 , f , 1  , 0.d0 , d1f , 1)


        return

    contains

        subroutine create_bs_matrix(A, matrix, kl, ku, N, M)
            implicit none
            integer :: kl, ku, N, M
            real*8, dimension(kl+ku+1, N)       :: A
            real*8, dimension(M, N)       :: matrix
            integer     :: i, j, k

            A=0.d0
            do J = 1, N
                K = KU + 1 - J
                do I = MAX( 1, J - KU ), MIN( M, J + KL )
                    A( K + I, J ) = matrix( I, J )
                enddo
            enddo

        end subroutine create_bs_matrix

    end subroutine

    subroutine D1c_CTR_CLOSURE_INF(f, d1f, h, imax)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(imax)

        implicit none

        integer :: i, n_max, imax
        real(kind=8) h
        real(kind=8) d1f(:)
        real(kind=8) f(:)

        real(kind=8) d1f_temp(CPT_MIN)


        if (treat_boundary_by_cpt_scheme) then

            call CPT(  &
            f(1:CPT_MIN), d1f_temp(1:CPT_MIN), CPT_MIN, h, .true.,      &
            Dirichlet)

            d1f(2)=d1f_temp(2)
            if (imax.eq.2) return

            d1f(3)=d1f_temp(3)
            if (imax.eq.3) return

        else

            ! i=2---------------------------------------------------
            i=2
            d1f(i)= (f(i+1)-f(i-1)) /   (2.d0*h)

            if (i.eq.imax) return
            ! i=3---------------------------------------------------
            i=3
            d1f(i)    =   2.d0*(f(i+1)-f(i-1))  / (3.d0*h)      &
            -   (f(i+2) - f(i-2))         / (12.d0*h)

            if (i.eq.imax) return

        endif

        ! i=4---------------------------------------------------
        i=4
        d1f(i)  =   0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h

        if (i.eq.imax) return

        ! i=5---------------------------------------------------
        i=5
        d1f(i)  =   0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h

        if (i.eq.imax) return

        ! i=6---SEE D1_EXP_5pts.wxm-----------------------------
        i=6
        d1f(i)  =   0.88149995153015d0*(f(i+1)-f(i-1))      / h     &
        -   0.29786830893263d0*(f(i+2)-f(i-2))              / h     &
        +   0.098184076271814d0* (f(i+3) - f(i-3))          / h     &
        -   0.02388959221513d0* (f(i+4) - f(i-4))           / h     &
        +   0.0030486998975861d0* (f(i+5) - f(i-5))         / h

        if (i.eq.imax) return

        ! i=7---SEE D1_EXP_6pts.wxm-----------------------------
        i=7
        d1f(i)  =   0.91595792650492d0  * (f(i+1)-f(i-1))   / h     &
        -   0.34922251636223d0  * (f(i+2)-f(i-2))           / h     &
        +   0.14398145036906d0  * (f(i+3) - f(i-3))         / h     &
        -   0.051236991729043d0 * (f(i+4) - f(i-4))         / h     &
        +   0.013273181125903d0 * (f(i+5) - f(i-5))         / h     &
        -   0.0018126562894445d0* (f(i+6) - f(i-6))         / h

        return

    end subroutine

    subroutine D1c_CTR_CLOSURE_SUP(f, d1f, n_max, h, imin, shifted)

        implicit none

        integer :: i, n_max, imin, j
        real(kind=8) h, fn
        real(kind=8) d1f(:)
        real(kind=8) f(:)
        logical :: shifted

        real(kind=8) d1f_temp(n_max-CPT_MIN:n_max-1)

        if (shifted) then
            i=n_max-2
        else
            i=n_max-1
        end if

        if (treat_boundary_by_cpt_scheme) then

            call CPT(f(n_max-CPT_MIN:n_max), d1f_temp, CPT_MIN+1, h, shifted, Dirichlet)

            d1f(i)=d1f_temp(i)
            if (i.eq.imin) return

            i=i-1
            d1f(i)=d1f_temp(i)
            if (i.eq.imin) return

        else

            ! First point-------------------------------------------
            d1f(i)= (f(i+1)-f(i-1)) /   (2.d0*h)

            if (i.eq.imin) return
            ! 2nd point----------------------------------------------
            i=i-1
            d1f(i)    =   2.d0*(f(i+1)-f(i-1))  / (3.d0*h)      &
            -   (f(i+2) - f(i-2))         / (12.d0*h)

            if (i.eq.imin) return

        endif

        ! 3rd point---------------------------------------------
        i=i-1
        d1f(i)  =   0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h

        if (i.eq.imin) return

        ! 4th point---------------------------------------------
        i=i-1
        d1f(i)  =   0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h

        if (i.eq.imin)return

        ! 5th point ----D1_EXP_5pts.wxm------------------------
        i=i-1
        d1f(i)  =   0.88149995153015d0*(f(i+1)-f(i-1))      / h     &
        -   0.29786830893263d0*(f(i+2)-f(i-2))              / h     &
        +   0.098184076271814d0* (f(i+3) - f(i-3))          / h     &
        -   0.02388959221513d0* (f(i+4) - f(i-4))           / h     &
        +   0.0030486998975861d0* (f(i+5) - f(i-5))         / h

        if (i.eq.imin) return

        ! 6th point ----D1_EXP_6pts.wxm------------------------
        i=i-1
        d1f(i)  =   0.91595792650492d0  * (f(i+1)-f(i-1))   / h     &
        -   0.34922251636223d0  * (f(i+2)-f(i-2))           / h     &
        +   0.14398145036906d0  * (f(i+3) - f(i-3))         / h     &
        -   0.051236991729043d0 * (f(i+4) - f(i-4))         / h     &
        +   0.013273181125903d0 * (f(i+5) - f(i-5))         / h     &
        -   0.0018126562894445d0* (f(i+6) - f(i-6))         / h

        if (i.eq.imin) return

        ! 7th point ---------------------------------------------
        i=i-1
        d1f(i)  =   0.95396219562045d0  * (f(i+1)-f(i-1))   / h     &
        -   0.41234590494721d0  * (f(i+2)-f(i-2))           / h     &
        +   0.21233981563217d0  * (f(i+3) - f(i-3))         / h     &
        -   0.10672135533957d0 * (f(i+4) - f(i-4))          / h     &
        +   0.046775295199319d0 * (f(i+5) - f(i-5))         / h     &
        -   0.015323662211638d0* (f(i+6) - f(i-6))          / h     &
        +   0.0026664396358809d0* (f(i+7) - f(i-7))          / h


        return

    end subroutine

    subroutine TS_init_loc(A, B, C, n_max)

        ! a - sub-diagonal (means it is the diagonal below the main diagonal)
        ! b - the main diagonal
        ! c - sup-diagonal (means it is the diagonal above the main diagonal)
        ! v - right part
        ! x - the answer
        ! n_max - number of equations

        implicit none
        integer,intent(in)                      :: n_max
        real(8),dimension(n_max),intent(in)     :: A, B, C
        integer i

        if (allocated(m_loc)) then
            deallocate(m_loc, bp_loc, upper_loc)
        end if

        allocate(m_loc(n_max), bp_loc(n_max), upper_loc(n_max))

        bp_loc(1) = 1.d0/B(1)
        upper_loc(1)=C(1)

        !   The first pass (setting coefficients):
        do i = 2,n_max
            m_loc(i) = A(i)*bp_loc(i-1)
            bp_loc(i) = 1.d0/(B(i) - m_loc(i)*C(i-1))
            upper_loc(i)=C(i)
        enddo

    end subroutine

    subroutine fill_matrix_opt(n_max, bc_type, shifted, h)
        implicit none
        integer :: i, n_max, bc_type
        logical :: shifted
        real(kind=8)  :: h

        if (allocated(A)) then
            deallocate(A)
        end if

        allocate(A(-1:1, n_max))

        do i=3, n_max-3
            A(-1,i) = 0.3793894912d0
            A(0,i) = 1.d0
            A(1,i) = 0.3793894912d0
        enddo


        A(-1,2)=0.d0
        A(0,2)=1.d0
        A(1,2)=0.d0

        B(1)=1.57557379d0       / (2.d0*h)
        B(2)= 0.1832051925d0    / (4.d0*h)

        Bs(-1)=-0.35d0             / h
        Bs(0)=-0.43333333333333d0  / h
        Bs(1)=0.9d0                / h
        Bs(2)=-0.1d0               / h
        Bs(3)=-0.016666666666667d0 / h

        ! End coefficients -----------------------
        Be(-1) =   -0.27333333333333d0 / h
        Be(0) =    -0.74d0             / h
        Be(1) =    1.36d0              / h
        Be(2) =    -0.40666666666667d0 / h
        Be(3) =    0.06d0              / h

        A(-1,n_max-2) = 0.3793894912d0
        A(0,n_max-2) = 1.d0
        A(1,n_max-2) = 0.3793894912d0

        A(-1,n_max-1)=0.d0
        A(0,n_max-1)=1.d0
        A(1,n_max-1)=0.d0

        call TS_init_loc(A(-1,2:n_max-1), A(0,2:n_max-1), A(1,2:n_max-1), n_max-2)


    end subroutine fill_matrix_opt

end module optimized_schemes


module optimized_schemes_multiple_RHS

    use boundaries_types
    use schemes_settings
    implicit none
    real(kind=8), allocatable, dimension(:) :: m_loc, bp_loc, upper_loc

    real(kind=8), dimension(:,:), allocatable :: A
    real(kind=8)    :: B(3), Bs(-1:3), Be(-1:3)

    integer :: last_nmax=0, last_bc_type1, last_bc_type2
    logical   :: last_shifted, update_matrix
    real(kind=8)   :: last_h

contains

    ! Schéma le plus simple: d1f = B*f ou B est une matrice tridiagonale
    subroutine mO2(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d1f(:)

        real(kind=8) B
        integer:: i

        B=0.5d0/h

        do i=2,n_max-2
            d1f(i)= ( f(i+1) - f(i-1) ) * B
        enddo

        d1f(n_max-1)= ( f(n_max) - f(n_max-2) ) * B

        return

    end subroutine

    ! Schéma compact: d1f est donné par [A]*d1f=[B]*f où A est une matrice tridiagonale

    subroutine mCPT(f, d1f, n_max, h, shifted, bc_type, nrhs)

        implicit none

        integer, intent(in)     :: n_max, bc_type, nrhs
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:, :)
        real(kind=8), intent(out)     :: d1f(:, :)

        integer :: i, r
        real(kind=8)  :: temp

        ! Calcul des matrices necessaires à l'algorithme -----------------------------------------------------------
        ! S'execute uniquement lorsque le nombre de points sur lequel on applique le schéma varie d'une exécution à l'autre...
        if (last_nmax/=n_max) then
            call fill_matrix_opt(n_max, bc_type, shifted, h)
            write(*,*)'REMPLISSAGE DES MATRICES'
        end if

        last_nmax=n_max

        ! SCHEME CALCULATION ---------------------------------------------------------------
        ! [A]*d1f = [B]*f : d1f est obtenu par application d'une méthode LU constituée de deux étapes (forward et backward)

        do r = 1, nrhs

            d1f(2, r)= Bs(-1)*f(1, r) + Bs(0) * f(2, r)          &
            + Bs(1) * f(3, r) + Bs(2) * f(4, r)        &
            + Bs(3) * f(5, r)


            ! forward step
            temp=d1f(2, r)
            do i=3, n_max-2
                d1f(i, r)=     B(2)*( f(i+2, r) - f(i-2, r))   &
                +   B(1)*( f(i+1, r) - f(i-1, r) ) - m_loc(i-1)*temp ! ici remplacer d1f par f multiplie la vitesse par 2! (0.32 ->0.14)
                temp=d1f(i, r)
            enddo

            d1f(n_max-1, r)= (-( + Be(-1)*f(n_max, r) + Be(0)*f(n_max-1, r)     &
            + Be(1)*f(n_max-2, r)+ Be(2)*f(n_max-3, r)      &
            + Be(3)*f(n_max-4, r) ) - m_loc(n_max-2)*temp)*bp_loc(n_max-2)



            ! backward step
            temp=d1f(n_max-1, r)
            do i = n_max-2, 2, -1
                d1f(i, r) = (d1f(i, r) - upper_loc(i-1)*temp)*bp_loc(i-1)
                temp=d1f(i, r)
            enddo

        end do


        return

    end subroutine

    ! UN SEUL RHS dans LAPACK

#ifdef IFORT

    subroutine mCPT_DDTTRFB(f, d1f, n_max, h, shifted, bc_type, nrhs)

        implicit none

        integer, intent(in)     :: n_max, bc_type, nrhs
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h,f(:,:)
        real(kind=8), intent(out)     :: d1f(:,:)
        real*8, dimension(:), save, allocatable    :: DL, D, DU
        real(kind=8),dimension (:), save, allocatable    :: ipiv, du2
        integer     :: info

        integer :: i, r
        real(kind=8)  :: temp


        !write(*,*)'D1_OUCS3'

        if (last_nmax/=n_max) then
            if (allocated(D)) deallocate(DL, D, DU, DU2, ipiv)
            allocate(DL(n_max-2), D(n_max-2), DU(n_max-2), DU2(n_max-2), ipiv(n_max-2))

            ! Definition de la matrice B -----------------------------------------------------------
            ! interieur de la matrice
            B(1)=1.57557379d0 / (2.d0*h)
            B(2)= 0.1832051925d0 / (4.d0*h)

            ! 1 ere ligne
            Bs(-1)=-0.35d0             / h
            Bs(0)=-0.43333333333333d0  / h
            Bs(1)=0.9d0                / h
            Bs(2)=-0.1d0               / h
            Bs(3)=-0.016666666666667d0 / h

            ! dernière ligne -----------------------
            Be(-1) =   -0.27333333333333d0 / h
            Be(0) =    -0.74d0             / h
            Be(1) =    1.36d0              / h
            Be(2) =    -0.40666666666667d0 / h
            Be(3) =    0.06d0              / h

            ! Definition de la matrice A -----------------------------------------------------------
                ! DL est la diagonale inférieure de A
                ! D la diagonale
                ! et DU la diagonale supérieure
            DL=0.3793894912d0
            D=1.d0
            DU=0.3793894912d0

            DU(1)=0.d0
            DL(n_max-3)=0.d0

            call ddttrfb( n_max-2, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), info)
            !call DGTTRF( n_max-2, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), du2(1:n_max-4), ipiv, info )

        end if

        last_nmax=n_max



        ! Calcul du terme droit d=B*f ---------------------------------------------------------------

        do r=1, nrhs

            do i=3, n_max-2
                d1f(i,r)=     B(2)*( f(i+2,r) - f(i-2,r))   &
                +   B(1)*( f(i+1,r) - f(i-1,r) )
            enddo

        enddo

        do r=1, nrhs

            d1f(2,r)= Bs(-1)*f(1,r) + Bs(0) * f(2,r)          &
            + Bs(1) * f(3,r) + Bs(2) * f(4,r)        &
            + Bs(3) * f(5,r)

            d1f(n_max-1,r)= -( + Be(-1)*f(n_max,r) + Be(0)*f(n_max-1,r)     &
            + Be(1)*f(n_max-2,r)+ Be(2)*f(n_max-3,r)      &
            + Be(3)*f(n_max-4,r) )


            call ddttrsb( 'N', n_max-2, 1, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), d1f(2:n_max-1, r), n_max-2, info )

        enddo
        ! Resolution du systeme [A] d1f = d
        !call DGTTS2( 0, n_max-2, 1, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), du2(1:n_max-4), IPIV, d1f(2:n_max-1), n_max-2 )
        !call ddttrsb( 'N', n_max-2, nrhs, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), d1f(2:n_max-1, 1:nrhs), n_max-2, info )
            !call DGTTRS( 'N', n_max-2, 1, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), du2(1:n_max-4), ipiv, d1f(2:n_max-1), n_max-2, info )
            !call Lapack_solver(DL, D, DU, d1f(2:n_max-1), d1f(2:n_max-1), n_max-2)


        return

    end subroutine


    ! PLUSIEURS RHS dans LAPACK
    subroutine mCPT_DDTTRFB2(f, d1f, n_max, h, shifted, bc_type, nrhs)

        implicit none

        integer, intent(in)     :: n_max, bc_type, nrhs
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h,f(:,:)
        real(kind=8), intent(out)     :: d1f(:,:)
        real*8, dimension(:), save, allocatable    :: DL, D, DU
        real(kind=8),dimension (:), save, allocatable    :: ipiv, du2
        integer     :: info

        integer :: i, r
        real(kind=8)  :: temp
        real*8  :: B1(-1:3)


        !write(*,*)'D1_OUCS3'

        if (last_nmax/=n_max) then
            if (allocated(D)) deallocate(DL, D, DU, DU2, ipiv)
            allocate(DL(n_max-2), D(n_max-2), DU(n_max-2), DU2(n_max-2), ipiv(n_max-2))

            ! Definition de la matrice B -----------------------------------------------------------
            ! interieur de la matrice
            B(1)=1.57557379d0 / (2.d0*h)
            B(2)= 0.1832051925d0 / (4.d0*h)

            ! 1 ere ligne
            B1(-1)=-0.35d0             / h
            B1(0)=-0.43333333333333d0  / h
            B1(1)=0.9d0                / h
            B1(2)=-0.1d0               / h
            B1(3)=-0.016666666666667d0 / h


            ! Definition de la matrice A -----------------------------------------------------------
                ! DL est la diagonale inférieure de A
                ! D la diagonale
                ! et DU la diagonale supérieure
            DL=0.3793894912d0
            D=1.d0
            DU=0.3793894912d0

            DU(1)=0.d0
            DL(n_max-3)=0.d0

            call ddttrfb( n_max-2, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), info)
            !call DGTTRF( n_max-2, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), du2(1:n_max-4), ipiv, info )

        end if

        last_nmax=n_max



        ! Calcul du terme droit d=B*f ---------------------------------------------------------------

        do r=1, nrhs

            do i=3, n_max-2
                d1f(i,r)=     B(2)*( f(i+2,r) - f(i-2,r))   &
                +   B(1)*( f(i+1,r) - f(i-1,r) )
            enddo

            !call ddttrsb( 'N', n_max-2, 1, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), d1f(2:n_max-1, r), n_max-2, info )

        enddo

        ! Plus optimisé si on fait une boucle pour les conditions aux limites
        do r=1, nrhs

            d1f(2,r)= B1(-1)*f(1,r) + B1(0) * f(2,r)          &
            + B1(1) * f(3,r) + B1(2) * f(4,r)        &
            + B1(3) * f(5,r)

            d1f(n_max-1,r)= -( B1(-1)*f(n_max,r) + B1(0)*f(n_max-1,r)     &
            + B1(1)*f(n_max-2,r)+ B1(2)*f(n_max-3,r)      &
            + B1(3)*f(n_max-4,r) )


            !call ddttrsb( 'N', n_max-2, 1, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), d1f(2:n_max-1, r), n_max-2, info )

        enddo
        ! Resolution du systeme [A] d1f = d
        !call DGTTS2( 0, n_max-2, 1, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), du2(1:n_max-4), IPIV, d1f(2:n_max-1), n_max-2 )
        call ddttrsb( 'N', n_max-2, nrhs, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), d1f(2:n_max-1, 1:nrhs), n_max-2, info )
            !call DGTTRS( 'N', n_max-2, 1, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), du2(1:n_max-4), ipiv, d1f(2:n_max-1), n_max-2, info )
            !call Lapack_solver(DL, D, DU, d1f(2:n_max-1), d1f(2:n_max-1), n_max-2)


        return

    end subroutine


    ! PLUSIEURS RHS dans LAPACK
    subroutine mKIM_DDTTRFB2(f, d1f, n_max, h, shifted, bc_type, nrhs)

        implicit none

        integer, intent(in)     :: n_max, bc_type, nrhs
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h,f(:,:)
        real(kind=8), intent(out)     :: d1f(:,:)
        real*8, dimension(:), save, allocatable    :: DL, D, DU
        real(kind=8),dimension (:), save, allocatable    :: ipiv, du2
        integer     :: info

        integer :: i, r
        real(kind=8)  :: temp
        real*8  :: B1(-1:3), B2(2)


        !write(*,*)'D1_OUCS3'

        if (last_nmax/=n_max) then
            if (allocated(D)) deallocate(DL, D, DU, DU2, ipiv)
            allocate(DL(n_max-2), D(n_max-2), DU(n_max-2), DU2(n_max-2), ipiv(n_max-2))

            ! Definition de la matrice B -----------------------------------------------------------
            ! interieur de la matrice

            B(1)=1.568098212d0 / (2.d0*h)
            B(2)= 0.271657107d0 / (4.d0*h)
            B(3)= -0.022576781d0 / (6.d0*h)

            ! 1 ere ligne
            B1(-1)=-0.35d0             / h
            B1(0)=-0.43333333333333d0  / h
            B1(1)=0.9d0                / h
            B1(2)=-0.1d0               / h
            B1(3)=-0.016666666666667d0 / h


            B2(1)=14.0/ (18.d0*h)
            B2(2)= 1.d0 / (36.d0*h)

            ! Definition de la matrice A -----------------------------------------------------------
                ! DL est la diagonale inférieure de A
                ! D la diagonale
                ! et DU la diagonale supérieure
            DL=0.408589269d0
            D=1.d0
            DU=0.408589269d0

            DU(1)=0.d0
            DU(2)=0.3333333d0
            DL(n_max-4)=0.3333333d0
            DL(n_max-3)=0.d0

            call ddttrfb( n_max-2, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), info)
            !call DGTTRF( n_max-2, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), du2(1:n_max-4), ipiv, info )

        end if

        last_nmax=n_max



        ! Calcul du terme droit d=B*f ---------------------------------------------------------------

        do r=1, nrhs

            do i=4, n_max-3
                d1f(i,r)= B(3)*( f(i+3, r) - f(i-3, r) )   &
                +   B(2)*( f(i+2,r) - f(i-2,r))   &
                +   B(1)*( f(i+1,r) - f(i-1,r) )
            enddo

            !call ddttrsb( 'N', n_max-2, 1, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), d1f(2:n_max-1, r), n_max-2, info )

        enddo

        ! Plus optimisé si on fait une boucle pour les conditions aux limites
        do r=1, nrhs

            d1f(2,r)= B1(-1)*f(1,r) + B1(0) * f(2,r)          &
            + B1(1) * f(3,r) + B1(2) * f(4,r)        &
            + B1(3) * f(5,r)

            d1f(3, r)=     B2(1)*( f(4, r) - f(2, r) )  &
            +   B2(2)*( f(5, r) - f(1, r) )

            d1f(n_max-2, r)= B2(1)*( f(n_max-1, r) - f(n_max-3, r) )  &
                + B2(2)*( f(n_max, r)   - f(n_max-4, r) )

            d1f(n_max-1,r)= -( B1(-1)*f(n_max,r) + B1(0)*f(n_max-1,r)     &
            + B1(1)*f(n_max-2,r)+ B1(2)*f(n_max-3,r)      &
            + B1(3)*f(n_max-4,r) )


            !call ddttrsb( 'N', n_max-2, 1, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), d1f(2:n_max-1, r), n_max-2, info )

        enddo
        ! Resolution du systeme [A] d1f = d
        !call DGTTS2( 0, n_max-2, 1, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), du2(1:n_max-4), IPIV, d1f(2:n_max-1), n_max-2 )
        call ddttrsb( 'N', n_max-2, nrhs, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), d1f(2:n_max-1, 1:nrhs), n_max-2, info )
            !call DGTTRS( 'N', n_max-2, 1, DL(1:n_max-3), D(1:n_max-2), DU(1:n_max-3), du2(1:n_max-4), ipiv, d1f(2:n_max-1), n_max-2, info )
            !call Lapack_solver(DL, D, DU, d1f(2:n_max-1), d1f(2:n_max-1), n_max-2)


        return

    end subroutine

#endif

    subroutine mDRP7(f, d1f, n_max, h, shifted, bc_type, nrhs)

        implicit none

        integer, intent(in)     :: n_max, bc_type, nrhs
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:,:)
        real(kind=8), intent(out)     :: d1f(:,:)

        integer :: i, r
        real(kind=8), allocatable, save :: A(:)
        integer, save   :: last_nmax=0

        if (last_nmax/=n_max) then

            if (allocated(A)) deallocate(A)
            allocate(A(7))
            A(1) =  0.95396219562045d0    /h
            A(2) =  -0.41234590494721d0   /h
            A(3) =  0.21233981563217d0    /h
            A(4) =  -0.10672135533957d0   /h
            A(5) =  0.046775295199319d0   /h
            A(6) = -0.015323662211638d0   /h
            A(7) = 0.0026664396358809d0   /h
        end if

        last_nmax=n_max

        !write(*,*)'LOULOULOU'

        do r = 1, nrhs



            do i=8,n_max-8
                d1f(i, r)= A(7)*(f(i+7,r) - f(i-7,r))                  &
                + A(6)*(f(i+6, r) - f(i-6, r))     &
                + A(5)*(f(i+5, r) - f(i-5, r))     &
                + A(4)*(f(i+4, r) - f(i-4, r))     &
                + A(3)*(f(i+3, r) - f(i-3, r))     &
                + A(2)*(f(i+2, r) - f(i-2, r))     &
                + A(1)*(f(i+1, r) - f(i-1, r))
            enddo


            if (bc_type.eq.Dirichlet) then
                call D1c_CTR_CLOSURE_INF(f(:,r), d1f(:,r), h, 7)
                call D1c_CTR_CLOSURE_SUP(f(:,r), d1f(:,r), n_max, h, n_max-7, shifted)
            end if

        end do

        return

    end subroutine

    subroutine mDRP7_2(f, d1f, n_max, h, shifted, bc_type, nrhs)

        implicit none

        integer, intent(in)     :: n_max, bc_type, nrhs
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:,:)
        real(kind=8), intent(out)     :: d1f(:,:)

        integer :: i, r
        real*8  :: A1(2), A2(3), A3(4), A4(5), A5(6)
        real(kind=8), allocatable, save :: A(:)
        integer, save   :: last_nmax=0

        if (last_nmax/=n_max) then

            if (allocated(A)) deallocate(A)
            allocate(A(7))
            A(1) =  0.95396219562045d0    /h
            A(2) =  -0.41234590494721d0   /h
            A(3) =  0.21233981563217d0    /h
            A(4) =  -0.10672135533957d0   /h
            A(5) =  0.046775295199319d0   /h
            A(6) = -0.015323662211638d0   /h
            A(7) = 0.0026664396358809d0   /h
        end if

        last_nmax=n_max

        !write(*,*)'LOULOULOU'

        do r = 1, nrhs


            do i=8,n_max-8
                d1f(i, r)= A(7)*(f(i+7,r) - f(i-7,r))                  &
                + A(6)*(f(i+6, r) - f(i-6, r))     &
                + A(5)*(f(i+5, r) - f(i-5, r))     &
                + A(4)*(f(i+4, r) - f(i-4, r))     &
                + A(3)*(f(i+3, r) - f(i-3, r))     &
                + A(2)*(f(i+2, r) - f(i-2, r))     &
                + A(1)*(f(i+1, r) - f(i-1, r))
            enddo


        end do

        do r = 1, nrhs

            call D1c_CTR_CLOSURE_INF(f(:,r), d1f(:,r), h, 7)
            call D1c_CTR_CLOSURE_SUP(f(:,r), d1f(:,r), n_max, h, n_max-7, shifted)

        end do

        return

    end subroutine

    subroutine mDRP6(f, d1f, n_max, h, shifted, bc_type, nrhs)

        implicit none

        integer, intent(in)     :: n_max, bc_type, nrhs
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:,:)
        real(kind=8), intent(out)     :: d1f(:,:)

        integer :: i, r
        real*8  :: A1(2), A2(3), A3(4), A4(5), A5(6)
        real(kind=8), allocatable, save :: A(:)
        integer, save   :: last_nmax=0

        if (last_nmax/=n_max) then

            if (allocated(A)) deallocate(A)
            allocate(A(6))
            A(1) =  0.91595792650492d0    /h
            A(2) =  -0.34922251636223d0   /h
            A(3) =  0.14398145036906d0    /h
            A(4) =  -0.051236991729043d0  /h
            A(5) =  0.013273181125903d0   /h
            A(6) = -0.0018126562894445d0  /h
        end if

        last_nmax=n_max

        !write(*,*)'LOULOULOU'

        do r = 1, nrhs


            do i=7, n_max-7
                d1f(i, r)=   A(6)*(f(i+6, r) - f(i-6, r))    &
                + A(5)*(f(i+5, r) - f(i-5, r))    &
                + A(4)*(f(i+4, r) - f(i-4, r))    &
                + A(3)*(f(i+3, r) - f(i-3, r))    &
                + A(2)*(f(i+2, r) - f(i-2, r))    &
                + A(1)*(f(i+1, r) - f(i-1, r))
            enddo


        end do

        return

    end subroutine

    subroutine mTamm(f, d1f, n_max, h, shifted, bc_type, nrhs)

        implicit none

        integer, intent(in)     :: n_max, bc_type, nrhs
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:,:)
        real(kind=8), intent(out)     :: d1f(:,:)

        integer :: i, r
        real*8  :: A1(2), A2(3), A3(4), A4(5), A5(6)
        real(kind=8), allocatable, save :: A(:)
        integer, save   :: last_nmax=0

        if (last_nmax/=n_max) then

            if (allocated(A)) deallocate(A)
            allocate(A(3))

            A(1) =  0.79926643d0 / h
            A(2) =  -0.18941314d0 / h
            A(3) =  0.02651995d0  / h
        end if

        last_nmax=n_max

        !write(*,*)'LOULOULOU'

        do r = 1, nrhs


            do i=4, n_max-4
                d1f(i, r)=   A(3)*(f(i+3, r) - f(i-3, r))    &
                + A(2)*(f(i+2, r) - f(i-2, r))    &
                + A(1)*(f(i+1, r) - f(i-1, r))
            enddo


        end do

        do r = 1, nrhs

            call D1c_CTR_CLOSURE_INF(f(:,r), d1f(:,r), h, 3)
            call D1c_CTR_CLOSURE_SUP(f(:,r), d1f(:,r), n_max, h, n_max-3, shifted)

        end do

        return

    end subroutine



    subroutine mDRP7_blas(f, d1f, n_max, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n_max, bc_type
        logical, intent(in)     :: shifted
        real(kind=8), intent(in)      :: h, f(:)
        real(kind=8), intent(out)     :: d1f(:)

        integer :: j
        integer, save   :: last_nmax=0

        ! BLAS data
        integer, parameter :: ml=7, mu=7
        real*8, dimension(:,:), allocatable, save       :: A, matrix

        integer, save     :: N, M
        real*8  :: alpha=1.d0, beta=0.d0
        real*8  :: a1, a2, a3, a4, a5, a6, a7
        real*8  :: b11, b21, b22, b31, b32, b33, b41, b42, b43, b44, b51, b52, b53, b54, b55
        real*8  :: b61, b62, b63, b64, b65, b66


        !call example1
        !call example2
        !call example3
        !call example4
        !call example5
        !call example6
        !call exit




        if (last_nmax/=n_max) then

            if (allocated(A)) deallocate(A, matrix)
            allocate(A(ml+mu+1, n_max))
            allocate(matrix(n_max, n_max))

            N=n_max
            M=n_max

            b11=0.5d0 /h; b21=2.d0/(3.d0) /h; b22=-1.d0/(12.d0) /h
            b31=0.79926643d0 /h; b32=-0.18941314d0 /h; b33=0.02651995d0 /h
            b51=0.88149995153015d0 /h; b52=-0.29786830893263d0 /h; b53= 0.098184076271814d0 /h
            b54=-0.02388959221513d0 /h; b55=0.0030486998975861d0 /h
            b61= 0.91595792650492d0/h; b62=- 0.34922251636223d0/h; b63= 0.14398145036906d0/h
            b64=- 0.051236991729043d0/h; b65= 0.013273181125903d0/h; b66=- 0.0018126562894445d0/h;


            a1= 0.95396219562045d0/h; a2=- 0.41234590494721d0/h; a3= 0.21233981563217d0/h
            a4=- 0.10672135533957d0/h; a5= 0.046775295199319d0/h; a6=- 0.015323662211638d0/h; a7= 0.0026664396358809d0/h;

            matrix=0.d0

            j=2
            matrix(j,j+1)=b11
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!

            j=3
            matrix(j,j+1)=b21
            matrix(j,j+2)=b22
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!

            j=4
            matrix(j,j+1)=b31
            matrix(j,j+2)=b32
            matrix(j,j+3)=b33
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!

            j=5
            matrix(j,j+1)=b31
            matrix(j,j+2)=b32
            matrix(j,j+3)=b33
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!

            j=6
            matrix(j,j+1)=b51
            matrix(j,j+2)=b52
            matrix(j,j+3)=b53
            matrix(j,j+4)=b54
            matrix(j,j+5)=b55
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-4)=-matrix(j,j+4) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-5)=-matrix(j,j+5) !!!!!!!!!!!!!!!!!!!

            j=7
            matrix(j,j+1)=b61
            matrix(j,j+2)=b62
            matrix(j,j+3)=b63
            matrix(j,j+4)=b64
            matrix(j,j+5)=b65
            matrix(j,j+6)=b66
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-4)=-matrix(j,j+4) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-5)=-matrix(j,j+5) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-6)=-matrix(j,j+6) !!!!!!!!!!!!!!!!!!!


            do j = 8, M-7
                matrix(j,j+1)=a1
                matrix(j,j+2)=a2
                matrix(j,j+3)=a3
                matrix(j,j+4)=a4
                matrix(j,j+5)=a5
                matrix(j,j+6)=a6
                matrix(j,j+7)=a7
                matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
                matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
                matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!
                matrix(j,j-4)=-matrix(j,j+4) !!!!!!!!!!!!!!!!!!!
                matrix(j,j-5)=-matrix(j,j+5) !!!!!!!!!!!!!!!!!!!
                matrix(j,j-6)=-matrix(j,j+6) !!!!!!!!!!!!!!!!!!!
                matrix(j,j-7)=-matrix(j,j+7) !!!!!!!!!!!!!!!!!!!

            end do

            j=M-6
            matrix(j,j+1)=b61
            matrix(j,j+2)=b62
            matrix(j,j+3)=b63
            matrix(j,j+4)=b64
            matrix(j,j+5)=b65
            matrix(j,j+6)=b66
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-4)=-matrix(j,j+4) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-5)=-matrix(j,j+5) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-6)=-matrix(j,j+6) !!!!!!!!!!!!!!!!!!!

            j=M-5
            matrix(j,j+1)=b51
            matrix(j,j+2)=b52
            matrix(j,j+3)=b53
            matrix(j,j+4)=b54
            matrix(j,j+5)=b55
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-4)=-matrix(j,j+4) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-5)=-matrix(j,j+5) !!!!!!!!!!!!!!!!!!!

            j=M-4
            matrix(j,j+1)=b31
            matrix(j,j+2)=b32
            matrix(j,j+3)=b33
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!

            j=M-3
            matrix(j,j+1)=b31
            matrix(j,j+2)=b32
            matrix(j,j+3)=b33
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-3)=-matrix(j,j+3) !!!!!!!!!!!!!!!!!!!

            j=M-2
            matrix(j,j+1)=b21
            matrix(j,j+2)=b22
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!
            matrix(j,j-2)=-matrix(j,j+2) !!!!!!!!!!!!!!!!!!!

            j=M-1
            matrix(j,j+1)=b11
            matrix(j,j-1)=-matrix(j,j+1) !!!!!!!!!!!!!!!!!!!

            j=M
            matrix(j,j)=0.d0

            call create_bs_matrix(A, matrix, ml, mu, N, M)

        end if

        last_nmax=n_max

        !f=1.d0

        !CALL DGBMV( 'N' , M , N , ML , MU ,  alpha , A , ml+mu+1 , f , 1  , beta , d1f , 1)
        CALL DGBMV( 'N' , N , N , ML , MU ,  1.d0 , A , ml+mu+1 , f , 1  , 0.d0 , d1f , 1)


        return

    contains

        subroutine create_bs_matrix(A, matrix, kl, ku, N, M)
            implicit none
            integer :: kl, ku, N, M
            real*8, dimension(kl+ku+1, N)       :: A
            real*8, dimension(M, N)       :: matrix
            integer     :: i, j, k

            A=0.d0
            do J = 1, N
                K = KU + 1 - J
                do I = MAX( 1, J - KU ), MIN( M, J + KL )
                    A( K + I, J ) = matrix( I, J )
                enddo
            enddo

        end subroutine create_bs_matrix

    end subroutine

    subroutine D1c_CTR_CLOSURE_INF(f, d1f, h, imax)

        ! Values location

        ! |------x------o------x------o----......---o------x---
        ! fs    df1     f1    df2     f2              df(imax)

        implicit none

        integer :: i, n_max, imax
        real(kind=8) h
        real(kind=8) d1f(:)
        real(kind=8) f(:)


        ! i=2---------------------------------------------------
        i=2
        d1f(i)= (f(i+1)-f(i-1)) /   (2.d0*h)

        if (i.eq.imax) return
        ! i=3---------------------------------------------------
        i=3
        d1f(i)    =   2.d0*(f(i+1)-f(i-1))  / (3.d0*h)      &
        -   (f(i+2) - f(i-2))         / (12.d0*h)

        if (i.eq.imax) return

        ! i=4---------------------------------------------------
        i=4
        d1f(i)  =   0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h

        if (i.eq.imax) return

        ! i=5---------------------------------------------------
        i=5
        d1f(i)  =   0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h

        if (i.eq.imax) return

        ! i=6---SEE D1_EXP_5pts.wxm-----------------------------
        i=6
        d1f(i)  =   0.88149995153015d0*(f(i+1)-f(i-1))      / h     &
        -   0.29786830893263d0*(f(i+2)-f(i-2))              / h     &
        +   0.098184076271814d0* (f(i+3) - f(i-3))          / h     &
        -   0.02388959221513d0* (f(i+4) - f(i-4))           / h     &
        +   0.0030486998975861d0* (f(i+5) - f(i-5))         / h

        if (i.eq.imax) return

        ! i=7---SEE D1_EXP_6pts.wxm-----------------------------
        i=7
        d1f(i)  =   0.91595792650492d0  * (f(i+1)-f(i-1))   / h     &
        -   0.34922251636223d0  * (f(i+2)-f(i-2))           / h     &
        +   0.14398145036906d0  * (f(i+3) - f(i-3))         / h     &
        -   0.051236991729043d0 * (f(i+4) - f(i-4))         / h     &
        +   0.013273181125903d0 * (f(i+5) - f(i-5))         / h     &
        -   0.0018126562894445d0* (f(i+6) - f(i-6))         / h

        return

    end subroutine

    subroutine D1c_CTR_CLOSURE_SUP(f, d1f, n_max, h, imin, shifted)

        implicit none

        integer :: i, n_max, imin, j
        real(kind=8) h, fn
        real(kind=8) d1f(:)
        real(kind=8) f(:)
        logical :: shifted

        if (shifted) then
            i=n_max-2
        else
            i=n_max-1
        end if

        ! First point-------------------------------------------
        d1f(i)= (f(i+1)-f(i-1)) /   (2.d0*h)

        if (i.eq.imin) return
        ! 2nd point----------------------------------------------
        i=i-1
        d1f(i)    =   2.d0*(f(i+1)-f(i-1))  / (3.d0*h)      &
        -   (f(i+2) - f(i-2))         / (12.d0*h)

        if (i.eq.imin) return


        ! 3rd point---------------------------------------------
        i=i-1
        d1f(i)  =   0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h

        if (i.eq.imin) return

        ! 4th point---------------------------------------------
        i=i-1
        d1f(i)  =   0.79926643d0*(f(i+1)-f(i-1))    / h     &
        -   0.18941314d0*(f(i+2)-f(i-2))            / h     &
        +   0.02651995d0* (f(i+3) - f(i-3))         / h

        if (i.eq.imin)return

        ! 5th point ----D1_EXP_5pts.wxm------------------------
        i=i-1
        d1f(i)  =   0.88149995153015d0*(f(i+1)-f(i-1))      / h     &
        -   0.29786830893263d0*(f(i+2)-f(i-2))              / h     &
        +   0.098184076271814d0* (f(i+3) - f(i-3))          / h     &
        -   0.02388959221513d0* (f(i+4) - f(i-4))           / h     &
        +   0.0030486998975861d0* (f(i+5) - f(i-5))         / h

        if (i.eq.imin) return

        ! 6th point ----D1_EXP_6pts.wxm------------------------
        i=i-1
        d1f(i)  =   0.91595792650492d0  * (f(i+1)-f(i-1))   / h     &
        -   0.34922251636223d0  * (f(i+2)-f(i-2))           / h     &
        +   0.14398145036906d0  * (f(i+3) - f(i-3))         / h     &
        -   0.051236991729043d0 * (f(i+4) - f(i-4))         / h     &
        +   0.013273181125903d0 * (f(i+5) - f(i-5))         / h     &
        -   0.0018126562894445d0* (f(i+6) - f(i-6))         / h

        if (i.eq.imin) return

        ! 7th point ---------------------------------------------
        i=i-1
        d1f(i)  =   0.95396219562045d0  * (f(i+1)-f(i-1))   / h     &
        -   0.41234590494721d0  * (f(i+2)-f(i-2))           / h     &
        +   0.21233981563217d0  * (f(i+3) - f(i-3))         / h     &
        -   0.10672135533957d0 * (f(i+4) - f(i-4))          / h     &
        +   0.046775295199319d0 * (f(i+5) - f(i-5))         / h     &
        -   0.015323662211638d0* (f(i+6) - f(i-6))          / h     &
        +   0.0026664396358809d0* (f(i+7) - f(i-7))          / h


        return

    end subroutine

    subroutine TS_init_loc(A, B, C, n_max)

        ! a - sub-diagonal (means it is the diagonal below the main diagonal)
        ! b - the main diagonal
        ! c - sup-diagonal (means it is the diagonal above the main diagonal)
        ! v - right part
        ! x - the answer
        ! n_max - number of equations

        implicit none
        integer,intent(in)                      :: n_max
        real(8),dimension(n_max),intent(in)     :: A, B, C
        integer i

        if (allocated(m_loc)) then
            deallocate(m_loc, bp_loc, upper_loc)
        end if

        allocate(m_loc(n_max), bp_loc(n_max), upper_loc(n_max))

        bp_loc(1) = 1.d0/B(1)
        upper_loc(1)=C(1)

        !   The first pass (setting coefficients):
        do i = 2,n_max
            m_loc(i) = A(i)*bp_loc(i-1)
            bp_loc(i) = 1.d0/(B(i) - m_loc(i)*C(i-1))
            upper_loc(i)=C(i)
        enddo

    end subroutine

    subroutine fill_matrix_opt(n_max, bc_type, shifted, h)
        implicit none
        integer :: i, n_max, bc_type
        logical :: shifted
        real(kind=8)  :: h

        if (allocated(A)) then
            deallocate(A)
        end if

        allocate(A(-1:1, n_max))

        do i=3, n_max-3
            A(-1,i) = 0.3793894912d0
            A(0,i) = 1.d0
            A(1,i) = 0.3793894912d0
        enddo


        A(-1,2)=0.d0
        A(0,2)=1.d0
        A(1,2)=0.d0

        B(1)=1.57557379d0       / (2.d0*h)
        B(2)= 0.1832051925d0    / (4.d0*h)

        Bs(-1)=-0.35d0             / h
        Bs(0)=-0.43333333333333d0  / h
        Bs(1)=0.9d0                / h
        Bs(2)=-0.1d0               / h
        Bs(3)=-0.016666666666667d0 / h

        ! End coefficients -----------------------
        Be(-1) =   -0.27333333333333d0 / h
        Be(0) =    -0.74d0             / h
        Be(1) =    1.36d0              / h
        Be(2) =    -0.40666666666667d0 / h
        Be(3) =    0.06d0              / h

        A(-1,n_max-2) = 0.3793894912d0
        A(0,n_max-2) = 1.d0
        A(1,n_max-2) = 0.3793894912d0

        A(-1,n_max-1)=0.d0
        A(0,n_max-1)=1.d0
        A(1,n_max-1)=0.d0

        call TS_init_loc(A(-1,2:n_max-1), A(0,2:n_max-1), A(1,2:n_max-1), n_max-2)


    end subroutine fill_matrix_opt

end module optimized_schemes_multiple_RHS
