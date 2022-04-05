module DRP_D1c
    use boundaries_types
    implicit none

contains

    subroutine D1c_DRP6_3Dx(f, d1f, n1,n2,n2e,n3,n3e, h, shifted, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: d1f

        integer :: i,j,k
        real(kind=8) A(7)

        A(1) =  0.91567612963343d0    /h
        A(2) =  -0.34876133244745d0   /h
        A(3) =  0.14348458980167d0    /h
        A(4) =  -0.050850207601385d0  /h
        A(5) =  0.013051374066223d0   /h
        A(6) = -0.0017438790115182d0  /h

        do k=1,n3e
            do j=1,n2e
                do i=7,n1-7
                    d1f(i, j, k)= A(6)*(f(i+6, j, k) - f(i-6, j, k))     &
                    + A(5)*(f(i+5, j, k) - f(i-5, j, k))     &
                    + A(4)*(f(i+4, j, k) - f(i-4, j, k))     &
                    + A(3)*(f(i+3, j, k) - f(i-3, j, k))     &
                    + A(2)*(f(i+2, j, k) - f(i-2, j, k))     &
                    + A(1)*(f(i+1, j, k) - f(i-1, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do k=1,n3e
                do j=1,n2e

                    do i=1, 6
                        d1f(i, j, k)=   A(6)*(f(i+6, j, k) - f(mod(n1+i-8,n1-1)+1, j, k))  &
                        + A(5)*(f(i+5, j, k) - f(mod(n1+i-7,n1-1)+1, j, k))  &
                        + A(4)*(f(i+4, j, k) - f(mod(n1+i-6,n1-1)+1, j, k))  &
                        + A(3)*(f(i+3, j, k) - f(mod(n1+i-5,n1-1)+1, j, k))  &
                        + A(2)*(f(i+2, j, k) - f(mod(n1+i-4,n1-1)+1, j, k))  &
                        + A(1)*(f(i+1, j, k) - f(mod(n1+i-3,n1-1)+1, j, k))
                    enddo

                    do i=n1-6, n1-1
                        d1f(i, j, k)= A(6)*(f(mod(i+5,n1-1)+1, j, k) - f(i-6, j, k))             &
                        + A(5)*(f(mod(i+4,n1-1)+1, j, k) - f(i-5, j, k))             &
                        + A(4)*(f(mod(i+3,n1-1)+1, j, k) - f(i-4, j, k))             &
                        + A(3)*(f(mod(i+2,n1-1)+1, j, k) - f(i-3, j, k))             &
                        + A(2)*(f(mod(i+1,n1-1)+1, j, k) - f(i-2, j, k))             &
                        + A(1)*(f(mod(i,n1-1)+1, j, k) - f(i-1, j, k))
                    enddo
                enddo

            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1,n3e
                do j=1,n2e

                    d1f(2,j,k)= 0.5d0*(f(3,j,k)-f(1,j,k)) /   h

                    ! i=3---------------------------------------------------
                    d1f(3,j,k)    =   2.d0*(f(4,j,k)-f(2,j,k))  / (3.d0*h)      &
                    -   (f(5,j,k) - f(1,j,k))         / (12.d0*h)

                    ! i=4---------------------------------------------------
                    d1f(4,j,k)  =   0.79926643d0*(f(5,j,k)-f(3,j,k))    / h     &
                    -   0.18941314d0*(f(6,j,k)-f(2,j,k))            / h     &
                    +   0.02651995d0* (f(7,j,k) - f(1,j,k))         / h

                    ! i=5---------------------------------------------------
                    d1f(5,j,k)  =   0.79926643d0*(f(6,j,k)-f(4,j,k))    / h     &
                    -   0.18941314d0*(f(7,j,k)-f(3,j,k))            / h     &
                    +   0.02651995d0* (f(8,j,k) - f(2,j,k))         / h


                    ! i=6---SEE D1_EXP_5pts.wxm-----------------------------
                    d1f(6,j,k)  =   0.88084033745666d0*(f(7,j,k)-f(5,j,k))      / h     &
                    -   0.29689284187948d0*(f(8,j,k)-f(4,j,k))              / h     &
                    +   0.097346322640388d0* (f(9,j,k) - f(3,j,k))          / h     &
                    -   0.023424784076677d0* (f(10,j,k) - f(2,j,k))           / h     &
                    +   0.0029211029375669d0* (f(11,j,k) - f(1,j,k))         / h

                enddo
            enddo



            if (shifted) then
                i=n1-2
            else
                i=n1-1
            end if

            do k=1,n3e
                do j=1,n2e
                    d1f(i,j,k)= (f(i+1,j,k)-f(i-1,j,k)) /   (2.d0*h)
                enddo
            enddo

            ! 2nd point----------------------------------------------
            i=i-1

            do k=1,n3e
                do j=1,n2e
                    d1f(i,j,k)    =   2.d0*(f(i+1,j,k)-f(i-1,j,k))  / (3.d0*h)      &
                    -   (f(i+2,j,k) - f(i-2,j,k))         / (12.d0*h)
                enddo
            enddo


            ! 3rd point---------------------------------------------
            i=i-1
            do k=1,n3e
                do j=1,n2e
                    d1f(i,j,k)  =   0.79926643d0*(f(i+1,j,k)-f(i-1,j,k))    / h     &
                    -   0.18941314d0*(f(i+2,j,k)-f(i-2,j,k))            / h     &
                    +   0.02651995d0* (f(i+3,j,k) - f(i-3,j,k))         / h
                enddo
            enddo

            ! 4th point---------------------------------------------
            i=i-1
            do k=1,n3e
                do j=1,n2e
                    d1f(i,j,k)  =   0.79926643d0*(f(i+1,j,k)-f(i-1,j,k))    / h     &
                    -   0.18941314d0*(f(i+2,j,k)-f(i-2,j,k))            / h     &
                    +   0.02651995d0* (f(i+3,j,k) - f(i-3,j,k))         / h
                enddo
            enddo

            ! 5th point ----D1_EXP_5pts.wxm------------------------
            i=i-1
            do k=1,n3e
                do j=1,n2e
                    d1f(i,j,k)  =   0.88084033745666d0*(f(i+1,j,k)-f(i-1,j,k))      / h     &
                    -   0.29689284187948d0*(f(i+2,j,k)-f(i-2,j,k))              / h     &
                    +   0.097346322640388d0* (f(i+3,j,k) - f(i-3,j,k))          / h     &
                    -   0.023424784076677d0* (f(i+4,j,k) - f(i-4,j,k))           / h     &
                    +   0.0029211029375669d0* (f(i+5,j,k) - f(i-5,j,k))         / h
                enddo
            enddo

            ! 6th point ----D1_EXP_6pts.wxm------------------------
            i=i-1
            do k=1,n3e
                do j=1,n2e
                    d1f(i,j,k)  =   A(1)  * (f(i+1,j,k)-f(i-1,j,k))     &
                    +   A(2)  * (f(i+2,j,k)-f(i-2,j,k))                 &
                    +   A(3)  * (f(i+3,j,k) - f(i-3,j,k))               &
                    +   A(4) * (f(i+4,j,k) - f(i-4,j,k))                &
                    +   A(5) * (f(i+5,j,k) - f(i-5,j,k))                &
                    +   A(6)* (f(i+6,j,k) - f(i-6,j,k))

                enddo
            enddo


        endif

        return

    end subroutine D1c_DRP6_3Dx



    subroutine D1c_DRP6_3Dy(f, d1f, n1,n1e,n2,n3,n3e, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: d1f

        integer :: i,j,k
        real(kind=8) A(7)

        A(1) =  0.91567612963343d0    /h
        A(2) =  -0.34876133244745d0   /h
        A(3) =  0.14348458980167d0    /h
        A(4) =  -0.050850207601385d0  /h
        A(5) =  0.013051374066223d0   /h
        A(6) = -0.0017438790115182d0  /h

        do k=1, n3e
            do i=1, n1e

                do j=7,n2-7
                    d1f(i,j,k)= A(6)*(f(i,j+6,k) - f(i,j-6,k))     &
                    + A(5)*(f(i,j+5,k) - f(i,j-5,k))     &
                    + A(4)*(f(i,j+4,k) - f(i,j-4,k))     &
                    + A(3)*(f(i,j+3,k) - f(i,j-3,k))     &
                    + A(2)*(f(i,j+2,k) - f(i,j-2,k))     &
                    + A(1)*(f(i,j+1,k) - f(i,j-1,k))
                enddo

            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=1, 6
                        d1f(i,j,k)=   A(6)*(f(i,j+6,k) - f(i,mod(n2+j-8,n2-1)+1,k))  &
                        + A(5)*(f(i,j+5,k) - f(i,mod(n2+j-7,n2-1)+1,k))  &
                        + A(4)*(f(i,j+4,k) - f(i,mod(n2+j-6,n2-1)+1,k))  &
                        + A(3)*(f(i,j+3,k) - f(i,mod(n2+j-5,n2-1)+1,k))  &
                        + A(2)*(f(i,j+2,k) - f(i,mod(n2+j-4,n2-1)+1,k))  &
                        + A(1)*(f(i,j+1,k) - f(i,mod(n2+j-3,n2-1)+1,k))
                    enddo

                    do j=n2-6, n2-1
                        d1f(i,j,k)= A(6)*(f(i,mod(j+5,n2-1)+1,k) - f(i,j-6,k))             &
                        + A(5)*(f(i,mod(j+4,n2-1)+1,k) - f(i,j-5,k))             &
                        + A(4)*(f(i,mod(j+3,n2-1)+1,k) - f(i,j-4,k))             &
                        + A(3)*(f(i,mod(j+2,n2-1)+1,k) - f(i,j-3,k))             &
                        + A(2)*(f(i,mod(j+1,n2-1)+1,k) - f(i,j-2,k))             &
                        + A(1)*(f(i,mod(j,n2-1)+1,k) - f(i,j-1,k))
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e

                    d1f(i,2,k)= 0.5d0*(f(i,3,k)-f(i,1,k)) /   h
                    d1f(i,3,k)    =   2.d0*(f(i,4,k)-f(i,2,k))  / (3.d0*h)      &
                    -   (f(i,5,k) - f(i,1,k))         / (12.d0*h)
                    d1f(i,4,k)  =   0.79926643d0*(f(i,5,k)-f(i,3,k))    / h     &
                    -   0.18941314d0*(f(i,6,k)-f(i,2,k))            / h     &
                    +   0.02651995d0* (f(i,7,k) - f(i,1,k))         / h
                    d1f(i,5,k)  =   0.79926643d0*(f(i,6,k)-f(i,4,k))    / h     &
                    -   0.18941314d0*(f(i,7,k)-f(i,3,k))            / h     &
                    +   0.02651995d0* (f(i,8,k) - f(i,2,k))         / h
                    d1f(i,6,k)  =   0.88084033745666d0*(f(i,7,k)-f(i,5,k))      / h     &
                    -   0.29689284187948d0*(f(i,8,k)-f(i,4,k))              / h     &
                    +   0.097346322640388d0* (f(i,9,k) - f(i,3,k))          / h     &
                    -   0.023424784076677d0* (f(i,10,k) - f(i,2,k))           / h     &
                    +   0.0029211029375669d0* (f(i,11,k) - f(i,1,k))         / h

                enddo
            enddo



            if (shifted) then

                do k=1, n3e
                    do i=1, n1e
                        d1f(i,n2-2,k)= (f(i,n2-1,k)-f(i,n2-3,k))                            / (2.d0*h)
                        d1f(i,n2-3,k)    =     2.d0*(f(i,n2-2,k)-f(i,n2-4,k))               / (3.d0*h)      &
                        -   (f(i,n2-1,k) - f(i,n2-5,k))                                     / (12.d0*h)
                        d1f(i,n2-4,k)  =   0.79926643d0*(f(i,n2-3,k)-f(i,n2-5,k))           / h     &
                        -   0.18941314d0*(f(i,n2-2,k)-f(i,n2-6,k))                          / h     &
                        +   0.02651995d0* (f(i,n2-1,k) - f(i,n2-7,k))                       / h
                        d1f(i,n2-5,k)  =   0.79926643d0*(f(i,n2-4,k)-f(i,n2-6,k))           / h     &
                        -   0.18941314d0*(f(i,n2-3,k)-f(i,n2-7,k))                          / h     &
                        +   0.02651995d0* (f(i,n2-2,k) - f(i,n2-8,k))                       / h
                        d1f(i,n2-6,k)  =   0.88084033745666d0*(f(i,n2-5,k)-f(i,n2-7,k))     / h     &
                        -   0.29689284187948d0*(f(i,n2-4,k)-f(i,n2-8,k))                    / h     &
                        +   0.097346322640388d0* (f(i,n2-3,k) - f(i,n2-9,k))                / h     &
                        -   0.023424784076677d0* (f(i,n2-2,k) - f(i,n2-10,k))                / h     &
                        +   0.0029211029375669d0* (f(i,n2-1,k) - f(i,n2-11,k))              / h
                    enddo
                enddo

            else

                do k=1, n3e
                    do i=1, n1e
                        d1f(i,n2-1,k)= (f(i,n2,k)-f(i,n2-2,k)) /   (2.d0*h)
                        d1f(i,n2-2,k)    =     2.d0*(f(i,n2-1,k)-f(i,n2-3,k))               / (3.d0*h)      &
                        -   (f(i,n2,k) - f(i,n2-4,k))                                       / (12.d0*h)
                        d1f(i,n2-3,k)  =   0.79926643d0*(f(i,n2-2,k)-f(i,n2-4,k))           / h     &
                        -   0.18941314d0*(f(i,n2-1,k)-f(i,n2-5,k))                          / h     &
                        +   0.02651995d0* (f(i,n2,k) - f(i,n2-6,k))                         / h
                        d1f(i,n2-4,k)  =   0.79926643d0*(f(i,n2-3,k)-f(i,n2-5,k))           / h     &
                        -   0.18941314d0*(f(i,n2-2,k)-f(i,n2-6,k))                          / h     &
                        +   0.02651995d0* (f(i,n2-1,k) - f(i,n2-7,k))                       / h
                        d1f(i,n2-5,k)  =   0.88084033745666d0*(f(i,n2-4,k)-f(i,n2-6,k))     / h     &
                        -   0.29689284187948d0*(f(i,n2-3,k)-f(i,n2-7,k))                    / h     &
                        +   0.097346322640388d0* (f(i,n2-2,k) - f(i,n2-8,k))                / h     &
                        -   0.023424784076677d0* (f(i,n2-1,k) - f(i,n2-9,k))                 / h     &
                        +   0.0029211029375669d0* (f(i,n2,k) - f(i,n2-10,k))                / h
                        d1f(i,n2-6,k)  =   A(1)  * (f(i,n2-5,k)-f(i,n2-7,k))    &
                        +   A(2)  * (f(i,n2-4,k)-f(i,n2-8,k))                   &
                        +   A(3)  * (f(i,n2-3,k) - f(i,n2-9,k))                 &
                        +   A(4) * (f(i,n2-2,k) - f(i,n2-10,k))                 &
                        +   A(5) * (f(i,n2-1,k) - f(i,n2-11,k))                 &
                        +   A(6)* (f(i,n2,k) - f(i,n2-12,k))
                    enddo
                enddo

            end if

        endif

        return

    end subroutine D1c_DRP6_3Dy



    subroutine D1c_DRP6_MULT_3Dy(f, d1f,  n1,n1e,n2,n3,n3e, h, shifted, bc_type, g)

        implicit none

        integer, intent(in)     ::  n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: d1f
        real(kind=8), dimension(:), intent(in)                       :: g

        integer :: i,j,k
        real(kind=8) A(7)

        A(1) =  0.91567612963343d0    /h
        A(2) =  -0.34876133244745d0   /h
        A(3) =  0.14348458980167d0    /h
        A(4) =  -0.050850207601385d0  /h
        A(5) =  0.013051374066223d0   /h
        A(6) = -0.0017438790115182d0  /h

        do k=1, n3e
            do i=1, n1e

                do j=7,n2-7
                    d1f(i,j,k)= (A(6)*(f(i,j+6,k) - f(i,j-6,k))     &
                    + A(5)*(f(i,j+5,k) - f(i,j-5,k))     &
                    + A(4)*(f(i,j+4,k) - f(i,j-4,k))     &
                    + A(3)*(f(i,j+3,k) - f(i,j-3,k))     &
                    + A(2)*(f(i,j+2,k) - f(i,j-2,k))     &
                    + A(1)*(f(i,j+1,k) - f(i,j-1,k)))*g(j)
                enddo

            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=1, 6
                        d1f(i,j,k)=   (A(6)*(f(i,j+6,k) - f(i,mod(n2+j-8,n2-1)+1,k))  &
                        + A(5)*(f(i,j+5,k) - f(i,mod(n2+j-7,n2-1)+1,k))  &
                        + A(4)*(f(i,j+4,k) - f(i,mod(n2+j-6,n2-1)+1,k))  &
                        + A(3)*(f(i,j+3,k) - f(i,mod(n2+j-5,n2-1)+1,k))  &
                        + A(2)*(f(i,j+2,k) - f(i,mod(n2+j-4,n2-1)+1,k))  &
                        + A(1)*(f(i,j+1,k) - f(i,mod(n2+j-3,n2-1)+1,k)))*g(j)
                    enddo

                    do j=n2-6, n2-1
                        d1f(i,j,k)= (A(6)*(f(i,mod(j+5,n2-1)+1,k) - f(i,j-6,k))             &
                        + A(5)*(f(i,mod(j+4,n2-1)+1,k) - f(i,j-5,k))             &
                        + A(4)*(f(i,mod(j+3,n2-1)+1,k) - f(i,j-4,k))             &
                        + A(3)*(f(i,mod(j+2,n2-1)+1,k) - f(i,j-3,k))             &
                        + A(2)*(f(i,mod(j+1,n2-1)+1,k) - f(i,j-2,k))             &
                        + A(1)*(f(i,mod(j,n2-1)+1,k) - f(i,j-1,k)))*g(j)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e

                    d1f(i,2,k)= (0.5d0*(f(i,3,k)-f(i,1,k)) /   h)*g(2)

                    ! i=3---------------------------------------------------
                    d1f(i,3,k)    =   (2.d0*(f(i,4,k)-f(i,2,k))  / (3.d0*h)      &
                    -   (f(i,5,k) - f(i,1,k))         / (12.d0*h))*g(3)

                    ! i=4---------------------------------------------------
                    d1f(i,4,k)  =   (0.79926643d0*(f(i,5,k)-f(i,3,k))    / h     &
                    -   0.18941314d0*(f(i,6,k)-f(i,2,k))            / h     &
                    +   0.02651995d0* (f(i,7,k) - f(i,1,k))         / h)*g(4)

                    ! i=5---------------------------------------------------
                    d1f(i,5,k)  =   (0.79926643d0*(f(i,6,k)-f(i,4,k))    / h     &
                    -   0.18941314d0*(f(i,7,k)-f(i,3,k))            / h     &
                    +   0.02651995d0* (f(i,8,k) - f(i,2,k))         / h)*g(5)

                    ! i=6---SEE D1_EXP_5pts.wxm-----------------------------
                    d1f(i,6,k)  =   (0.88084033745666d0*(f(i,7,k)-f(i,5,k))      / h     &
                    -   0.29689284187948d0*(f(i,8,k)-f(i,4,k))              / h     &
                    +   0.097346322640388d0* (f(i,9,k) - f(i,3,k))          / h     &
                    -   0.023424784076677d0* (f(i,10,k) - f(i,2,k))           / h     &
                    +   0.0029211029375669d0* (f(i,11,k) - f(i,1,k))         / h)*g(6)

                enddo
            enddo



            if (shifted) then

                do k=1, n3e
                    do i=1, n1e
                        d1f(i,n2-2,k)= ((f(i,n2-1,k)-f(i,n2-3,k))                           / (2.d0*h))*g(n2-2)
                        d1f(i,n2-3,k)    =     (2.d0*(f(i,n2-2,k)-f(i,n2-4,k))               / (3.d0*h)      &
                        -   (f(i,n2-1,k) - f(i,n2-5,k))                                     / (12.d0*h))*g(n2-3)
                        d1f(i,n2-4,k)  =   (0.79926643d0*(f(i,n2-3,k)-f(i,n2-5,k))          / h     &
                        -   0.18941314d0*(f(i,n2-2,k)-f(i,n2-6,k))                          / h     &
                        +   0.02651995d0* (f(i,n2-1,k) - f(i,n2-7,k))                       / h)*g(n2-4)
                        d1f(i,n2-5,k)  =   (0.79926643d0*(f(i,n2-4,k)-f(i,n2-6,k))          / h     &
                        -   0.18941314d0*(f(i,n2-3,k)-f(i,n2-7,k))                          / h     &
                        +   0.02651995d0* (f(i,n2-2,k) - f(i,n2-8,k))                       / h)*g(n2-5)
                        d1f(i,n2-6,k)  =   (0.88084033745666d0*(f(i,n2-5,k)-f(i,n2-7,k))    / h     &
                        -   0.29689284187948d0*(f(i,n2-4,k)-f(i,n2-8,k))                    / h     &
                        +   0.097346322640388d0* (f(i,n2-3,k) - f(i,n2-9,k))                / h     &
                        -   0.023424784076677d0* (f(i,n2-2,k) - f(i,n2-10,k))                / h     &
                        +   0.0029211029375669d0* (f(i,n2-1,k) - f(i,n2-11,k))              / h)*g(n2-6)
                    enddo
                enddo

            else

                do k=1, n3e
                    do i=1, n1e
                        d1f(i,n2-1,k)= ((f(i,n2,k)-f(i,n2-2,k))                             / (2.d0*h))*g(n2-1)
                        d1f(i,n2-2,k)    =     (2.d0*(f(i,n2-1,k)-f(i,n2-3,k))               / (3.d0*h)      &
                        -   (f(i,n2,k) - f(i,n2-4,k))                                       / (12.d0*h))*g(n2-2)
                        d1f(i,n2-3,k)  =   (0.79926643d0*(f(i,n2-2,k)-f(i,n2-4,k))          / h     &
                        -   0.18941314d0*(f(i,n2-1,k)-f(i,n2-5,k))                          / h     &
                        +   0.02651995d0* (f(i,n2,k) - f(i,n2-6,k))                         / h)*g(n2-3)
                        d1f(i,n2-4,k)  =   (0.79926643d0*(f(i,n2-3,k)-f(i,n2-5,k))          / h     &
                        -   0.18941314d0*(f(i,n2-2,k)-f(i,n2-6,k))                          / h     &
                        +   0.02651995d0* (f(i,n2-1,k) - f(i,n2-7,k))                       / h)*g(n2-4)
                        d1f(i,n2-5,k)  =   (0.88084033745666d0*(f(i,n2-4,k)-f(i,n2-6,k))    / h     &
                        -   0.29689284187948d0*(f(i,n2-3,k)-f(i,n2-7,k))                    / h     &
                        +   0.097346322640388d0* (f(i,n2-2,k) - f(i,n2-8,k))                / h     &
                        -   0.023424784076677d0* (f(i,n2-1,k) - f(i,n2-9,k))                 / h     &
                        +   0.0029211029375669d0* (f(i,n2,k) - f(i,n2-10,k))              / h)*g(n2-5)
                        d1f(i,n2-6,k)  =   (A(1)  * (f(i,n2-5,k)-f(i,n2-7,k))   &
                        +   A(2)  * (f(i,n2-4,k)-f(i,n2-8,k))                   &
                        +   A(3)  * (f(i,n2-3,k) - f(i,n2-9,k))                 &
                        +   A(4) * (f(i,n2-2,k) - f(i,n2-10,k))                 &
                        +   A(5) * (f(i,n2-1,k) - f(i,n2-11,k))                 &
                        +   A(6)* (f(i,n2,k) - f(i,n2-12,k))                    )*g(n2-6)
                    enddo
                enddo

            end if

        endif

        return

    end subroutine D1c_DRP6_MULT_3Dy



    subroutine D1c_DRP6_MULTACC_3Dy(f, d1f, n1,n1e,n2,n3,n3e, h, shifted, bc_type, g)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: d1f
        real(kind=8), dimension(:), intent(in)                       :: g

        integer :: i,j,k
        real(kind=8) A(7)

        A(1) =  0.91567612963343d0    /h
        A(2) =  -0.34876133244745d0   /h
        A(3) =  0.14348458980167d0    /h
        A(4) =  -0.050850207601385d0  /h
        A(5) =  0.013051374066223d0   /h
        A(6) = -0.0017438790115182d0  /h

        do k=1, n3e
            do i=1, n1e

                do j=7,n2-7
                    d1f(i,j,k)= d1f(i,j,k)+(A(6)*(f(i,j+6,k) - f(i,j-6,k))     &
                    + A(5)*(f(i,j+5,k) - f(i,j-5,k))     &
                    + A(4)*(f(i,j+4,k) - f(i,j-4,k))     &
                    + A(3)*(f(i,j+3,k) - f(i,j-3,k))     &
                    + A(2)*(f(i,j+2,k) - f(i,j-2,k))     &
                    + A(1)*(f(i,j+1,k) - f(i,j-1,k)))*g(j)
                enddo

            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=1, 6
                        d1f(i,j,k)=   d1f(i,j,k)+(A(6)*(f(i,j+6,k) - f(i,mod(n2+j-8,n2-1)+1,k))  &
                        + A(5)*(f(i,j+5,k) - f(i,mod(n2+j-7,n2-1)+1,k))  &
                        + A(4)*(f(i,j+4,k) - f(i,mod(n2+j-6,n2-1)+1,k))  &
                        + A(3)*(f(i,j+3,k) - f(i,mod(n2+j-5,n2-1)+1,k))  &
                        + A(2)*(f(i,j+2,k) - f(i,mod(n2+j-4,n2-1)+1,k))  &
                        + A(1)*(f(i,j+1,k) - f(i,mod(n2+j-3,n2-1)+1,k)))*g(j)
                    enddo

                    do j=n2-6, n2-1
                        d1f(i,j,k)= d1f(i,j,k)+(A(6)*(f(i,mod(j+5,n2-1)+1,k) - f(i,j-6,k))             &
                        + A(5)*(f(i,mod(j+4,n2-1)+1,k) - f(i,j-5,k))             &
                        + A(4)*(f(i,mod(j+3,n2-1)+1,k) - f(i,j-4,k))             &
                        + A(3)*(f(i,mod(j+2,n2-1)+1,k) - f(i,j-3,k))             &
                        + A(2)*(f(i,mod(j+1,n2-1)+1,k) - f(i,j-2,k))             &
                        + A(1)*(f(i,mod(j,n2-1)+1,k) - f(i,j-1,k)))*g(j)
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e

                    d1f(i,2,k)= d1f(i,2,k)+(0.5d0*(f(i,3,k)-f(i,1,k)) /   h)*g(2)

                    ! i=3---------------------------------------------------
                    d1f(i,3,k)    =   d1f(i,3,k)+(2.d0*(f(i,4,k)-f(i,2,k))  / (3.d0*h)      &
                    -   (f(i,5,k) - f(i,1,k))         / (12.d0*h))*g(3)

                    ! i=4---------------------------------------------------
                    d1f(i,4,k)  =   d1f(i,4,k)+(0.79926643d0*(f(i,5,k)-f(i,3,k))    / h     &
                    -   0.18941314d0*(f(i,6,k)-f(i,2,k))            / h     &
                    +   0.02651995d0* (f(i,7,k) - f(i,1,k))         / h)*g(4)

                    ! i=5---------------------------------------------------
                    d1f(i,5,k)  =   d1f(i,5,k)+(0.79926643d0*(f(i,6,k)-f(i,4,k))    / h     &
                    -   0.18941314d0*(f(i,7,k)-f(i,3,k))            / h     &
                    +   0.02651995d0* (f(i,8,k) - f(i,2,k))         / h)*g(5)

                    ! i=6---SEE D1_EXP_5pts.wxm-----------------------------
                    d1f(i,6,k)  =   d1f(i,6,k)+(0.88084033745666d0*(f(i,7,k)-f(i,5,k))      / h     &
                    -   0.29689284187948d0*(f(i,8,k)-f(i,4,k))              / h     &
                    +   0.097346322640388d0* (f(i,9,k) - f(i,3,k))          / h     &
                    -   0.023424784076677d0* (f(i,10,k) - f(i,2,k))           / h     &
                    +   0.0029211029375669d0* (f(i,11,k) - f(i,1,k))         / h)*g(6)

                enddo
            enddo



            if (shifted) then

                do k=1, n3e
                    do i=1, n1e
                        d1f(i,n2-2,k)= d1f(i,n2-2,k)+((f(i,n2-1,k)-f(i,n2-3,k))                             / (2.d0*h))*g(n2-2)
                        d1f(i,n2-3,k)    =     d1f(i,n2-3,k)+(2.d0*(f(i,n2-2,k)-f(i,n2-4,k))                / (3.d0*h)      &
                        -   (f(i,n2-1,k) - f(i,n2-5,k))                                                     / (12.d0*h))*g(n2-3)
                        d1f(i,n2-4,k)  =   d1f(i,n2-4,k)+(0.79926643d0*(f(i,n2-3,k)-f(i,n2-5,k))            / h     &
                        -   0.18941314d0*(f(i,n2-2,k)-f(i,n2-6,k))                                          / h     &
                        +   0.02651995d0* (f(i,n2-1,k) - f(i,n2-7,k))                                       / h)*g(n2-4)
                        d1f(i,n2-5,k)  =   d1f(i,n2-5,k)+(0.79926643d0*(f(i,n2-4,k)-f(i,n2-6,k))            / h     &
                        -   0.18941314d0*(f(i,n2-3,k)-f(i,n2-7,k))                                          / h     &
                        +   0.02651995d0* (f(i,n2-2,k) - f(i,n2-8,k))                                       / h)*g(n2-5)
                        d1f(i,n2-6,k)  =   d1f(i,n2-6,k)+(0.88084033745666d0*(f(i,n2-5,k)-f(i,n2-7,k))      / h     &
                        -   0.29689284187948d0*(f(i,n2-4,k)-f(i,n2-8,k))                                    / h     &
                        +   0.097346322640388d0* (f(i,n2-3,k) - f(i,n2-9,k))                                / h     &
                        -   0.023424784076677d0* (f(i,n2-2,k) - f(i,n2-10,k))                                / h     &
                        +   0.0029211029375669d0* (f(i,n2-1,k) - f(i,n2-11,k))                              / h)*g(n2-6)
                    enddo
                enddo

            else

                do k=1, n3e
                    do i=1, n1e
                        d1f(i,n2-1,k)= d1f(i,n2-1,k)+((f(i,n2,k)-f(i,n2-2,k))                               / (2.d0*h))*g(n2-1)
                        d1f(i,n2-2,k)    =     d1f(i,n2-2,k)+(2.d0*(f(i,n2-1,k)-f(i,n2-3,k))                / (3.d0*h)      &
                        -   (f(i,n2,k) - f(i,n2-4,k))                                                       / (12.d0*h))*g(n2-2)
                        d1f(i,n2-3,k)  =   d1f(i,n2-3,k)+(0.79926643d0*(f(i,n2-2,k)-f(i,n2-4,k))            / h     &
                        -   0.18941314d0*(f(i,n2-1,k)-f(i,n2-5,k))                                          / h     &
                        +   0.02651995d0* (f(i,n2,k) - f(i,n2-6,k))                                         / h)*g(n2-3)
                        d1f(i,n2-4,k)  =   d1f(i,n2-4,k)+(0.79926643d0*(f(i,n2-3,k)-f(i,n2-5,k))            / h     &
                        -   0.18941314d0*(f(i,n2-2,k)-f(i,n2-6,k))                                          / h     &
                        +   0.02651995d0* (f(i,n2-1,k) - f(i,n2-7,k))                                       / h)*g(n2-4)
                        d1f(i,n2-5,k)  =   d1f(i,n2-5,k)+(0.88084033745666d0*(f(i,n2-4,k)-f(i,n2-6,k))      / h     &
                        -   0.29689284187948d0*(f(i,n2-3,k)-f(i,n2-7,k))                                    / h     &
                        +   0.097346322640388d0* (f(i,n2-2,k) - f(i,n2-8,k))                                / h     &
                        -   0.023424784076677d0* (f(i,n2-1,k) - f(i,n2-9,k))                                 / h     &
                        +   0.0029211029375669d0* (f(i,n2,k) - f(i,n2-10,k))                                / h)*g(n2-5)
                        d1f(i,n2-6,k)  =   d1f(i,n2-6,k)+(A(1)  * (f(i,n2-5,k)-f(i,n2-7,k))     &
                        +   A(2)  * (f(i,n2-4,k)-f(i,n2-8,k))                                   &
                        +   A(3)  * (f(i,n2-3,k) - f(i,n2-9,k))                                 &
                        +   A(4) * (f(i,n2-2,k) - f(i,n2-10,k))                                 &
                        +   A(5) * (f(i,n2-1,k) - f(i,n2-11,k))                                 &
                        +   A(6)* (f(i,n2,k) - f(i,n2-12,k))                                    )*g(n2-6)
                    enddo
                enddo

            end if

        endif

        return

    end subroutine D1c_DRP6_MULTACC_3Dy
    subroutine D1c_DRP6_3Dz(f, d1f, n1,n1e,n2,n2e,n3, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: d1f

        integer :: i,j,k
        real(kind=8) A(7)

        A(1) =  0.91567612963343d0    /h
        A(2) =  -0.34876133244745d0   /h
        A(3) =  0.14348458980167d0    /h
        A(4) =  -0.050850207601385d0  /h
        A(5) =  0.013051374066223d0   /h
        A(6) = -0.0017438790115182d0  /h

        do j=1, n2e
            do i=1, n1e
                do k=7,n3-7
                    d1f(i,j,k)= A(6)*(f(i,j,k+6) - f(i,j,k-6))     &
                    + A(5)*(f(i,j,k+5) - f(i,j,k-5))     &
                    + A(4)*(f(i,j,k+4) - f(i,j,k-4))     &
                    + A(3)*(f(i,j,k+3) - f(i,j,k-3))     &
                    + A(2)*(f(i,j,k+2) - f(i,j,k-2))     &
                    + A(1)*(f(i,j,k+1) - f(i,j,k-1))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e

                    do k=1, 6
                        d1f(i,j,k)=   A(6)*(f(i,j,k+6) - f(i,j,mod(n3+k-8,n3-1)+1))  &
                        + A(5)*(f(i,j,k+5) - f(i,j,mod(n3+k-7,n3-1)+1))  &
                        + A(4)*(f(i,j,k+4) - f(i,j,mod(n3+k-6,n3-1)+1))  &
                        + A(3)*(f(i,j,k+3) - f(i,j,mod(n3+k-5,n3-1)+1))  &
                        + A(2)*(f(i,j,k+2) - f(i,j,mod(n3+k-4,n3-1)+1))  &
                        + A(1)*(f(i,j,k+1) - f(i,j,mod(n3+k-3,n3-1)+1))
                    enddo

                    do k=n3-6, n3-1
                        d1f(i,j,k)= A(6)*(f(i,j,mod(k+5,n3-1)+1) - f(i,j,k-6))             &
                        + A(5)*(f(i,j,mod(k+4,n3-1)+1) - f(i,j,k-5))             &
                        + A(4)*(f(i,j,mod(k+3,n3-1)+1) - f(i,j,k-4))             &
                        + A(3)*(f(i,j,mod(k+2,n3-1)+1) - f(i,j,k-3))             &
                        + A(2)*(f(i,j,mod(k+1,n3-1)+1) - f(i,j,k-2))             &
                        + A(1)*(f(i,j,mod(k,n3-1)+1) - f(i,j,k-1))
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e

                    d1f(i,j,2)= 0.5d0*(f(i,j,3)-f(i,j,1)) /   h

                    ! i=3---------------------------------------------------
                    d1f(i,j,3)    =   2.d0*(f(i,j,4)-f(i,j,2))  / (3.d0*h)      &
                    -   (f(i,j,5) - f(i,j,1))         / (12.d0*h)

                    ! i=4---------------------------------------------------
                    d1f(i,j,4)  =   0.79926643d0*(f(i,j,5)-f(i,j,3))    / h     &
                    -   0.18941314d0*(f(i,j,6)-f(i,j,2))            / h     &
                    +   0.02651995d0* (f(i,j,7) - f(i,j,1))         / h

                    ! i=5---------------------------------------------------
                    d1f(i,j,5)  =   0.79926643d0*(f(i,j,6)-f(i,j,4))    / h     &
                    -   0.18941314d0*(f(i,j,7)-f(i,j,3))            / h     &
                    +   0.02651995d0* (f(i,j,8) - f(i,j,2))         / h

                    ! i=6---SEE D1_EXP_5pts.wxm-----------------------------
                    d1f(i,j,6)  =   0.88084033745666d0*(f(i,j,7)-f(i,j,5))      / h     &
                    -   0.29689284187948d0*(f(i,j,8)-f(i,j,4))              / h     &
                    +   0.097346322640388d0* (f(i,j,9) - f(i,j,3))          / h     &
                    -   0.023424784076677d0* (f(i,j,10) - f(i,j,2))           / h     &
                    +   0.0029211029375669d0* (f(i,j,11) - f(i,j,1))         / h

                enddo
            enddo

            if (shifted) then
                k=n3-2
            else
                k=n3-1
            end if

            ! First point-------------------------------------------
            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,k)= (f(i,j,k+1)-f(i,j,k-1)) /   (2.d0*h)
                enddo
            enddo

            ! 2nd point----------------------------------------------
            k=k-1
            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,k)    =   2.d0*(f(i,j,k+1)-f(i,j,k-1))  / (3.d0*h)      &
                    -   (f(i,j,k+2) - f(i,j,k-2))         / (12.d0*h)
                enddo
            enddo

            ! 3rd point---------------------------------------------
            k=k-1
            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,k)  =   0.79926643d0*(f(i,j,k+1)-f(i,j,k-1))    / h     &
                    -   0.18941314d0*(f(i,j,k+2)-f(i,j,k-2))            / h     &
                    +   0.02651995d0* (f(i,j,k+3) - f(i,j,k-3))         / h
                enddo
            enddo

            ! 4th point---------------------------------------------
            k=k-1
            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,k)  =   0.79926643d0*(f(i,j,k+1)-f(i,j,k-1))    / h     &
                    -   0.18941314d0*(f(i,j,k+2)-f(i,j,k-2))            / h     &
                    +   0.02651995d0* (f(i,j,k+3) - f(i,j,k-3))         / h
                enddo
            enddo

            ! 5th point ----D1_EXP_5pts.wxm------------------------
            k=k-1
            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,k)  =   0.88084033745666d0*(f(i,j,k+1)-f(i,j,k-1))      / h     &
                    -   0.29689284187948d0*(f(i,j,k+2)-f(i,j,k-2))              / h     &
                    +   0.097346322640388d0* (f(i,j,k+3) - f(i,j,k-3))          / h     &
                    -   0.023424784076677d0* (f(i,j,k+4) - f(i,j,k-4))           / h     &
                    +   0.0029211029375669d0* (f(i,j,k+5) - f(i,j,k-5))         / h
                enddo
            enddo

            ! 6th point ----D1_EXP_6pts.wxm------------------------
            k=k-1
            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,k)  =   A(1)  * (f(i,j,k+1)-f(i,j,k-1))     &
                    +   A(2)  * (f(i,j,k+2)-f(i,j,k-2))                 &
                    +   A(3)  * (f(i,j,k+3) - f(i,j,k-3))               &
                    +   A(4) * (f(i,j,k+4) - f(i,j,k-4))                &
                    +   A(5) * (f(i,j,k+5) - f(i,j,k-5))                &
                    +   A(6)* (f(i,j,k+6) - f(i,j,k-6))
                enddo
            enddo

        endif

        return

    end subroutine D1c_DRP6_3Dz

end module DRP_D1c
