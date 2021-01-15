module DRP_D1s
    use boundaries_types
    implicit none

contains

    subroutine D1s_DRP5_3Dx(f, d1f, n1,n2,n2e,n3,n3e, h, shifted, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: d1f

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        do k=1,n3e
            do j=1,n2e
                do i=5,n1-6
                    d1f(i, j, k)= A(5)*(f(i+5, j, k) - f(i-4, j, k))     &
                    + A(4)*(f(i+4, j, k) - f(i-3, j, k))     &
                    + A(3)*(f(i+3, j, k) - f(i-2, j, k))     &
                    + A(2)*(f(i+2, j, k) - f(i-1, j, k))     &
                    + A(1)*(f(i+1, j, k) - f(i, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do k=1,n3e
                do j=1,n2e

                    do i=1, 4
                        d1f(i, j, k)=   A(5)*(f(i+5, j, k) - f(mod(n1+i-6,n1-1)+1, j, k))  &
                        + A(4)*(f(i+4, j, k) - f(mod(n1+i-5,n1-1)+1, j, k))  &
                        + A(3)*(f(i+3, j, k) - f(mod(n1+i-4,n1-1)+1, j, k))  &
                        + A(2)*(f(i+2, j, k) - f(mod(n1+i-3,n1-1)+1, j, k))  &
                        + A(1)*(f(i+1, j, k) - f(i, j, k))
                    enddo

                    do i=n1-5,n1-1
                        d1f(i, j, k)= A(5)*(f(mod(i+4,n1-1)+1, j, k) - f(i-4, j, k))             &
                        + A(4)*(f(mod(i+3,n1-1)+1, j, k) - f(i-3, j, k))             &
                        + A(3)*(f(mod(i+2,n1-1)+1, j, k) - f(i-2, j, k))             &
                        + A(2)*(f(mod(i+1,n1-1)+1, j, k) - f(i-1, j, k))             &
                        + A(1)*(f(mod(i,n1-1)+1, j, k) - f(i, j, k))
                    enddo
                enddo

            enddo

        endif

        if (bc_type.eq.antisymetric) then

            do k=1, n3e
                do j=1, n2e

                    d1f(1,j,k)=   A(5)*(f(6,j,k) + f(5,j,k))   &
                    + A(4)*(f(5,j,k) + f(4,j,k))           &
                    + A(3)*(f(4,j,k) + f(3,j,k))           &
                    + A(2)*(f(3,j,k) + f(2,j,k))           &
                    + A(1)*(f(2,j,k) - f(1,j,k))

                    d1f(1,j,k)=d1f(1,j,k)-2.d0*(A(2)+A(3)+A(4)+A(5))*f(1,j,k)

                    d1f(2,j,k)=   A(5)*(f(7,j,k) + f(4,j,k))   &
                    + A(4)*(f(6,j,k) + f(3,j,k))           &
                    + A(3)*(f(5,j,k) + f(2,j,k))           &
                    + A(2)*(f(4,j,k) - f(1,j,k))           &
                    + A(1)*(f(3,j,k) - f(2,j,k))

                    d1f(2,j,k)=d1f(2,j,k)-2.d0*(A(3)+A(4)+A(5))*f(1,j,k)

                    d1f(3,j,k)=   A(5)*(f(8,j,k) + f(3,j,k))   &
                    + A(4)*(f(7,j,k) + f(2,j,k))           &
                    + A(3)*(f(6,j,k) - f(1,j,k))           &
                    + A(2)*(f(5,j,k) - f(2,j,k))           &
                    + A(1)*(f(4,j,k) - f(3,j,k))

                    d1f(3,j,k)=d1f(3,j,k)-2.d0*(A(4)+A(5))*f(1,j,k)

                    d1f(4,j,k)=   A(5)*(f(9,j,k) + f(2,j,k))   &
                    + A(4)*(f(8,j,k) - f(1,j,k))           &
                    + A(3)*(f(7,j,k) - f(2,j,k))           &
                    + A(2)*(f(6,j,k) - f(3,j,k))           &
                    + A(1)*(f(5,j,k) - f(4,j,k))

                    d1f(4,j,k)=d1f(4,j,k)-2.d0*(A(5))*f(1,j,k)


                    d1f(n1-1,j,k)=   A(5)*(-f(n1-4,j,k) - f(n1-5,j,k))     &
                    + A(4)*(-f(n1-3,j,k) - f(n1-4,j,k))             &
                    + A(3)*(-f(n1-2,j,k) - f(n1-3,j,k))             &
                    + A(2)*(-f(n1-1,j,k) - f(n1-2,j,k))             &
                    + A(1)*(f(n1,j,k) - f(n1-1,j,k))

                    d1f(n1-1,j,k)=d1f(n1-1,j,k)+2.d0*(A(2)+A(3)+A(4)+A(5))*f(n1,j,k)


                    d1f(n1-2,j,k)=   A(5)*(-f(n1-3,j,k) - f(n1-6,j,k))     &
                    + A(4)*(-f(n1-2,j,k) - f(n1-5,j,k))             &
                    + A(3)*(-f(n1-1,j,k) - f(n1-4,j,k))             &
                    + A(2)*(f(n1,j,k) - f(n1-3,j,k))             &
                    + A(1)*(f(n1-1,j,k) - f(n1-2,j,k))

                    d1f(n1-2,j,k)=d1f(n1-2,j,k)+2.d0*(A(3)+A(4)+A(5))*f(n1,j,k)


                    d1f(n1-3,j,k)=   A(5)*(-f(n1-2,j,k) - f(n1-7,j,k))     &
                    + A(4)*(-f(n1-1,j,k) - f(n1-6,j,k))             &
                    + A(3)*(f(n1,j,k) - f(n1-5,j,k))             &
                    + A(2)*(f(n1-1,j,k) - f(n1-4,j,k))             &
                    + A(1)*(f(n1-2,j,k) - f(n1-3,j,k))

                    d1f(n1-3,j,k)=d1f(n1-3,j,k)+2.d0*(A(4)+A(5))*f(n1,j,k)


                    d1f(n1-4,j,k)=   A(5)*(-f(n1-1,j,k) - f(n1-8,j,k))     &
                    + A(4)*(f(n1,j,k) - f(n1-7,j,k))             &
                    + A(3)*(f(n1-1,j,k) - f(n1-6,j,k))             &
                    + A(2)*(f(n1-2,j,k) - f(n1-5,j,k))             &
                    + A(1)*(f(n1-3,j,k) - f(n1-4,j,k))

                    d1f(n1-4,j,k)=d1f(n1-4,j,k)+2.d0*(A(5))*f(n1,j,k)


                    d1f(n1-5,j,k)=   A(5)*(f(n1,j,k) - f(n1-9,j,k))     &
                    + A(4)*(f(n1-1,j,k) - f(n1-8,j,k))             &
                    + A(3)*(f(n1-2,j,k) - f(n1-7,j,k))             &
                    + A(2)*(f(n1-3,j,k) - f(n1-6,j,k))             &
                    + A(1)*(f(n1-4,j,k) - f(n1-5,j,k))

                enddo
            enddo

        end if


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do j=1,n2e
                    d1f(1,j,k)= (f(2,j,k)-f(1,j,k))/h
                    d1f(2,j,k)  =  9.0d0/8.d0  * (f(3,j,k)-f(2,j,k))     / h     &
                    -   1.d0/24.d0  * (f(4,j,k)-f(1,j,k))           / h
                    d1f(3,j,k)  = 1.189101951200031d0  *   (f(4,j,k)-f(3,j,k))  / h         &
                    -   0.073717642266682d0  *   (f(5,j,k)-f(2,j,k))        / h         &
                    +   0.006410195120003d0 *   (f(6,j,k) - f(1,j,k))       / h
                    d1f(4,j,k)  = 1.189101951200031d0  *   (f(5,j,k)-f(4,j,k))  / h         &
                    -   0.073717642266682d0  *   (f(6,j,k)-f(3,j,k))        / h         &
                    +   0.006410195120003d0 *   (f(7,j,k) - f(2,j,k))       / h
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    d1f(n1-1,j,k)=(f(n1,j,k)-f(n1-1,j,k))/h
                    d1f(n1-2,j,k)  =   9.0d0/8.d0 * (f(n1-1,j,k)-f(n1-2,j,k))  /h    &
                    -   1.d0/24.d0    * (f(n1,j,k)-f(n1-3,j,k))   /h
                    d1f(n1-3,j,k)  =   1.189101951200031d0  *   (f(n1-2,j,k)-f(n1-3,j,k))   /h   &
                    -   0.073717642266682d0  *   (f(n1-1,j,k)-f(n1-4,j,k))                   /h  &
                    +   0.006410195120003d0 *   (f(n1,j,k) - f(n1-5,j,k))                   /h
                    d1f(n1-4,j,k)  =   1.189101951200031d0  *   (f(n1-3,j,k)-f(n1-4,j,k))    /h  &
                    -   0.073717642266682d0  *   (f(n1-2,j,k)-f(n1-5,j,k))                   /h  &
                    +   0.006410195120003d0 *   (f(n1-1,j,k) - f(n1-6,j,k))                 /h
                    d1f(n1-5,j,k)  =   A(1)  *   (f(n1-4,j,k)-f(n1-5,j,k))  &
                    +   A(2)  *   (f(n1-3,j,k)-f(n1-6,j,k))     &
                    +   A(3) *   (f(n1-2,j,k) - f(n1-7,j,k))    &
                    +   A(4) *   (f(n1-1,j,k) - f(n1-8,j,k))    &
                    +   A(5)*   (f(n1,j,k) - f(n1-9,j,k))
                enddo
            enddo


        endif

        return

    end subroutine D1s_DRP5_3Dx

    subroutine D1s_DRP5_ACC_3Dx(f, d1f, n1,n2,n2e,n3,n3e, h, shifted, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: d1f

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        do k=1,n3e
            do j=1,n2e
                do i=5,n1-6
                    d1f(i, j, k)= d1f(i, j, k)+A(5)*(f(i+5, j, k) - f(i-4, j, k))     &
                    + A(4)*(f(i+4, j, k) - f(i-3, j, k))     &
                    + A(3)*(f(i+3, j, k) - f(i-2, j, k))     &
                    + A(2)*(f(i+2, j, k) - f(i-1, j, k))     &
                    + A(1)*(f(i+1, j, k) - f(i, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do k=1,n3e
                do j=1,n2e

                    do i=1, 4
                        d1f(i, j, k)=   d1f(i, j, k)+A(5)*(f(i+5, j, k) - f(mod(n1+i-6,n1-1)+1, j, k))  &
                        + A(4)*(f(i+4, j, k) - f(mod(n1+i-5,n1-1)+1, j, k))  &
                        + A(3)*(f(i+3, j, k) - f(mod(n1+i-4,n1-1)+1, j, k))  &
                        + A(2)*(f(i+2, j, k) - f(mod(n1+i-3,n1-1)+1, j, k))  &
                        + A(1)*(f(i+1, j, k) - f(i, j, k))
                    enddo

                    do i=n1-5,n1-1
                        d1f(i, j, k)= d1f(i, j, k)+A(5)*(f(mod(i+4,n1-1)+1, j, k) - f(i-4, j, k))             &
                        + A(4)*(f(mod(i+3,n1-1)+1, j, k) - f(i-3, j, k))             &
                        + A(3)*(f(mod(i+2,n1-1)+1, j, k) - f(i-2, j, k))             &
                        + A(2)*(f(mod(i+1,n1-1)+1, j, k) - f(i-1, j, k))             &
                        + A(1)*(f(mod(i,n1-1)+1, j, k) - f(i, j, k))
                    enddo
                enddo

            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do j=1,n2e
                    d1f(1,j,k)= d1f(1,j,k)+(f(2,j,k)-f(1,j,k))/h
                    d1f(2,j,k)  =  d1f(2,j,k)+9.0d0/8.d0  * (f(3,j,k)-f(2,j,k))     / h     &
                    -   1.d0/24.d0  * (f(4,j,k)-f(1,j,k))           / h
                    d1f(3,j,k)  = d1f(3,j,k)+1.189101951200031d0  *   (f(4,j,k)-f(3,j,k))  / h         &
                    -   0.073717642266682d0  *   (f(5,j,k)-f(2,j,k))        / h         &
                    +   0.006410195120003d0 *   (f(6,j,k) - f(1,j,k))       / h
                    d1f(4,j,k)  = d1f(4,j,k)+1.189101951200031d0  *   (f(5,j,k)-f(4,j,k))  / h         &
                    -   0.073717642266682d0  *   (f(6,j,k)-f(3,j,k))        / h         &
                    +   0.006410195120003d0 *   (f(7,j,k) - f(2,j,k))       / h
                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    d1f(n1-1,j,k)=d1f(n1-1,j,k)+(f(n1,j,k)-f(n1-1,j,k))/h
                    d1f(n1-2,j,k)  =   d1f(n1-2,j,k)+9.0d0/8.d0 * (f(n1-1,j,k)-f(n1-2,j,k))  /h    &
                    -   1.d0/24.d0    * (f(n1,j,k)-f(n1-3,j,k))   /h
                    d1f(n1-3,j,k)  =   d1f(n1-3,j,k)+1.189101951200031d0  *   (f(n1-2,j,k)-f(n1-3,j,k))   /h   &
                    -   0.073717642266682d0  *   (f(n1-1,j,k)-f(n1-4,j,k))                   /h  &
                    +   0.006410195120003d0 *   (f(n1,j,k) - f(n1-5,j,k))                   /h
                    d1f(n1-4,j,k)  =   d1f(n1-4,j,k)+1.189101951200031d0  *   (f(n1-3,j,k)-f(n1-4,j,k))    /h  &
                    -   0.073717642266682d0  *   (f(n1-2,j,k)-f(n1-5,j,k))                   /h  &
                    +   0.006410195120003d0 *   (f(n1-1,j,k) - f(n1-6,j,k))                 /h
                    d1f(n1-5,j,k)  =  d1f(n1-5,j,k)+A(1)  *   (f(n1-4,j,k)-f(n1-5,j,k))     &
                    +   A(2)  *   (f(n1-3,j,k)-f(n1-6,j,k))     &
                    +   A(3) *   (f(n1-2,j,k) - f(n1-7,j,k))    &
                    +   A(4) *   (f(n1-1,j,k) - f(n1-8,j,k))    &
                    +   A(5)*   (f(n1,j,k) - f(n1-9,j,k))
                enddo
            enddo


        endif

        if ((bc_type.eq.antisymetric)) then
            write(*,*)'D1s_DRP5_ACC_3Dx: ANTISYMETRIC NOT SUPPORTED'
        endif

        return

    end subroutine D1s_DRP5_ACC_3Dx


    subroutine D1s_DRP5_3Dy(f, d1f, n1,n1e,n2,n3,n3e, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: d1f

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        do k=1, n3e
            do i=1, n1e
                do j=5,n2-6
                    d1f(i, j, k)= A(5)*(f(i, j+5, k) - f(i, j-4, k))     &
                    + A(4)*(f(i, j+4, k) - f(i, j-3, k))     &
                    + A(3)*(f(i, j+3, k) - f(i, j-2, k))     &
                    + A(2)*(f(i, j+2, k) - f(i, j-1, k))     &
                    + A(1)*(f(i, j+1, k) - f(i, j, k))
                enddo
            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=1, 4
                        d1f(i, j, k)=   A(5)*(f(i, j+5, k) - f(i, mod(n2+j-6,n2-1)+1, k))  &
                        + A(4)*(f(i, j+4, k) - f(i, mod(n2+j-5,n2-1)+1, k))  &
                        + A(3)*(f(i, j+3, k) - f(i, mod(n2+j-4,n2-1)+1, k))  &
                        + A(2)*(f(i, j+2, k) - f(i, mod(n2+j-3,n2-1)+1, k))  &
                        + A(1)*(f(i, j+1, k) - f(i, j, k))
                    enddo

                    do j=n2-5,n2-1
                        d1f(i, j, k)= A(5)*(f(i, mod(j+4,n2-1)+1, k) - f(i, j-4, k))             &
                        + A(4)*(f(i, mod(j+3,n2-1)+1, k) - f(i, j-3, k))             &
                        + A(3)*(f(i, mod(j+2,n2-1)+1, k) - f(i, j-2, k))             &
                        + A(2)*(f(i, mod(j+1,n2-1)+1, k) - f(i, j-1, k))             &
                        + A(1)*(f(i, mod(j,n2-1)+1, k) - f(i, j, k))
                    enddo

                enddo
            enddo

        endif

        if (bc_type.eq.antisymetric) then

            do k=1, n3e
                do i=1, n1e

                    d1f(i,1,k)=   A(5)*(f(i,6,k) + f(i,5,k))   &
                    + A(4)*(f(i,5,k) + f(i,4,k))           &
                    + A(3)*(f(i,4,k) + f(i,3,k))           &
                    + A(2)*(f(i,3,k) + f(i,2,k))           &
                    + A(1)*(f(i,2,k) - f(i,1,k))

                    d1f(i,1,k)=d1f(i,1,k)-2.d0*(A(2)+A(3)+A(4)+A(5))*f(i,1,k)

                    d1f(i,2,k)=   A(5)*(f(i,7,k) + f(i,4,k))   &
                    + A(4)*(f(i,6,k) + f(i,3,k))           &
                    + A(3)*(f(i,5,k) + f(i,2,k))           &
                    + A(2)*(f(i,4,k) - f(i,1,k))           &
                    + A(1)*(f(i,3,k) - f(i,2,k))

                    d1f(i,2,k)=d1f(i,2,k)-2.d0*(A(3)+A(4)+A(5))*f(i,1,k)

                    d1f(i,3,k)=   A(5)*(f(i,8,k) + f(i,3,k))   &
                    + A(4)*(f(i,7,k) + f(i,2,k))           &
                    + A(3)*(f(i,6,k) - f(i,1,k))           &
                    + A(2)*(f(i,5,k) - f(i,2,k))           &
                    + A(1)*(f(i,4,k) - f(i,3,k))

                    d1f(i,3,k)=d1f(i,3,k)-2.d0*(A(4)+A(5))*f(i,1,k)

                    d1f(i,4,k)=   A(5)*(f(i,9,k) + f(i,2,k))   &
                    + A(4)*(f(i,8,k) - f(i,1,k))           &
                    + A(3)*(f(i,7,k) - f(i,2,k))           &
                    + A(2)*(f(i,6,k) - f(i,3,k))           &
                    + A(1)*(f(i,5,k) - f(i,4,k))

                    d1f(i,4,k)=d1f(i,4,k)-2.d0*(A(5))*f(i,1,k)


                    d1f(i,n2-1,k)=   A(5)*(-f(i,n2-4,k) - f(i,n2-5,k))     &
                    + A(4)*(-f(i,n2-3,k) - f(i,n2-4,k))             &
                    + A(3)*(-f(i,n2-2,k) - f(i,n2-3,k))             &
                    + A(2)*(-f(i,n2-1,k) - f(i,n2-2,k))             &
                    + A(1)*(f(i,n2,k) - f(i,n2-1,k))

                    d1f(i,n2-1,k)=d1f(i,n2-1,k)+2.d0*(A(2)+A(3)+A(4)+A(5))*f(i,n2,k)


                    d1f(i,n2-2,k)=   A(5)*(-f(i,n2-3,k) - f(i,n2-6,k))     &
                    + A(4)*(-f(i,n2-2,k) - f(i,n2-5,k))             &
                    + A(3)*(-f(i,n2-1,k) - f(i,n2-4,k))             &
                    + A(2)*(f(i,n2,k) - f(i,n2-3,k))             &
                    + A(1)*(f(i,n2-1,k) - f(i,n2-2,k))

                    d1f(i,n2-2,k)=d1f(i,n2-2,k)+2.d0*(A(3)+A(4)+A(5))*f(i,n2,k)


                    d1f(i,n2-3,k)=   A(5)*(-f(i,n2-2,k) - f(i,n2-7,k))     &
                    + A(4)*(-f(i,n2-1,k) - f(i,n2-6,k))             &
                    + A(3)*(f(i,n2,k) - f(i,n2-5,k))             &
                    + A(2)*(f(i,n2-1,k) - f(i,n2-4,k))             &
                    + A(1)*(f(i,n2-2,k) - f(i,n2-3,k))

                    d1f(i,n2-3,k)=d1f(i,n2-3,k)+2.d0*(A(4)+A(5))*f(i,n2,k)


                    d1f(i,n2-4,k)=   A(5)*(-f(i,n2-1,k) - f(i,n2-8,k))     &
                    + A(4)*(f(i,n2,k) - f(i,n2-7,k))             &
                    + A(3)*(f(i,n2-1,k) - f(i,n2-6,k))             &
                    + A(2)*(f(i,n2-2,k) - f(i,n2-5,k))             &
                    + A(1)*(f(i,n2-3,k) - f(i,n2-4,k))

                    d1f(i,n2-4,k)=d1f(i,n2-4,k)+2.d0*(A(5))*f(i,n2,k)


                    d1f(i,n2-5,k)=   A(5)*(f(i,n2,k) - f(i,n2-9,k))     &
                    + A(4)*(f(i,n2-1,k) - f(i,n2-8,k))             &
                    + A(3)*(f(i,n2-2,k) - f(i,n2-7,k))             &
                    + A(2)*(f(i,n2-3,k) - f(i,n2-6,k))             &
                    + A(1)*(f(i,n2-4,k) - f(i,n2-5,k))

                enddo
            enddo

        end if


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e
                    d1f(i,1,k)= (f(i,2,k)-f(i,1,k))/h
                    d1f(i,2,k)  =  9.0d0/8.d0  * (f(i,3,k)-f(i,2,k))            / h     &
                    -   1.d0/24.d0  * (f(i,4,k)-f(i,1,k))                       / h
                    d1f(i,3,k)  = 1.189101951200031d0  *   (f(i,4,k)-f(i,3,k))  / h     &
                    -   0.073717642266682d0  *   (f(i,5,k)-f(i,2,k))            / h     &
                    +   0.006410195120003d0 *   (f(i,6,k) - f(i,1,k))           / h
                    d1f(i,4,k)  = 1.189101951200031d0  *   (f(i,5,k)-f(i,4,k))  / h     &
                    -   0.073717642266682d0  *   (f(i,6,k)-f(i,3,k))            / h     &
                    +   0.006410195120003d0 *   (f(i,7,k) - f(i,2,k))           / h

                enddo
            enddo

            do k=1, n3e
                do i=1, n1e
                    d1f(i,n2-1,k)=(f(i,n2,k)-f(i,n2-1,k))/h
                    d1f(i,n2-2,k)  =   9.0d0/8.d0 * (f(i,n2-1,k)-f(i,n2-2,k))               /h  &
                    -   1.d0/24.d0    * (f(i,n2,k)-f(i,n2-3,k))   /h
                    d1f(i,n2-3,k)  =   1.189101951200031d0  *   (f(i,n2-2,k)-f(i,n2-3,k))   /h  &
                    -   0.073717642266682d0  *   (f(i,n2-1,k)-f(i,n2-4,k))                  /h  &
                    +   0.006410195120003d0 *   (f(i,n2,k) - f(i,n2-5,k))                   /h
                    d1f(i,n2-4,k)  =   1.189101951200031d0  *   (f(i,n2-3,k)-f(i,n2-4,k))   /h  &
                    -   0.073717642266682d0  *   (f(i,n2-2,k)-f(i,n2-5,k))                  /h  &
                    +   0.006410195120003d0 *   (f(i,n2-1,k) - f(i,n2-6,k))                 /h
                    d1f(i,n2-5,k)  =   A(1)  *   (f(i,n2-4,k)-f(i,n2-5,k))      &
                    +   A(2)  *   (f(i,n2-3,k)-f(i,n2-6,k))     &
                    +   A(3) *   (f(i,n2-2,k) - f(i,n2-7,k))    &
                    +   A(4) *   (f(i,n2-1,k) - f(i,n2-8,k))    &
                    +   A(5)*   (f(i,n2,k) - f(i,n2-9,k))

                enddo
            enddo

        endif

        return

    end subroutine D1s_DRP5_3Dy


    subroutine D1s_DRP5_MULT_3Dy(f, d1f, n1,n1e,n2,n3,n3e, h, shifted, bc_type, g)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)       :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)      :: d1f
        real(kind=8), dimension(:), intent(in)              :: g

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        do k=1, n3e
            do i=1, n1e
                do j=5,n2-6
                    d1f(i, j, k)= (A(5)*(f(i, j+5, k) - f(i, j-4, k))     &
                    + A(4)*(f(i, j+4, k) - f(i, j-3, k))     &
                    + A(3)*(f(i, j+3, k) - f(i, j-2, k))     &
                    + A(2)*(f(i, j+2, k) - f(i, j-1, k))     &
                    + A(1)*(f(i, j+1, k) - f(i, j, k)))*g(j)
                enddo
            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=1, 4
                        d1f(i, j, k)=   (A(5)*(f(i, j+5, k) - f(i, mod(n2+j-6,n2-1)+1, k))  &
                        + A(4)*(f(i, j+4, k) - f(i, mod(n2+j-5,n2-1)+1, k))  &
                        + A(3)*(f(i, j+3, k) - f(i, mod(n2+j-4,n2-1)+1, k))  &
                        + A(2)*(f(i, j+2, k) - f(i, mod(n2+j-3,n2-1)+1, k))  &
                        + A(1)*(f(i, j+1, k) - f(i, j, k)))*g(j)
                    enddo

                    do j=n2-5,n2-1
                        d1f(i, j, k)= (A(5)*(f(i, mod(j+4,n2-1)+1, k) - f(i, j-4, k))             &
                        + A(4)*(f(i, mod(j+3,n2-1)+1, k) - f(i, j-3, k))             &
                        + A(3)*(f(i, mod(j+2,n2-1)+1, k) - f(i, j-2, k))             &
                        + A(2)*(f(i, mod(j+1,n2-1)+1, k) - f(i, j-1, k))             &
                        + A(1)*(f(i, mod(j,n2-1)+1, k) - f(i, j, k)))*g(j)
                    enddo

                enddo
            enddo

        endif

        if (bc_type.eq.antisymetric) then

            do k=1, n3e
                do i=1, n1e

                    d1f(i,1,k)=  (A(5)*(f(i,6,k) + f(i,5,k))   &
                    + A(4)*(f(i,5,k) + f(i,4,k))           &
                    + A(3)*(f(i,4,k) + f(i,3,k))           &
                    + A(2)*(f(i,3,k) + f(i,2,k))           &
                    + A(1)*(f(i,2,k) - f(i,1,k)) ) * g(1)

                    d1f(i,1,k)=d1f(i,1,k)-2.d0*(A(2)+A(3)+A(4)+A(5))*f(i,1,k)* g(1)

                    d1f(i,2,k)=  (A(5)*(f(i,7,k) + f(i,4,k))   &
                    + A(4)*(f(i,6,k) + f(i,3,k))           &
                    + A(3)*(f(i,5,k) + f(i,2,k))           &
                    + A(2)*(f(i,4,k) - f(i,1,k))           &
                    + A(1)*(f(i,3,k) - f(i,2,k)) ) * g(2)

                    d1f(i,2,k)=d1f(i,2,k)-2.d0*(A(3)+A(4)+A(5))*f(i,1,k)* g(2)

                    d1f(i,3,k)=  (A(5)*(f(i,8,k) + f(i,3,k))   &
                    + A(4)*(f(i,7,k) + f(i,2,k))           &
                    + A(3)*(f(i,6,k) - f(i,1,k))           &
                    + A(2)*(f(i,5,k) - f(i,2,k))           &
                    + A(1)*(f(i,4,k) - f(i,3,k)) ) * g(3)

                    d1f(i,3,k)=d1f(i,3,k)-2.d0*(A(4)+A(5))*f(i,1,k)* g(3)

                    d1f(i,4,k)=  (A(5)*(f(i,9,k) + f(i,2,k))   &
                    + A(4)*(f(i,8,k) - f(i,1,k))           &
                    + A(3)*(f(i,7,k) - f(i,2,k))           &
                    + A(2)*(f(i,6,k) - f(i,3,k))           &
                    + A(1)*(f(i,5,k) - f(i,4,k)) ) * g(4)

                    d1f(i,4,k)=d1f(i,4,k)-2.d0*(A(5))*f(i,1,k)* g(4)


                    d1f(i,n2-1,k)=  (A(5)*(-f(i,n2-4,k) - f(i,n2-5,k))     &
                    + A(4)*(-f(i,n2-3,k) - f(i,n2-4,k))             &
                    + A(3)*(-f(i,n2-2,k) - f(i,n2-3,k))             &
                    + A(2)*(-f(i,n2-1,k) - f(i,n2-2,k))             &
                    + A(1)*(f(i,n2,k) - f(i,n2-1,k)) ) * g(n2-1)

                    d1f(i,n2-1,k)=d1f(i,n2-1,k)+2.d0*(A(2)+A(3)+A(4)+A(5))*f(i,n2,k)* g(n2-1)


                    d1f(i,n2-2,k)=  (A(5)*(-f(i,n2-3,k) - f(i,n2-6,k))     &
                    + A(4)*(-f(i,n2-2,k) - f(i,n2-5,k))             &
                    + A(3)*(-f(i,n2-1,k) - f(i,n2-4,k))             &
                    + A(2)*(f(i,n2,k) - f(i,n2-3,k))             &
                    + A(1)*(f(i,n2-1,k) - f(i,n2-2,k)) ) * g(n2-2)

                    d1f(i,n2-2,k)=d1f(i,n2-2,k)+2.d0*(A(3)+A(4)+A(5))*f(i,n2,k)* g(n2-2)


                    d1f(i,n2-3,k)=  (A(5)*(-f(i,n2-2,k) - f(i,n2-7,k))     &
                    + A(4)*(-f(i,n2-1,k) - f(i,n2-6,k))             &
                    + A(3)*(f(i,n2,k) - f(i,n2-5,k))             &
                    + A(2)*(f(i,n2-1,k) - f(i,n2-4,k))             &
                    + A(1)*(f(i,n2-2,k) - f(i,n2-3,k)) ) * g(n2-3)

                    d1f(i,n2-3,k)=d1f(i,n2-3,k)+2.d0*(A(4)+A(5))*f(i,n2,k)* g(n2-3)


                    d1f(i,n2-4,k)=  (A(5)*(-f(i,n2-1,k) - f(i,n2-8,k))     &
                    + A(4)*(f(i,n2,k) - f(i,n2-7,k))             &
                    + A(3)*(f(i,n2-1,k) - f(i,n2-6,k))             &
                    + A(2)*(f(i,n2-2,k) - f(i,n2-5,k))             &
                    + A(1)*(f(i,n2-3,k) - f(i,n2-4,k)) ) * g(n2-4)

                    d1f(i,n2-4,k)=d1f(i,n2-4,k)+2.d0*(A(5))*f(i,n2,k)* g(n2-4)


                    d1f(i,n2-5,k)=  (A(5)*(f(i,n2,k) - f(i,n2-9,k))     &
                    + A(4)*(f(i,n2-1,k) - f(i,n2-8,k))             &
                    + A(3)*(f(i,n2-2,k) - f(i,n2-7,k))             &
                    + A(2)*(f(i,n2-3,k) - f(i,n2-6,k))             &
                    + A(1)*(f(i,n2-4,k) - f(i,n2-5,k)) ) * g(n2-5)

                enddo
            enddo

        end if


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e
                    d1f(i,1,k)= ((f(i,2,k)-f(i,1,k))/h)*g(1)
                    d1f(i,2,k)  =  (9.0d0/8.d0  * (f(i,3,k)-f(i,2,k))            / h     &
                    -   1.d0/24.d0  * (f(i,4,k)-f(i,1,k))                       / h)*g(2)
                    d1f(i,3,k)  = (1.189101951200031d0  *   (f(i,4,k)-f(i,3,k))  / h     &
                    -   0.073717642266682d0  *   (f(i,5,k)-f(i,2,k))            / h     &
                    +   0.006410195120003d0 *   (f(i,6,k) - f(i,1,k))           / h)*g(3)
                    d1f(i,4,k)  = (1.189101951200031d0  *   (f(i,5,k)-f(i,4,k))  / h     &
                    -   0.073717642266682d0  *   (f(i,6,k)-f(i,3,k))            / h     &
                    +   0.006410195120003d0 *   (f(i,7,k) - f(i,2,k))           / h)*g(4)

                enddo
            enddo

            do k=1, n3e
                do i=1, n1e
                    d1f(i,n2-1,k)=((f(i,n2,k)-f(i,n2-1,k))/h)*g(n2-1)
                    d1f(i,n2-2,k)  =   (9.0d0/8.d0 * (f(i,n2-1,k)-f(i,n2-2,k))               /h  &
                    -   1.d0/24.d0    * (f(i,n2,k)-f(i,n2-3,k))   /h)*g(n2-2)
                    d1f(i,n2-3,k)  =   (1.189101951200031d0  *   (f(i,n2-2,k)-f(i,n2-3,k))   /h  &
                    -   0.073717642266682d0  *   (f(i,n2-1,k)-f(i,n2-4,k))                  /h  &
                    +   0.006410195120003d0 *   (f(i,n2,k) - f(i,n2-5,k))                   /h)*g(n2-3)
                    d1f(i,n2-4,k)  =   (1.189101951200031d0  *   (f(i,n2-3,k)-f(i,n2-4,k))   /h  &
                    -   0.073717642266682d0  *   (f(i,n2-2,k)-f(i,n2-5,k))                  /h  &
                    +   0.006410195120003d0 *   (f(i,n2-1,k) - f(i,n2-6,k))                 /h)*g(n2-4)
                    d1f(i,n2-5,k)  =   (A(1)  *   (f(i,n2-4,k)-f(i,n2-5,k))         &
                    +   A(2)  *   (f(i,n2-3,k)-f(i,n2-6,k))     &
                    +   A(3) *   (f(i,n2-2,k) - f(i,n2-7,k))    &
                    +   A(4) *   (f(i,n2-1,k) - f(i,n2-8,k))    &
                    +   A(5)*   (f(i,n2,k) - f(i,n2-9,k))       )*g(n2-5)

                enddo
            enddo

        endif

        return

    end subroutine D1s_DRP5_MULT_3Dy


    subroutine D1s_DRP5_MULTACC_3Dy(f, d1f, n1,n1e,n2,n3,n3e, h, shifted, bc_type, g)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)       :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)      :: d1f
        real(kind=8), dimension(:), intent(in)              :: g

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        do k=1, n3e
            do i=1, n1e
                do j=5,n2-6
                    d1f(i, j, k)= d1f(i, j, k)+(A(5)*(f(i, j+5, k) - f(i, j-4, k))     &
                    + A(4)*(f(i, j+4, k) - f(i, j-3, k))     &
                    + A(3)*(f(i, j+3, k) - f(i, j-2, k))     &
                    + A(2)*(f(i, j+2, k) - f(i, j-1, k))     &
                    + A(1)*(f(i, j+1, k) - f(i, j, k)))*g(j)
                enddo
            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=1, 4
                        d1f(i, j, k)=   d1f(i, j, k)+(A(5)*(f(i, j+5, k) - f(i, mod(n2+j-6,n2-1)+1, k))  &
                        + A(4)*(f(i, j+4, k) - f(i, mod(n2+j-5,n2-1)+1, k))  &
                        + A(3)*(f(i, j+3, k) - f(i, mod(n2+j-4,n2-1)+1, k))  &
                        + A(2)*(f(i, j+2, k) - f(i, mod(n2+j-3,n2-1)+1, k))  &
                        + A(1)*(f(i, j+1, k) - f(i, j, k)))*g(j)
                    enddo

                    do j=n2-5,n2-1
                        d1f(i, j, k)= d1f(i, j, k)+(A(5)*(f(i, mod(j+4,n2-1)+1, k) - f(i, j-4, k))             &
                        + A(4)*(f(i, mod(j+3,n2-1)+1, k) - f(i, j-3, k))             &
                        + A(3)*(f(i, mod(j+2,n2-1)+1, k) - f(i, j-2, k))             &
                        + A(2)*(f(i, mod(j+1,n2-1)+1, k) - f(i, j-1, k))             &
                        + A(1)*(f(i, mod(j,n2-1)+1, k) - f(i, j, k)))*g(j)
                    enddo

                enddo
            enddo

        endif

        if (bc_type.eq.antisymetric) then

            do k=1, n3e
                do i=1, n1e

                    d1f(i,1,k)=  d1f(i,1,k) + (A(5)*(f(i,6,k) + f(i,5,k))   &
                    + A(4)*(f(i,5,k) + f(i,4,k))           &
                    + A(3)*(f(i,4,k) + f(i,3,k))           &
                    + A(2)*(f(i,3,k) + f(i,2,k))           &
                    + A(1)*(f(i,2,k) - f(i,1,k)) ) * g(1)

                    d1f(i,1,k)=d1f(i,1,k)-2.d0*(A(2)+A(3)+A(4)+A(5))*f(i,1,k)* g(1)

                    d1f(i,2,k)=  d1f(i,2,k) + (A(5)*(f(i,7,k) + f(i,4,k))   &
                    + A(4)*(f(i,6,k) + f(i,3,k))           &
                    + A(3)*(f(i,5,k) + f(i,2,k))           &
                    + A(2)*(f(i,4,k) - f(i,1,k))           &
                    + A(1)*(f(i,3,k) - f(i,2,k)) ) * g(2)

                    d1f(i,2,k)=d1f(i,2,k)-2.d0*(A(3)+A(4)+A(5))*f(i,1,k)* g(2)

                    d1f(i,3,k)=  d1f(i,3,k) + (A(5)*(f(i,8,k) + f(i,3,k))   &
                    + A(4)*(f(i,7,k) + f(i,2,k))           &
                    + A(3)*(f(i,6,k) - f(i,1,k))           &
                    + A(2)*(f(i,5,k) - f(i,2,k))           &
                    + A(1)*(f(i,4,k) - f(i,3,k)) ) * g(3)

                    d1f(i,3,k)=d1f(i,3,k)-2.d0*(A(4)+A(5))*f(i,1,k)* g(3)

                    d1f(i,4,k)=  d1f(i,4,k) + (A(5)*(f(i,9,k) + f(i,2,k))   &
                    + A(4)*(f(i,8,k) - f(i,1,k))           &
                    + A(3)*(f(i,7,k) - f(i,2,k))           &
                    + A(2)*(f(i,6,k) - f(i,3,k))           &
                    + A(1)*(f(i,5,k) - f(i,4,k)) ) * g(4)

                    d1f(i,4,k)=d1f(i,4,k)-2.d0*(A(5))*f(i,1,k)* g(4)


                    d1f(i,n2-1,k)=  d1f(i,n2-1,k) + (A(5)*(-f(i,n2-4,k) - f(i,n2-5,k))     &
                    + A(4)*(-f(i,n2-3,k) - f(i,n2-4,k))             &
                    + A(3)*(-f(i,n2-2,k) - f(i,n2-3,k))             &
                    + A(2)*(-f(i,n2-1,k) - f(i,n2-2,k))             &
                    + A(1)*(f(i,n2,k) - f(i,n2-1,k)) ) * g(n2-1)

                    d1f(i,n2-1,k)=d1f(i,n2-1,k)+2.d0*(A(2)+A(3)+A(4)+A(5))*f(i,n2,k)* g(n2-1)


                    d1f(i,n2-2,k)=  d1f(i,n2-2,k) + (A(5)*(-f(i,n2-3,k) - f(i,n2-6,k))     &
                    + A(4)*(-f(i,n2-2,k) - f(i,n2-5,k))             &
                    + A(3)*(-f(i,n2-1,k) - f(i,n2-4,k))             &
                    + A(2)*(f(i,n2,k) - f(i,n2-3,k))             &
                    + A(1)*(f(i,n2-1,k) - f(i,n2-2,k)) ) * g(n2-2)

                    d1f(i,n2-2,k)=d1f(i,n2-2,k)+2.d0*(A(3)+A(4)+A(5))*f(i,n2,k)* g(n2-2)


                    d1f(i,n2-3,k)=  d1f(i,n2-3,k) + (A(5)*(-f(i,n2-2,k) - f(i,n2-7,k))     &
                    + A(4)*(-f(i,n2-1,k) - f(i,n2-6,k))             &
                    + A(3)*(f(i,n2,k) - f(i,n2-5,k))             &
                    + A(2)*(f(i,n2-1,k) - f(i,n2-4,k))             &
                    + A(1)*(f(i,n2-2,k) - f(i,n2-3,k)) ) * g(n2-3)

                    d1f(i,n2-3,k)=d1f(i,n2-3,k)+2.d0*(A(4)+A(5))*f(i,n2,k)* g(n2-3)


                    d1f(i,n2-4,k)=  d1f(i,n2-4,k) + (A(5)*(-f(i,n2-1,k) - f(i,n2-8,k))     &
                    + A(4)*(f(i,n2,k) - f(i,n2-7,k))             &
                    + A(3)*(f(i,n2-1,k) - f(i,n2-6,k))             &
                    + A(2)*(f(i,n2-2,k) - f(i,n2-5,k))             &
                    + A(1)*(f(i,n2-3,k) - f(i,n2-4,k)) ) * g(n2-4)

                    d1f(i,n2-4,k)=d1f(i,n2-4,k)+2.d0*(A(5))*f(i,n2,k)* g(n2-4)


                    d1f(i,n2-5,k)=  d1f(i,n2-5,k) + (A(5)*(f(i,n2,k) - f(i,n2-9,k))     &
                    + A(4)*(f(i,n2-1,k) - f(i,n2-8,k))             &
                    + A(3)*(f(i,n2-2,k) - f(i,n2-7,k))             &
                    + A(2)*(f(i,n2-3,k) - f(i,n2-6,k))             &
                    + A(1)*(f(i,n2-4,k) - f(i,n2-5,k)) ) * g(n2-5)

                enddo
            enddo

        end if


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e
                    d1f(i,1,k)= d1f(i,1,k)+((f(i,2,k)-f(i,1,k))/h)*g(1)
                    d1f(i,2,k)  =  d1f(i,2,k)+(9.0d0/8.d0  * (f(i,3,k)-f(i,2,k))            / h     &
                    -   1.d0/24.d0  * (f(i,4,k)-f(i,1,k))                       / h)*g(2)
                    d1f(i,3,k)  = d1f(i,3,k)+(1.189101951200031d0  *   (f(i,4,k)-f(i,3,k))  / h     &
                    -   0.073717642266682d0  *   (f(i,5,k)-f(i,2,k))            / h     &
                    +   0.006410195120003d0 *   (f(i,6,k) - f(i,1,k))           / h)*g(3)
                    d1f(i,4,k)  = d1f(i,4,k)+(1.189101951200031d0  *   (f(i,5,k)-f(i,4,k))  / h     &
                    -   0.073717642266682d0  *   (f(i,6,k)-f(i,3,k))            / h     &
                    +   0.006410195120003d0 *   (f(i,7,k) - f(i,2,k))           / h)*g(4)

                enddo
            enddo

            do k=1, n3e
                do i=1, n1e
                    d1f(i,n2-1,k)=d1f(i,n2-1,k)+((f(i,n2,k)-f(i,n2-1,k))/h)*g(n2-1)
                    d1f(i,n2-2,k)  =   d1f(i,n2-2,k)+(9.0d0/8.d0 * (f(i,n2-1,k)-f(i,n2-2,k))               /h  &
                    -   1.d0/24.d0    * (f(i,n2,k)-f(i,n2-3,k))   /h)*g(n2-2)
                    d1f(i,n2-3,k)  =   d1f(i,n2-3,k)+(1.189101951200031d0  *   (f(i,n2-2,k)-f(i,n2-3,k))   /h  &
                    -   0.073717642266682d0  *   (f(i,n2-1,k)-f(i,n2-4,k))                  /h  &
                    +   0.006410195120003d0 *   (f(i,n2,k) - f(i,n2-5,k))                   /h)*g(n2-3)
                    d1f(i,n2-4,k)  =   d1f(i,n2-4,k)+(1.189101951200031d0  *   (f(i,n2-3,k)-f(i,n2-4,k))   /h  &
                    -   0.073717642266682d0  *   (f(i,n2-2,k)-f(i,n2-5,k))                  /h  &
                    +   0.006410195120003d0 *   (f(i,n2-1,k) - f(i,n2-6,k))                 /h)*g(n2-4)
                    d1f(i,n2-5,k)  =   d1f(i,n2-5,k)+(A(1)  *   (f(i,n2-4,k)-f(i,n2-5,k))   &
                    +   A(2)  *   (f(i,n2-3,k)-f(i,n2-6,k))         &
                    +   A(3) *   (f(i,n2-2,k) - f(i,n2-7,k))        &
                    +   A(4) *   (f(i,n2-1,k) - f(i,n2-8,k))        &
                    +   A(5)*   (f(i,n2,k) - f(i,n2-9,k))           )*g(n2-5)

                enddo
            enddo

        endif

        return

    end subroutine D1s_DRP5_MULTACC_3Dy

    subroutine D1s_DRP5_3Dz(f, d1f, n1,n1e,n2,n2e,n3, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: d1f

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        do j=1, n2e
            do i=1, n1e
                do k=5,n3-6
                    d1f(i, j, k)= A(5)*(f(i, j, k+5) - f(i, j, k-4))     &
                    + A(4)*(f(i, j, k+4) - f(i, j, k-3))     &
                    + A(3)*(f(i, j, k+3) - f(i, j, k-2))     &
                    + A(2)*(f(i, j, k+2) - f(i, j, k-1))     &
                    + A(1)*(f(i, j, k+1) - f(i, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e

                    do k=1, 4
                        d1f(i, j, k)=   A(5)*(f(i, j, k+5) - f(i, j, mod(n3+k-6,n3-1)+1))  &
                        + A(4)*(f(i, j, k+4) - f(i, j, mod(n3+k-5,n3-1)+1))  &
                        + A(3)*(f(i, j, k+3) - f(i, j, mod(n3+k-4,n3-1)+1))  &
                        + A(2)*(f(i, j, k+2) - f(i, j, mod(n3+k-3,n3-1)+1))  &
                        + A(1)*(f(i, j, k+1) - f(i, j, k))
                    enddo

                    do k=n3-5,n3-1
                        d1f(i, j, k)= A(5)*(f(i, j, mod(k+4,n3-1)+1) - f(i, j, k-4))             &
                        + A(4)*(f(i, j, mod(k+3,n3-1)+1) - f(i, j, k-3))             &
                        + A(3)*(f(i, j, mod(k+2,n3-1)+1) - f(i, j, k-2))             &
                        + A(2)*(f(i, j, mod(k+1,n3-1)+1) - f(i, j, k-1))             &
                        + A(1)*(f(i, j, mod(k,n3-1)+1) - f(i, j, k))
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,1)= (f(i,j,2)-f(i,j,1))/h
                    d1f(i,j,2)  =  9.0d0/8.d0  * (f(i,j,3)-f(i,j,2))     / h     &
                    -   1.d0/24.d0  * (f(i,j,4)-f(i,j,1))           / h
                    d1f(i,j,3)  = 1.189101951200031d0  *   (f(i,j,4)-f(i,j,3))  / h         &
                    -   0.073717642266682d0  *   (f(i,j,5)-f(i,j,2))        / h         &
                    +   0.006410195120003d0 *   (f(i,j,6) - f(i,j,1))       / h
                    d1f(i,j,4)  = 1.189101951200031d0  *   (f(i,j,5)-f(i,j,4))  / h         &
                    -   0.073717642266682d0  *   (f(i,j,6)-f(i,j,3))        / h         &
                    +   0.006410195120003d0 *   (f(i,j,7) - f(i,j,2))       / h
                    d1f(i,j,5)  = A(1)  *   (f(i,j,6)-f(i,j,5))     &
                    +   A(2)  *   (f(i,j,7)-f(i,j,4))               &
                    +   A(3) *   (f(i,j,8) - f(i,j,3))              &
                    +   A(4) *   (f(i,j,9) - f(i,j,2))              &
                    +   A(5)*   (f(i,j,10) - f(i,j,1))
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,n3-1)=(f(i,j,n3)-f(i,j,n3-1))/h
                    d1f(i,j,n3-2)  =   9.0d0/8.d0 * (f(i,j,n3-1)-f(i,j,n3-2))  /h    &
                    -   1.d0/24.d0    * (f(i,j,n3)-f(i,j,n3-3))   /h
                    d1f(i,j,n3-3)  =   1.189101951200031d0  *   (f(i,j,n3-2)-f(i,j,n3-3))   /h   &
                    -   0.073717642266682d0  *   (f(i,j,n3-1)-f(i,j,n3-4))                   /h  &
                    +   0.006410195120003d0 *   (f(i,j,n3) - f(i,j,n3-5))                   /h
                    d1f(i,j,n3-4)  =   1.189101951200031d0  *   (f(i,j,n3-3)-f(i,j,n3-4))    /h  &
                    -   0.073717642266682d0  *   (f(i,j,n3-2)-f(i,j,n3-5))                   /h  &
                    +   0.006410195120003d0 *   (f(i,j,n3-1) - f(i,j,n3-6))                 /h
                    d1f(i,j,n3-5)  =   A(1)  *   (f(i,j,n3-4)-f(i,j,n3-5))      &
                    +   A(2)  *   (f(i,j,n3-3)-f(i,j,n3-6))                     &
                    +   A(3) *   (f(i,j,n3-2) - f(i,j,n3-7))                    &
                    +   A(4) *   (f(i,j,n3-1) - f(i,j,n3-8))                    &
                    +   A(5)*   (f(i,j,n3) - f(i,j,n3-9))
                enddo
            enddo

        endif

        return

    end subroutine D1s_DRP5_3Dz

    subroutine D1s_DRP5_ACC_3Dz(f, d1f, n1,n1e,n2,n2e,n3, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: d1f

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        do j=1, n2e
            do i=1, n1e
                do k=5,n3-6
                    d1f(i, j, k)= d1f(i, j, k)+A(5)*(f(i, j, k+5) - f(i, j, k-4))     &
                    + A(4)*(f(i, j, k+4) - f(i, j, k-3))     &
                    + A(3)*(f(i, j, k+3) - f(i, j, k-2))     &
                    + A(2)*(f(i, j, k+2) - f(i, j, k-1))     &
                    + A(1)*(f(i, j, k+1) - f(i, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e

                    do k=1, 4
                        d1f(i, j, k)=   d1f(i, j, k)+A(5)*(f(i, j, k+5) - f(i, j, mod(n3+k-6,n3-1)+1))  &
                        + A(4)*(f(i, j, k+4) - f(i, j, mod(n3+k-5,n3-1)+1))  &
                        + A(3)*(f(i, j, k+3) - f(i, j, mod(n3+k-4,n3-1)+1))  &
                        + A(2)*(f(i, j, k+2) - f(i, j, mod(n3+k-3,n3-1)+1))  &
                        + A(1)*(f(i, j, k+1) - f(i, j, k))
                    enddo

                    do k=n3-5,n3-1
                        d1f(i, j, k)= d1f(i, j, k)+A(5)*(f(i, j, mod(k+4,n3-1)+1) - f(i, j, k-4))             &
                        + A(4)*(f(i, j, mod(k+3,n3-1)+1) - f(i, j, k-3))             &
                        + A(3)*(f(i, j, mod(k+2,n3-1)+1) - f(i, j, k-2))             &
                        + A(2)*(f(i, j, mod(k+1,n3-1)+1) - f(i, j, k-1))             &
                        + A(1)*(f(i, j, mod(k,n3-1)+1) - f(i, j, k))
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,1)= d1f(i,j,1)+(f(i,j,2)-f(i,j,1))/h
                    d1f(i,j,2)  =  d1f(i,j,2)+9.0d0/8.d0  * (f(i,j,3)-f(i,j,2))     / h     &
                    -   1.d0/24.d0  * (f(i,j,4)-f(i,j,1))           / h
                    d1f(i,j,3)  = d1f(i,j,3)+1.189101951200031d0  *   (f(i,j,4)-f(i,j,3))  / h         &
                    -   0.073717642266682d0  *   (f(i,j,5)-f(i,j,2))        / h         &
                    +   0.006410195120003d0 *   (f(i,j,6) - f(i,j,1))       / h
                    d1f(i,j,4)  = d1f(i,j,4)+1.189101951200031d0  *   (f(i,j,5)-f(i,j,4))  / h         &
                    -   0.073717642266682d0  *   (f(i,j,6)-f(i,j,3))        / h         &
                    +   0.006410195120003d0 *   (f(i,j,7) - f(i,j,2))       / h
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,n3-1)=d1f(i,j,n3-1)+(f(i,j,n3)-f(i,j,n3-1))/h
                    d1f(i,j,n3-2)  =   d1f(i,j,n3-2)+9.0d0/8.d0 * (f(i,j,n3-1)-f(i,j,n3-2))  /h    &
                    -   1.d0/24.d0    * (f(i,j,n3)-f(i,j,n3-3))   /h
                    d1f(i,j,n3-3)  =   d1f(i,j,n3-3)+1.189101951200031d0  *   (f(i,j,n3-2)-f(i,j,n3-3))   /h   &
                    -   0.073717642266682d0  *   (f(i,j,n3-1)-f(i,j,n3-4))                   /h  &
                    +   0.006410195120003d0 *   (f(i,j,n3) - f(i,j,n3-5))                   /h
                    d1f(i,j,n3-4)  =   d1f(i,j,n3-4)+1.189101951200031d0  *   (f(i,j,n3-3)-f(i,j,n3-4))    /h  &
                    -   0.073717642266682d0  *   (f(i,j,n3-2)-f(i,j,n3-5))                   /h  &
                    +   0.006410195120003d0 *   (f(i,j,n3-1) - f(i,j,n3-6))                 /h
                    d1f(i,j,n3-5)  =   d1f(i,j,n3-5)+A(1)  *   (f(i,j,n3-4)-f(i,j,n3-5))    &
                    +   A(2)  *   (f(i,j,n3-3)-f(i,j,n3-6))                     &
                    +   A(3) *   (f(i,j,n3-2) - f(i,j,n3-7))                    &
                    +   A(4) *   (f(i,j,n3-1) - f(i,j,n3-8))                    &
                    +   A(5)*   (f(i,j,n3) - f(i,j,n3-9))
                enddo
            enddo

        endif

        return

    end subroutine D1s_DRP5_ACC_3Dz

end module DRP_D1s

module DRP_D1s_sh
    use boundaries_types
    implicit none

contains

    subroutine D1ssh_DRP5_3Dx(f, d1f, n1,n2,n2e,n3,n3e, h, shifted, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: d1f

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        do k=1,n3e
            do j=1,n2e
                do i=6,n1-5
                    d1f(i, j, k)= A(5)*(f(i+4, j, k) - f(i-5, j, k))     &
                    + A(4)*(f(i+3, j, k) - f(i-4, j, k))     &
                    + A(3)*(f(i+2, j, k) - f(i-3, j, k))     &
                    + A(2)*(f(i+1, j, k) - f(i-2, j, k))     &
                    + A(1)*(f(i, j, k) - f(i-1, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do k=1,n3e
                do j=1,n2e

                    do i=1, 5
                        d1f(i, j, k)=   A(5)*(f(i+4, j, k) - f(mod(n1+i-7,n1-1)+1, j, k))  &
                        + A(4)*(f(i+3, j, k) - f(mod(n1+i-6,n1-1)+1, j, k))  &
                        + A(3)*(f(i+2, j, k) - f(mod(n1+i-5,n1-1)+1, j, k))  &
                        + A(2)*(f(i+1, j, k) - f(mod(n1+i-4,n1-1)+1, j, k))  &
                        + A(1)*(f(i, j, k) - f(mod(n1+i-3,n1-1)+1, j, k))
                    enddo

                    do i=n1-4,n1-1
                        d1f(i, j, k)= A(5)*(f(mod(i+3,n1-1)+1, j, k) - f(i-5, j, k))             &
                        + A(4)*(f(mod(i+2,n1-1)+1, j, k) - f(i-4, j, k))             &
                        + A(3)*(f(mod(i+1,n1-1)+1, j, k) - f(i-3, j, k))             &
                        + A(2)*(f(mod(i,n1-1)+1, j, k) - f(i-2, j, k))             &
                        + A(1)*(f(i, j, k) - f(i-1, j, k))
                    enddo
                enddo

            enddo

        endif

        if (bc_type==symetric) then
            do k=1, n3e
                do j=1, n2e

                    do i=1, 5
                        d1f(i,j,k)=   A(5)*(f(i+4,j,k) - f(max(i-5,abs(i-6)),j,k))   &
                        + A(4)*(f(i+3,j,k) - f(max(i-4,abs(i-5)),j,k))           &
                        + A(3)*(f(i+2,j,k) - f(max(i-3,abs(i-4)),j,k))           &
                        + A(2)*(f(i+1,j,k) - f(max(i-2,abs(i-3)),j,k))           &
                        + A(1)*(f(i,j,k) - f(max(i-1,abs(i-2)),j,k))
                    enddo

                    do i=n1-4,n1-1
                        d1f(i,j,k)=   A(5)*(f(min(i+4, 2*n1-(i+5)),j,k) - f(i-5,j,k))     &
                        + A(4)*(f(min(i+3, 2*n1-(i+4)),j,k) - f(i-4,j,k))             &
                        + A(3)*(f(min(i+2, 2*n1-(i+3)),j,k) - f(i-3,j,k))             &
                        + A(2)*(f(min(i+1, 2*n1-(i+2)),j,k) - f(i-2,j,k))             &
                        + A(1)*(f(i,j,k) - f(i-1,j,k))
                    enddo

                enddo
            enddo

        end if


        if ((bc_type.eq.Dirichlet)) then

            do k=1,n3e
                do j=1,n2e
                    d1f(2,j,k)= (f(2,j,k)-f(1,j,k))                                 /h
                    d1f(3,j,k)  = 9.0d0/8.d0  * (f(3,j,k)-f(2,j,k))                 / h         &
                    -   1.d0/24.d0  * (f(4,j,k)-f(1,j,k))                           / h
                    d1f(4,j,k)  = 1.189101951200031d0  *   (f(4,j,k)-f(3,j,k))      / h         &
                    -   0.073717642266682d0  *   (f(5,j,k)-f(2,j,k))                / h         &
                    +   0.006410195120003d0 *   (f(6,j,k) - f(1,j,k))               / h
                    d1f(5,j,k)  = 1.189101951200031d0  *   (f(5,j,k)-f(4,j,k))      / h         &
                    -   0.073717642266682d0  *   (f(6,j,k)-f(3,j,k))                / h         &
                    +   0.006410195120003d0 *   (f(7,j,k) - f(2,j,k))               / h

                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    d1f(n1-1,j,k)=(f(n1-1,j,k)-f(n1-2,j,k))                                 /h
                    d1f(n1-2,j,k)  =   9.0d0/8.d0 * (f(n1-2,j,k)-f(n1-3,j,k))               /h      &
                    -   1.d0/24.d0    * (f(n1-1,j,k)-f(n1-4,j,k))                           /h
                    d1f(n1-3,j,k)  =   1.189101951200031d0  *   (f(n1-3,j,k)-f(n1-4,j,k))   /h      &
                    -   0.073717642266682d0  *   (f(n1-2,j,k)-f(n1-5,j,k))                  /h      &
                    +   0.006410195120003d0 *   (f(n1-1,j,k) - f(n1-6,j,k))                 /h
                    d1f(n1-4,j,k)  =   1.189101951200031d0  *   (f(n1-4,j,k)-f(n1-5,j,k))   /h      &
                    -   0.073717642266682d0  *   (f(n1-3,j,k)-f(n1-6,j,k))                  /h      &
                    +   0.006410195120003d0 *   (f(n1-2,j,k) - f(n1-7,j,k))                 /h
                enddo
            enddo


        endif

        return

    end subroutine D1ssh_DRP5_3Dx

    subroutine D1ssh_DRP5_ACC_3Dx(f, d1f, n1,n2,n2e,n3,n3e, h, shifted, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: d1f

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        do k=1,n3e
            do j=1,n2e
                do i=6,n1-5
                    d1f(i, j, k)= d1f(i, j, k)+A(5)*(f(i+4, j, k) - f(i-5, j, k))     &
                    + A(4)*(f(i+3, j, k) - f(i-4, j, k))     &
                    + A(3)*(f(i+2, j, k) - f(i-3, j, k))     &
                    + A(2)*(f(i+1, j, k) - f(i-2, j, k))     &
                    + A(1)*(f(i, j, k) - f(i-1, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do k=1,n3e
                do j=1,n2e

                    do i=1, 5
                        d1f(i, j, k)=   d1f(i, j, k)+A(5)*(f(i+4, j, k) - f(mod(n1+i-7,n1-1)+1, j, k))  &
                        + A(4)*(f(i+3, j, k) - f(mod(n1+i-6,n1-1)+1, j, k))  &
                        + A(3)*(f(i+2, j, k) - f(mod(n1+i-5,n1-1)+1, j, k))  &
                        + A(2)*(f(i+1, j, k) - f(mod(n1+i-4,n1-1)+1, j, k))  &
                        + A(1)*(f(i, j, k) - f(mod(n1+i-3,n1-1)+1, j, k))
                    enddo

                    do i=n1-4,n1-1
                        d1f(i, j, k)= d1f(i, j, k)+A(5)*(f(mod(i+3,n1-1)+1, j, k) - f(i-5, j, k))             &
                        + A(4)*(f(mod(i+2,n1-1)+1, j, k) - f(i-4, j, k))             &
                        + A(3)*(f(mod(i+1,n1-1)+1, j, k) - f(i-3, j, k))             &
                        + A(2)*(f(mod(i,n1-1)+1, j, k) - f(i-2, j, k))             &
                        + A(1)*(f(i, j, k) - f(i-1, j, k))
                    enddo
                enddo

            enddo

        endif

        if (bc_type==symetric) then
            do k=1, n3e
                do j=1, n2e

                    do i=1, 5
                        d1f(i,j,k)=   d1f(i,j,k)+A(5)*(f(i+4,j,k) - f(max(i-5,abs(i-6)),j,k))   &
                        + A(4)*(f(i+3,j,k) - f(max(i-4,abs(i-5)),j,k))           &
                        + A(3)*(f(i+2,j,k) - f(max(i-3,abs(i-4)),j,k))           &
                        + A(2)*(f(i+1,j,k) - f(max(i-2,abs(i-3)),j,k))           &
                        + A(1)*(f(i,j,k) - f(max(i-1,abs(i-2)),j,k))
                    enddo

                    do i=n1-4,n1-1
                        d1f(i,j,k)=   d1f(i,j,k)+A(5)*(f(min(i+4, 2*n1-(i+5)),j,k) - f(i-5,j,k))     &
                        + A(4)*(f(min(i+3, 2*n1-(i+4)),j,k) - f(i-4,j,k))             &
                        + A(3)*(f(min(i+2, 2*n1-(i+3)),j,k) - f(i-3,j,k))             &
                        + A(2)*(f(min(i+1, 2*n1-(i+2)),j,k) - f(i-2,j,k))             &
                        + A(1)*(f(i,j,k) - f(i-1,j,k))
                    enddo

                enddo
            enddo

        end if


        if ((bc_type.eq.Dirichlet)) then

            do k=1,n3e
                do j=1,n2e
                    d1f(2,j,k)= d1f(2,j,k)+(f(2,j,k)-f(1,j,k))                                 /h
                    d1f(3,j,k)  = d1f(3,j,k)+9.0d0/8.d0  * (f(3,j,k)-f(2,j,k))                 / h         &
                    -   1.d0/24.d0  * (f(4,j,k)-f(1,j,k))                           / h
                    d1f(4,j,k)  = d1f(4,j,k)+1.189101951200031d0  *   (f(4,j,k)-f(3,j,k))      / h         &
                    -   0.073717642266682d0  *   (f(5,j,k)-f(2,j,k))                / h         &
                    +   0.006410195120003d0 *   (f(6,j,k) - f(1,j,k))               / h
                    d1f(5,j,k)  = d1f(5,j,k)+1.189101951200031d0  *   (f(5,j,k)-f(4,j,k))      / h         &
                    -   0.073717642266682d0  *   (f(6,j,k)-f(3,j,k))                / h         &
                    +   0.006410195120003d0 *   (f(7,j,k) - f(2,j,k))               / h

                enddo
            enddo

            do k=1,n3e
                do j=1,n2e
                    d1f(n1-1,j,k)=d1f(n1-1,j,k)+(f(n1-1,j,k)-f(n1-2,j,k))                                 /h
                    d1f(n1-2,j,k)  =  d1f(n1-2,j,k)+ 9.0d0/8.d0 * (f(n1-2,j,k)-f(n1-3,j,k))               /h      &
                    -   1.d0/24.d0    * (f(n1-1,j,k)-f(n1-4,j,k))                           /h
                    d1f(n1-3,j,k)  =   d1f(n1-3,j,k)+1.189101951200031d0  *   (f(n1-3,j,k)-f(n1-4,j,k))   /h      &
                    -   0.073717642266682d0  *   (f(n1-2,j,k)-f(n1-5,j,k))                  /h      &
                    +   0.006410195120003d0 *   (f(n1-1,j,k) - f(n1-6,j,k))                 /h
                    d1f(n1-4,j,k)  =  d1f(n1-4,j,k)+ 1.189101951200031d0  *   (f(n1-4,j,k)-f(n1-5,j,k))   /h      &
                    -   0.073717642266682d0  *   (f(n1-3,j,k)-f(n1-6,j,k))                  /h      &
                    +   0.006410195120003d0 *   (f(n1-2,j,k) - f(n1-7,j,k))                 /h
                enddo
            enddo


        endif

        return

    end subroutine D1ssh_DRP5_ACC_3Dx



    subroutine D1ssh_DRP5_3Dy(f, d1f, n1,n1e,n2,n3,n3e, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: d1f

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        do k=1, n3e
            do i=1, n1e
                do j=5,n2-6
                    d1f(i, j+1, k)= A(5)*(f(i, j+5, k) - f(i, j-4, k))     &
                    + A(4)*(f(i, j+4, k) - f(i, j-3, k))     &
                    + A(3)*(f(i, j+3, k) - f(i, j-2, k))     &
                    + A(2)*(f(i, j+2, k) - f(i, j-1, k))     &
                    + A(1)*(f(i, j+1, k) - f(i, j, k))
                enddo
            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=1, 5
                        d1f(i, j, k)=   A(5)*(f(i, j+4, k) - f(i, mod(n2+j-7,n2-1)+1, k))  &
                        + A(4)*(f(i, j+3, k) - f(i, mod(n2+j-6,n2-1)+1, k))  &
                        + A(3)*(f(i, j+2, k) - f(i, mod(n2+j-5,n2-1)+1, k))  &
                        + A(2)*(f(i, j+1, k) - f(i, mod(n2+j-4,n2-1)+1, k))  &
                        + A(1)*(f(i, j, k) - f(i, mod(n2+j-3,n2-1)+1, k))
                    enddo

                    do j=n2-4,n2-1
                        d1f(i, j, k)= A(5)*(f(i, mod(j+3,n2-1)+1, k) - f(i, j-5, k))             &
                        + A(4)*(f(i, mod(j+2,n2-1)+1, k) - f(i, j-4, k))             &
                        + A(3)*(f(i, mod(j+1,n2-1)+1, k) - f(i, j-3, k))             &
                        + A(2)*(f(i, mod(j,n2-1)+1, k) - f(i, j-2, k))             &
                        + A(1)*(f(i, j, k) - f(i, j-1, k))
                    enddo

                enddo
            enddo

        endif

        if (bc_type==symetric) then
            do k=1, n3e
                do i=1, n1e

                    do j=1, 5
                        d1f(i,j,k)=   A(5)*(f(i,j+4,k) - f(i,max(j-5,abs(j-6)),k))   &
                        + A(4)*(f(i,j+3,k) - f(i,max(j-4,abs(j-5)),k))           &
                        + A(3)*(f(i,j+2,k) - f(i,max(j-3,abs(j-4)),k))           &
                        + A(2)*(f(i,j+1,k) - f(i,max(j-2,abs(j-3)),k))           &
                        + A(1)*(f(i,j,k) - f(i,max(j-1,abs(j-2)),k))
                    enddo

                    do j=n2-4,n2-1
                        d1f(i,j,k)=   A(5)*(f(i,min(j+4, 2*n2-(j+5)),k) - f(i,j-5,k))     &
                        + A(4)*(f(i,min(j+3, 2*n2-(j+4)),k) - f(i,j-4,k))             &
                        + A(3)*(f(i,min(j+2, 2*n2-(j+3)),k) - f(i,j-3,k))             &
                        + A(2)*(f(i,min(j+1, 2*n2-(j+2)),k) - f(i,j-2,k))             &
                        + A(1)*(f(i,j,k) - f(i,j-1,k))
                    enddo

                enddo
            enddo

        end if


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e
                    d1f(i,2,k)= (f(i,2,k)-f(i,1,k))                                 /h
                    d1f(i,3,k)  = 9.0d0/8.d0  * (f(i,3,k)-f(i,2,k))                 / h         &
                    -   1.d0/24.d0  * (f(i,4,k)-f(i,1,k))                           / h
                    d1f(i,4,k)  = 1.189101951200031d0  *   (f(i,4,k)-f(i,3,k))      / h         &
                    -   0.073717642266682d0  *   (f(i,5,k)-f(i,2,k))                / h         &
                    +   0.006410195120003d0 *   (f(i,6,k) - f(i,1,k))               / h
                    d1f(i,5,k)  = 1.189101951200031d0  *   (f(i,5,k)-f(i,4,k))      / h         &
                    -   0.073717642266682d0  *   (f(i,6,k)-f(i,3,k))                / h         &
                    +   0.006410195120003d0 *   (f(i,7,k) - f(i,2,k))               / h

                enddo
            enddo

            do k=1, n3e
                do i=1, n1e
                    d1f(i,n2-1,k)=(f(i,n2-1,k)-f(i,n2-2,k))                                 /h
                    d1f(i,n2-2,k)  =   9.0d0/8.d0 * (f(i,n2-2,k)-f(i,n2-3,k))               /h      &
                    -   1.d0/24.d0    * (f(i,n2-1,k)-f(i,n2-4,k))                           /h
                    d1f(i,n2-3,k)  =   1.189101951200031d0  *   (f(i,n2-3,k)-f(i,n2-4,k))   /h      &
                    -   0.073717642266682d0  *   (f(i,n2-2,k)-f(i,n2-5,k))                  /h      &
                    +   0.006410195120003d0 *   (f(i,n2-1,k) - f(i,n2-6,k))                 /h
                    d1f(i,n2-4,k)  =   1.189101951200031d0  *   (f(i,n2-4,k)-f(i,n2-5,k))   /h      &
                    -   0.073717642266682d0  *   (f(i,n2-3,k)-f(i,n2-6,k))                  /h      &
                    +   0.006410195120003d0 *   (f(i,n2-2,k) - f(i,n2-7,k))                 /h

                enddo
            enddo

        endif

        return

    end subroutine D1ssh_DRP5_3Dy

    subroutine D1ssh_DRP5_MULT_3Dy(f, d1f, n1,n1e,n2,n3,n3e, h, shifted, bc_type, g)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)       :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)      :: d1f
        real(kind=8), dimension(:), intent(in)              :: g

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        do k=1, n3e
            do i=1, n1e
                do j=5,n2-6
                    d1f(i, j+1, k)= (A(5)*(f(i, j+5, k) - f(i, j-4, k))     &
                    + A(4)*(f(i, j+4, k) - f(i, j-3, k))     &
                    + A(3)*(f(i, j+3, k) - f(i, j-2, k))     &
                    + A(2)*(f(i, j+2, k) - f(i, j-1, k))     &
                    + A(1)*(f(i, j+1, k) - f(i, j, k)))*g(j+1)
                enddo
            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=1, 5
                        d1f(i, j, k)=   (A(5)*(f(i, j+4, k) - f(i, mod(n2+j-7,n2-1)+1, k))  &
                        + A(4)*(f(i, j+3, k) - f(i, mod(n2+j-6,n2-1)+1, k))  &
                        + A(3)*(f(i, j+2, k) - f(i, mod(n2+j-5,n2-1)+1, k))  &
                        + A(2)*(f(i, j+1, k) - f(i, mod(n2+j-4,n2-1)+1, k))  &
                        + A(1)*(f(i, j, k) - f(i, mod(n2+j-3,n2-1)+1, k)))*g(j)
                    enddo

                    do j=n2-4,n2-1
                        d1f(i, j, k)= (A(5)*(f(i, mod(j+3,n2-1)+1, k) - f(i, j-5, k))             &
                        + A(4)*(f(i, mod(j+2,n2-1)+1, k) - f(i, j-4, k))             &
                        + A(3)*(f(i, mod(j+1,n2-1)+1, k) - f(i, j-3, k))             &
                        + A(2)*(f(i, mod(j,n2-1)+1, k) - f(i, j-2, k))             &
                        + A(1)*(f(i, j, k) - f(i, j-1, k)))*g(j)
                    enddo

                enddo
            enddo

        endif

        if (bc_type==symetric) then
            do k=1, n3e
                do i=1, n1e

                    do j=1, 5
                        d1f(i,j,k)=   (A(5)*(f(i,j+4,k) - f(i,max(j-5,abs(j-6)),k))   &
                        + A(4)*(f(i,j+3,k) - f(i,max(j-4,abs(j-5)),k))           &
                        + A(3)*(f(i,j+2,k) - f(i,max(j-3,abs(j-4)),k))           &
                        + A(2)*(f(i,j+1,k) - f(i,max(j-2,abs(j-3)),k))           &
                        + A(1)*(f(i,j,k) - f(i,max(j-1,abs(j-2)),k)))*g(j)
                    enddo

                    do j=n2-4,n2-1
                        d1f(i,j,k)=   (A(5)*(f(i,min(j+4, 2*n2-(j+5)),k) - f(i,j-5,k))     &
                        + A(4)*(f(i,min(j+3, 2*n2-(j+4)),k) - f(i,j-4,k))             &
                        + A(3)*(f(i,min(j+2, 2*n2-(j+3)),k) - f(i,j-3,k))             &
                        + A(2)*(f(i,min(j+1, 2*n2-(j+2)),k) - f(i,j-2,k))             &
                        + A(1)*(f(i,j,k) - f(i,j-1,k)))*g(j)
                    enddo

                enddo
            enddo

        end if


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e
                    d1f(i,2,k)= ((f(i,2,k)-f(i,1,k))                                 /h)*g(2)
                    d1f(i,3,k)  = (9.0d0/8.d0  * (f(i,3,k)-f(i,2,k))                 / h         &
                    -   1.d0/24.d0  * (f(i,4,k)-f(i,1,k))                           / h)*g(3)
                    d1f(i,4,k)  = (1.189101951200031d0  *   (f(i,4,k)-f(i,3,k))      / h         &
                    -   0.073717642266682d0  *   (f(i,5,k)-f(i,2,k))                / h         &
                    +   0.006410195120003d0 *   (f(i,6,k) - f(i,1,k))               / h)*g(4)
                    d1f(i,5,k)  = (1.189101951200031d0  *   (f(i,5,k)-f(i,4,k))      / h         &
                    -   0.073717642266682d0  *   (f(i,6,k)-f(i,3,k))                / h         &
                    +   0.006410195120003d0 *   (f(i,7,k) - f(i,2,k))               / h)*g(5)

                enddo
            enddo

            do k=1, n3e
                do i=1, n1e
                    d1f(i,n2-1,k)=((f(i,n2-1,k)-f(i,n2-2,k))                                 /h)*g(n2-1)
                    d1f(i,n2-2,k)  =   (9.0d0/8.d0 * (f(i,n2-2,k)-f(i,n2-3,k))               /h      &
                    -   1.d0/24.d0    * (f(i,n2-1,k)-f(i,n2-4,k))                           /h)*g(n2-2)
                    d1f(i,n2-3,k)  =   (1.189101951200031d0  *   (f(i,n2-3,k)-f(i,n2-4,k))   /h      &
                    -   0.073717642266682d0  *   (f(i,n2-2,k)-f(i,n2-5,k))                  /h      &
                    +   0.006410195120003d0 *   (f(i,n2-1,k) - f(i,n2-6,k))                 /h)*g(n2-3)
                    d1f(i,n2-4,k)  =   (1.189101951200031d0  *   (f(i,n2-4,k)-f(i,n2-5,k))   /h      &
                    -   0.073717642266682d0  *   (f(i,n2-3,k)-f(i,n2-6,k))                  /h      &
                    +   0.006410195120003d0 *   (f(i,n2-2,k) - f(i,n2-7,k))                 /h)*g(n2-4)

                enddo
            enddo

        endif

        return

    end subroutine D1ssh_DRP5_MULT_3Dy

    subroutine D1ssh_DRP5_MULTACC_3Dy(f, d1f, n1,n1e,n2,n3,n3e, h, shifted, bc_type, g)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)       :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)      :: d1f
        real(kind=8), dimension(:), intent(in)              :: g

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        do k=1, n3e
            do i=1, n1e
                do j=5,n2-6
                    d1f(i, j+1, k)= d1f(i, j+1, k)+(A(5)*(f(i, j+5, k) - f(i, j-4, k))     &
                    + A(4)*(f(i, j+4, k) - f(i, j-3, k))     &
                    + A(3)*(f(i, j+3, k) - f(i, j-2, k))     &
                    + A(2)*(f(i, j+2, k) - f(i, j-1, k))     &
                    + A(1)*(f(i, j+1, k) - f(i, j, k)))*g(j+1)
                enddo
            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e

                    do j=1, 5
                        d1f(i, j, k)=   d1f(i, j, k)+(A(5)*(f(i, j+4, k) - f(i, mod(n2+j-7,n2-1)+1, k))  &
                        + A(4)*(f(i, j+3, k) - f(i, mod(n2+j-6,n2-1)+1, k))  &
                        + A(3)*(f(i, j+2, k) - f(i, mod(n2+j-5,n2-1)+1, k))  &
                        + A(2)*(f(i, j+1, k) - f(i, mod(n2+j-4,n2-1)+1, k))  &
                        + A(1)*(f(i, j, k) - f(i, mod(n2+j-3,n2-1)+1, k)))*g(j)
                    enddo

                    do j=n2-4,n2-1
                        d1f(i, j, k)= d1f(i, j, k)+(A(5)*(f(i, mod(j+3,n2-1)+1, k) - f(i, j-5, k))             &
                        + A(4)*(f(i, mod(j+2,n2-1)+1, k) - f(i, j-4, k))             &
                        + A(3)*(f(i, mod(j+1,n2-1)+1, k) - f(i, j-3, k))             &
                        + A(2)*(f(i, mod(j,n2-1)+1, k) - f(i, j-2, k))             &
                        + A(1)*(f(i, j, k) - f(i, j-1, k)))*g(j)
                    enddo

                enddo
            enddo

        endif

        if (bc_type==symetric) then
            do k=1, n3e
                do i=1, n1e

                    do j=1, 5
                        d1f(i,j,k)=   d1f(i, j, k)+(A(5)*(f(i,j+4,k) - f(i,max(j-5,abs(j-6)),k))   &
                        + A(4)*(f(i,j+3,k) - f(i,max(j-4,abs(j-5)),k))           &
                        + A(3)*(f(i,j+2,k) - f(i,max(j-3,abs(j-4)),k))           &
                        + A(2)*(f(i,j+1,k) - f(i,max(j-2,abs(j-3)),k))           &
                        + A(1)*(f(i,j,k) - f(i,max(j-1,abs(j-2)),k)))*g(j)
                    enddo

                    do j=n2-4,n2-1
                        d1f(i,j,k)=   d1f(i, j, k)+(A(5)*(f(i,min(j+4, 2*n2-(j+5)),k) - f(i,j-5,k))     &
                        + A(4)*(f(i,min(j+3, 2*n2-(j+4)),k) - f(i,j-4,k))             &
                        + A(3)*(f(i,min(j+2, 2*n2-(j+3)),k) - f(i,j-3,k))             &
                        + A(2)*(f(i,min(j+1, 2*n2-(j+2)),k) - f(i,j-2,k))             &
                        + A(1)*(f(i,j,k) - f(i,j-1,k)))*g(j)
                    enddo

                enddo
            enddo

        end if


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e
                    d1f(i,2,k)= d1f(i,2,k)+((f(i,2,k)-f(i,1,k))                                 /h)*g(2)
                    d1f(i,3,k)  = d1f(i,3,k)+(9.0d0/8.d0  * (f(i,3,k)-f(i,2,k))                 / h         &
                    -   1.d0/24.d0  * (f(i,4,k)-f(i,1,k))                           / h)*g(3)
                    d1f(i,4,k)  = d1f(i,4,k)+(1.189101951200031d0  *   (f(i,4,k)-f(i,3,k))      / h         &
                    -   0.073717642266682d0  *   (f(i,5,k)-f(i,2,k))                / h         &
                    +   0.006410195120003d0 *   (f(i,6,k) - f(i,1,k))               / h)*g(4)
                    d1f(i,5,k)  = d1f(i,5,k)+(1.189101951200031d0  *   (f(i,5,k)-f(i,4,k))      / h         &
                    -   0.073717642266682d0  *   (f(i,6,k)-f(i,3,k))                / h         &
                    +   0.006410195120003d0 *   (f(i,7,k) - f(i,2,k))               / h)*g(5)

                enddo
            enddo

            do k=1, n3e
                do i=1, n1e
                    d1f(i,n2-1,k)=d1f(i,n2-1,k)+((f(i,n2-1,k)-f(i,n2-2,k))                                 /h)*g(n2-1)
                    d1f(i,n2-2,k)  =   d1f(i,n2-2,k)+(9.0d0/8.d0 * (f(i,n2-2,k)-f(i,n2-3,k))               /h      &
                    -   1.d0/24.d0    * (f(i,n2-1,k)-f(i,n2-4,k))                           /h)*g(n2-2)
                    d1f(i,n2-3,k)  =   d1f(i,n2-3,k)+(1.189101951200031d0  *   (f(i,n2-3,k)-f(i,n2-4,k))   /h      &
                    -   0.073717642266682d0  *   (f(i,n2-2,k)-f(i,n2-5,k))                  /h      &
                    +   0.006410195120003d0 *   (f(i,n2-1,k) - f(i,n2-6,k))                 /h)*g(n2-3)
                    d1f(i,n2-4,k)  =   d1f(i,n2-4,k)+(1.189101951200031d0  *   (f(i,n2-4,k)-f(i,n2-5,k))   /h      &
                    -   0.073717642266682d0  *   (f(i,n2-3,k)-f(i,n2-6,k))                  /h      &
                    +   0.006410195120003d0 *   (f(i,n2-2,k) - f(i,n2-7,k))                 /h)*g(n2-4)

                enddo
            enddo

        endif

        return

    end subroutine D1ssh_DRP5_MULTACC_3Dy

    subroutine D1ssh_DRP5_3Dz(f, d1f, n1,n1e,n2,n2e,n3, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: d1f

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        do j=1, n2e
            do i=1, n1e
                do k=5,n3-6
                    d1f(i, j, k+1)= A(5)*(f(i, j, k+5) - f(i, j, k-4))     &
                    + A(4)*(f(i, j, k+4) - f(i, j, k-3))     &
                    + A(3)*(f(i, j, k+3) - f(i, j, k-2))     &
                    + A(2)*(f(i, j, k+2) - f(i, j, k-1))     &
                    + A(1)*(f(i, j, k+1) - f(i, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e

                    do k=1, 5
                        d1f(i, j, k)=   A(5)*(f(i, j, k+4) - f(i, j, mod(n3+k-7,n3-1)+1))  &
                        + A(4)*(f(i, j, k+3) - f(i, j, mod(n3+k-6,n3-1)+1))  &
                        + A(3)*(f(i, j, k+2) - f(i, j, mod(n3+k-5,n3-1)+1))  &
                        + A(2)*(f(i, j, k+1) - f(i, j, mod(n3+k-4,n3-1)+1))  &
                        + A(1)*(f(i, j, k) - f(i, j, mod(n3+k-3,n3-1)+1))
                    enddo

                    do k=n3-4,n3-1
                        d1f(i, j, k)= A(5)*(f(i, j, mod(k+3,n3-1)+1) - f(i, j, k-5))             &
                        + A(4)*(f(i, j, mod(k+2,n3-1)+1) - f(i, j, k-4))             &
                        + A(3)*(f(i, j, mod(k+1,n3-1)+1) - f(i, j, k-3))             &
                        + A(2)*(f(i, j, mod(k,n3-1)+1) - f(i, j, k-2))             &
                        + A(1)*(f(i, j, k) - f(i, j, k-1))
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,2)= (f(i,j,2)-f(i,j,1))                                 /h
                    d1f(i,j,3)  = 9.0d0/8.d0  * (f(i,j,3)-f(i,j,2))                 / h         &
                    -   1.d0/24.d0  * (f(i,j,4)-f(i,j,1))                           / h
                    d1f(i,j,4)  = 1.189101951200031d0  *   (f(i,j,4)-f(i,j,3))      / h         &
                    -   0.073717642266682d0  *   (f(i,j,5)-f(i,j,2))                / h         &
                    +   0.006410195120003d0 *   (f(i,j,6) - f(i,j,1))               / h
                    d1f(i,j,5)  = 1.189101951200031d0  *   (f(i,j,5)-f(i,j,4))      / h         &
                    -   0.073717642266682d0  *   (f(i,j,6)-f(i,j,3))                / h         &
                    +   0.006410195120003d0 *   (f(i,j,7) - f(i,j,2))               / h
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,n3-1)=(f(i,j,n3-1)-f(i,j,n3-2))                                 /h
                    d1f(i,j,n3-2)  =   9.0d0/8.d0 * (f(i,j,n3-2)-f(i,j,n3-3))               /h      &
                    -   1.d0/24.d0    * (f(i,j,n3-1)-f(i,j,n3-4))                           /h
                    d1f(i,j,n3-3)  =   1.189101951200031d0  *   (f(i,j,n3-3)-f(i,j,n3-4))   /h      &
                    -   0.073717642266682d0  *   (f(i,j,n3-2)-f(i,j,n3-5))                  /h      &
                    +   0.006410195120003d0 *   (f(i,j,n3-1) - f(i,j,n3-6))                 /h
                    d1f(i,j,n3-4)  =   1.189101951200031d0  *   (f(i,j,n3-4)-f(i,j,n3-5))   /h      &
                    -   0.073717642266682d0  *   (f(i,j,n3-3)-f(i,j,n3-6))                  /h      &
                    +   0.006410195120003d0 *   (f(i,j,n3-2) - f(i,j,n3-7))                 /h
                enddo
            enddo

        endif

        return

    end subroutine D1ssh_DRP5_3Dz

    subroutine D1ssh_DRP5_ACC_3Dz(f, d1f, n1,n1e,n2,n2e,n3, h, shifted, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        logical, intent(in)     :: shifted
        real*8, intent(in)      :: h
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: d1f

        integer :: i,j,k
        real(kind=8) A(5)

        A(1) =  1.229770180607572d0      /   (h)
        A(2) =  -0.30840585496884d0     /   (3.d0*h)
        A(3) =  0.10067299049373d0      /   (5.d0*h)
        A(4) =  -0.025440011202968d0    /   (7.d0*h)
        A(5) =  0.0034026950705542d0    /   (9.d0*h)

        do j=1, n2e
            do i=1, n1e
                do k=5,n3-6
                    d1f(i, j, k+1)= d1f(i, j, k+1)+A(5)*(f(i, j, k+5) - f(i, j, k-4))     &
                    + A(4)*(f(i, j, k+4) - f(i, j, k-3))     &
                    + A(3)*(f(i, j, k+3) - f(i, j, k-2))     &
                    + A(2)*(f(i, j, k+2) - f(i, j, k-1))     &
                    + A(1)*(f(i, j, k+1) - f(i, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e

                    do k=1, 5
                        d1f(i, j, k)=   d1f(i, j, k)+A(5)*(f(i, j, k+4) - f(i, j, mod(n3+k-7,n3-1)+1))  &
                        + A(4)*(f(i, j, k+3) - f(i, j, mod(n3+k-6,n3-1)+1))  &
                        + A(3)*(f(i, j, k+2) - f(i, j, mod(n3+k-5,n3-1)+1))  &
                        + A(2)*(f(i, j, k+1) - f(i, j, mod(n3+k-4,n3-1)+1))  &
                        + A(1)*(f(i, j, k) - f(i, j, mod(n3+k-3,n3-1)+1))
                    enddo

                    do k=n3-4,n3-1
                        d1f(i, j, k)= d1f(i, j, k)+A(5)*(f(i, j, mod(k+3,n3-1)+1) - f(i, j, k-5))             &
                        + A(4)*(f(i, j, mod(k+2,n3-1)+1) - f(i, j, k-4))             &
                        + A(3)*(f(i, j, mod(k+1,n3-1)+1) - f(i, j, k-3))             &
                        + A(2)*(f(i, j, mod(k,n3-1)+1) - f(i, j, k-2))             &
                        + A(1)*(f(i, j, k) - f(i, j, k-1))
                    enddo

                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,2)= d1f(i,j,2)+(f(i,j,2)-f(i,j,1))                                 /h
                    d1f(i,j,3)  = d1f(i,j,3)+9.0d0/8.d0  * (f(i,j,3)-f(i,j,2))                 / h         &
                    -   1.d0/24.d0  * (f(i,j,4)-f(i,j,1))                           / h
                    d1f(i,j,4)  = d1f(i,j,4)+1.189101951200031d0  *   (f(i,j,4)-f(i,j,3))      / h         &
                    -   0.073717642266682d0  *   (f(i,j,5)-f(i,j,2))                / h         &
                    +   0.006410195120003d0 *   (f(i,j,6) - f(i,j,1))               / h
                    d1f(i,j,5)  = d1f(i,j,5)+1.189101951200031d0  *   (f(i,j,5)-f(i,j,4))      / h         &
                    -   0.073717642266682d0  *   (f(i,j,6)-f(i,j,3))                / h         &
                    +   0.006410195120003d0 *   (f(i,j,7) - f(i,j,2))               / h
                enddo
            enddo

            do j=1, n2e
                do i=1, n1e
                    d1f(i,j,n3-1)=d1f(i,j,n3-1)+(f(i,j,n3-1)-f(i,j,n3-2))                                 /h
                    d1f(i,j,n3-2)  =   d1f(i,j,n3-2)+9.0d0/8.d0 * (f(i,j,n3-2)-f(i,j,n3-3))               /h      &
                    -   1.d0/24.d0    * (f(i,j,n3-1)-f(i,j,n3-4))                           /h
                    d1f(i,j,n3-3)  =   d1f(i,j,n3-3)+1.189101951200031d0  *   (f(i,j,n3-3)-f(i,j,n3-4))   /h      &
                    -   0.073717642266682d0  *   (f(i,j,n3-2)-f(i,j,n3-5))                  /h      &
                    +   0.006410195120003d0 *   (f(i,j,n3-1) - f(i,j,n3-6))                 /h
                    d1f(i,j,n3-4)  =   d1f(i,j,n3-4)+1.189101951200031d0  *   (f(i,j,n3-4)-f(i,j,n3-5))   /h      &
                    -   0.073717642266682d0  *   (f(i,j,n3-3)-f(i,j,n3-6))                  /h      &
                    +   0.006410195120003d0 *   (f(i,j,n3-2) - f(i,j,n3-7))                 /h
                enddo
            enddo

        endif

        return

    end subroutine D1ssh_DRP5_ACC_3Dz

end module DRP_D1s_sh
