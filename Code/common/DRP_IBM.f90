module DRP_IBM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module that contains all schemes used with the DRP scheme and the IBM
! Basically, this is only the O2 scheme that is slightly modified
! to compute derivatives on some sections of the domain

! Assumptions :
! - Periodic conditions in X and Z directions
! - Dirichlet condition in Y direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!     D2c Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine D2c_IBM_O2_3Dx(f,d2f,n1size,n2size,n3size,h)

        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, intent(in)                                      :: h
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: d2f

        integer :: i,j,k
        real*8 A(0:1)

        A(1) =  1.d0 / h**2
        A(0) = -2.d0 / h**2

        do i=2, n1size-1
            do j=1, n2size
                do k=1, n3size

                    d2f(i, j, k)= A(1)*(f(i+1, j, k) + f(i-1, j, k))     &
                    + A(0)*f(i, j, k)

                enddo
            enddo
        enddo

        return

    end subroutine D2c_IBM_O2_3Dx



    subroutine D2c_IBM_O2_ACC_3Dx(f,d2f,n1size,n2size,n3size,h)

        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, intent(in)                                      :: h
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: d2f

        integer :: i,j,k
        real*8 A(0:1)

        A(1) =  1.d0 / h**2
        A(0) = -2.d0 / h**2

        do i=2, n1size-1
            do j=1, n2size
                do k=1, n3size

                    d2f(i, j, k)= d2f(i, j, k) + A(1)*(f(i+1, j, k) + f(i-1, j, k))     &
                    + A(0)*f(i, j, k)

                enddo
            enddo
        enddo

        return

    end subroutine D2c_IBM_O2_ACC_3Dx



    subroutine D2c_IBM_O2_MULTACC_3Dy(f,d2f,n1size,n2size,n3size,h,g)

        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, intent(in)                                      :: h
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: d2f
        real*8, dimension(:), intent(in)                        :: g

        integer :: i,j,k
        real*8 A(0:1)

        A(1) =  1.d0 / h**2
        A(0) = -2.d0*(A(1))

        do i=1, n1size
            do j=2, n2size-1
                do k=1, n3size

                    d2f(i,j,k)= d2f(i,j,k)+(A(1)*(f(i,j+1,k) + f(i,j-1,k))     &
                    + A(0)*f(i,j,k))*g(j)

                enddo

            enddo
        enddo

        return

    end subroutine D2c_IBM_O2_MULTACC_3Dy



    subroutine D2c_IBM_O2_3Dz(f,d2f,n1size,n2size,n3size,h)

        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, intent(in)                                      :: h
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: d2f

        integer :: i,j,k
        real*8 A(0:1)

        A(1) =  1.d0 / h**2
        A(0) = -2.d0 / h**2

        do i=1, n1size
            do j=1, n2size
                do k=2, n3size-1

                    d2f(i,j,k)= A(1)*(f(i,j,k+1) + f(i,j,k-1))        &
                    + A(0)*f(i,j,k)

                enddo
            enddo
        enddo

        return

    end subroutine D2c_IBM_O2_3Dz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!     D1c Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine D1c_IBM_O2_3Dx(f, d1f, n1size,n2size,n3size,h)

        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, intent(in)                                      :: h
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: d1f

        integer :: i,j,k
        real*8 A

        A =  0.5d0    /h

        do i=2, n1size-1
            do j=1, n2size
                do k=1, n3size

                    d1f(i,j,k)= A*(f(i+1,j,k) - f(i-1,j,k))

                enddo
            enddo
        enddo

        return

    end subroutine D1c_IBM_O2_3Dx

    subroutine D1c_IBM_O2_ACC_3Dx(f, d1f, n1size,n2size,n3size,h)

        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, intent(in)                                      :: h
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: d1f

        integer :: i,j,k
        real*8 A

        A =  0.5d0    /h

        do i=2, n1size-1
            do j=1, n2size
                do k=1, n3size

                    d1f(i,j,k)= d1f(i,j,k) + A*(f(i+1,j,k) - f(i-1,j,k))

                enddo
            enddo
        enddo

        return

    end subroutine D1c_IBM_O2_ACC_3Dx



    subroutine D1c_IBM_O2_MULT_3Dy(f,d1f,n1size,n2size,n3size,h,g)

        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, intent(in)                                      :: h
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: d1f
        real*8, dimension(:), intent(in)                        :: g

        integer :: i,j,k
        real*8 A

        A =  0.5d0    /h

        do i=1, n1size
            do j=2, n2size-1
                do k=1, n3size

                    d1f(i,j,k)= (A*(f(i,j+1,k) - f(i,j-1,k)))*g(j)

                enddo
            enddo
        enddo

        return

    end subroutine D1c_IBM_O2_MULT_3Dy



    subroutine D1c_IBM_O2_MULTACC_3Dy(f,d1f,n1size,n2size,n3size,h,g)

        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, intent(in)                                      :: h
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: d1f
        real*8, dimension(:), intent(in)                        :: g

        integer :: i,j,k
        real*8 A

        A =  0.5d0    /h

        do i=1, n1size
            do j=2, n2size-1
                do k=1, n3size

                    d1f(i,j,k) = d1f(i,j,k) + (A*(f(i,j+1,k) - f(i,j-1,k)))*g(j)

                enddo
            enddo
        enddo

        return

    end subroutine D1c_IBM_O2_MULTACC_3Dy



    subroutine D1c_IBM_O2_3Dz(f, d1f, n1size,n2size,n3size,h)

        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, intent(in)                                      :: h
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: d1f

        integer :: i,j,k
        real*8 A

        A =  0.5d0    /h

        do i=1, n1size
            do j=1, n2size
                do k=2, n3size-1

                    d1f(i,j,k)= A*(f(i,j,k+1) - f(i,j,k-1))

                enddo
            enddo
        enddo

        return

    end subroutine D1c_IBM_O2_3Dz

    subroutine D1c_IBM_O2_ACC_3Dz(f, d1f, n1size,n2size,n3size,h)

        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, intent(in)                                      :: h
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: d1f

        integer :: i,j,k
        real*8 A

        A =  0.5d0    /h

        do i=1, n1size
            do j=1, n2size
                do k=2, n3size-1

                    d1f(i,j,k)= d1f(i,j,k) + A*(f(i,j,k+1) - f(i,j,k-1))

                enddo
            enddo
        enddo

        return

    end subroutine D1c_IBM_O2_ACC_3Dz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!     D0ssh Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine D0ssh_IBM_O2_3Dx(f,ff,n1size,n2size,n3size)
        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: ff

        integer :: i,j,k

        do i=2, n1size
            do j=1, n2size
                do k=1, n3size

                    ff(i, j, k)= 0.5d0*(f(i, j, k) + f(i-1, j, k))

                enddo
            enddo
        enddo

        return

    end subroutine D0ssh_IBM_O2_3Dx



    subroutine D0ssh_IBM_O2_3Dy(f,ff,n1size,n2size,n3size)
        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: ff

        integer :: i,j,k

        do i=1, n1size
            do j=2, n2size
                do k=1, n3size

                    ff(i, j, k)= 0.5d0*(f(i, j, k) + f(i, j-1, k))

                enddo
            enddo
        enddo

        return

    end subroutine D0ssh_IBM_O2_3Dy



    subroutine D0ssh_IBM_O2_3Dz(f,ff,n1size,n2size,n3size)
        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: ff

        integer :: i,j,k

        do i=1, n1size
            do j=1, n2size
                do k=2, n3size

                    ff(i, j, k)= 0.5d0*(f(i, j, k) + f(i, j, k-1))

                enddo
            enddo
        enddo

        return

    end subroutine D0ssh_IBM_O2_3Dz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!     D0s Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine D0s_IBM_O2_3Dx(f,ff,n1size,n2size,n3size)
        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: ff

        integer :: i,j,k

        do i=1, n1size-1
            do j=1, n2size
                do k=1, n3size

                    ff(i, j, k)= 0.5d0*(f(i+1, j, k) + f(i, j, k))

                enddo
            enddo
        enddo

        return

    end subroutine D0s_IBM_O2_3Dx



    subroutine D0s_IBM_O2_3Dy(f,ff,n1size,n2size,n3size)
        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: ff

        integer :: i,j,k

        do i=1, n1size
            do j=1, n2size-1
                do k=1, n3size

                    ff(i, j, k)= 0.5d0*(f(i, j+1, k) + f(i, j, k))

                enddo
            enddo
        enddo

        return

    end subroutine D0s_IBM_O2_3Dy




    subroutine D0s_IBM_O2_3Dz(f,ff,n1size,n2size,n3size)
        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: ff

        integer :: i,j,k

        do i=1, n1size
            do j=1, n2size
                do k=1, n3size-1

                    ff(i, j, k)= 0.5d0*(f(i, j, k+1) + f(i, j, k))

                enddo
            enddo
        enddo

        return

    end subroutine D0s_IBM_O2_3Dz



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!     D1ssh Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine D1ssh_IBM_O2_ACC_3Dx(f,d1f,n1size,n2size,n3size,h)
        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, intent(in)                                      :: h
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: d1f

        integer :: i,j,k
        real(kind=8) A

        A =  1.0d0      /   (h)

        do i=2, n1size
            do j=1, n2size
                do k=1, n3size

                    d1f(i, j, k)= d1f(i, j, k)+A*(f(i, j, k) - f(i-1, j, k))

                enddo
            enddo
        enddo

        return

    end subroutine D1ssh_IBM_O2_ACC_3Dx



    subroutine D1ssh_IBM_O2_MULTACC_3Dy(f,d1f,n1size,n2size,n3size,h,g)
        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, intent(in)                                      :: h
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: d1f
        real*8, dimension(:), intent(in)                        :: g

        integer :: i,j,k
        real(kind=8) A

        A =  1.0d0      /   (h)

        do i=1, n1size
            do j=2, n2size
                do k=1, n3size

                    d1f(i, j, k)= d1f(i, j, k)+(A*(f(i, j, k) - f(i, j-1, k)))*g(j)

                enddo
            enddo
        enddo

        return

    end subroutine D1ssh_IBM_O2_MULTACC_3Dy



    subroutine D1ssh_IBM_O2_ACC_3Dz(f,d1f,n1size,n2size,n3size,h)
        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, intent(in)                                      :: h
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: d1f

        integer :: i,j,k
        real(kind=8) A

        A =  1.0d0      /   (h)

        do i=1, n1size
            do j=1, n2size
                do k=2, n3size

                    d1f(i, j, k)= A*(f(i, j, k) - f(i, j, k-1))

                enddo
            enddo
        enddo

        return

    end subroutine D1ssh_IBM_O2_ACC_3Dz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!     D1s Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine D1s_IBM_O2_ACC_3Dx(f,d1f,n1size,n2size,n3size,h)
        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: d1f
        real*8, intent(in)                                      :: h

        integer :: i,j,k
        real*8 A

        A =  1.0d0      /   (h)

        do i=1, n1size-1
            do j=1, n2size
                do k=1, n3size

                    d1f(i, j, k)= d1f(i, j, k) + A*(f(i+1, j, k) - f(i, j, k))

                enddo
            enddo
        enddo

        return

    end subroutine D1s_IBM_O2_ACC_3Dx



    subroutine D1s_IBM_O2_MULTACC_3Dy(f,d1f,n1size,n2size,n3size,h,g)
        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: d1f
        real*8, intent(in)                                      :: h
        real*8, dimension(:), intent(in)                        :: g

        integer :: i,j,k
        real*8 A

        A =  1.0d0      /   (h)

        do i=1, n1size
            do j=1, n2size-1
                do k=1, n3size

                    d1f(i, j, k)= d1f(i, j, k)+(A*(f(i, j+1, k) - f(i, j, k)))*g(j)

                enddo
            enddo
        enddo

        return

    end subroutine D1s_IBM_O2_MULTACC_3Dy



    subroutine D1s_IBM_O2_ACC_3Dz(f,d1f,n1size,n2size,n3size,h)
        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: d1f
        real*8, intent(in)                                      :: h

        integer :: i,j,k
        real*8 A

        A =  1.0d0      /   (h)

        do i=1, n1size
            do j=1, n2size
                do k=1, n3size-1

                    d1f(i, j, k)= d1f(i, j, k)+A*(f(i, j, k+1) - f(i, j, k))

                enddo
            enddo
        enddo

        return

    end subroutine D1s_IBM_O2_ACC_3Dz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!     shift Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Forward_shift_IBM_3Dx(f,ff,n1size,n2size,n3size)
        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: ff

        integer :: i,j,k

        do i=1, n1size-1
            do j=1, n2size
                do k=1, n3size

                    ff(i, j, k)= f(i+1, j, k)

                enddo
            enddo
        enddo

        return

    end subroutine Forward_shift_IBM_3Dx

    subroutine Forward_shift_IBM_3Dy(f,ff,n1size,n2size,n3size)
        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: ff

        integer :: i,j,k

        do i=1, n1size
            do j=1, n2size-1
                do k=1, n3size

                    ff(i, j, k)= f(i, j+1, k)

                enddo
            enddo
        enddo

        return

    end subroutine Forward_shift_IBM_3Dy

    subroutine Forward_shift_IBM_3Dz(f,ff,n1size,n2size,n3size)
        implicit none

        integer, intent(in)                                     :: n1size,n2size,n3size
        real*8, dimension(:,:,:), intent(in)    :: f
        real*8, dimension(:,:,:), intent(out)   :: ff

        integer :: i,j,k

        do i=1, n1size
            do j=1, n2size
                do k=1, n3size-1

                    ff(i, j, k)= f(i, j, k+1)

                enddo
            enddo
        enddo

        return

    end subroutine Forward_shift_IBM_3Dz

end module DRP_IBM