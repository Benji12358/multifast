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

    subroutine D2c_IBM_O2_3Dx(f,d2f,n1,i_IBM_start,i_IBM_end,h)

        implicit none

        integer, intent(in) :: n1,i_IBM_start,i_IBM_end
        logical                 :: onBounds
        real*8, intent(in)      :: h
        real*8, dimension(n1), intent(in):: f
        real*8, dimension(n1), intent(out):: d2f

        integer :: i
        real(kind=8) A(0:1)

        A(1) =  1.d0 / h**2
        A(0) = -2.d0 / h**2

        do i = i_IBM_start,i_IBM_end

            onBounds = (i.le.1).or.(i.ge.n1-1)
        
            if (.not.onBounds) d2f(i) = A(1)*(f(i+1) + f(i-1)) + A(0)*f(i)

            if (i.eq.1) d2f(1) = A(1)*(f(2) + f(n1-1)) + A(0)*f(1)

            if (i.eq.n1-1) d2f(n1-1) = A(1)*(f(1) + f(n1-2)) + A(0)*f(n1-1)

        enddo

        return

    end subroutine D2c_IBM_O2_3Dx



    subroutine D2c_IBM_O2_ACC_3Dx(f,d2f,n1,i_IBM_start,i_IBM_end,h)

        implicit none

        integer, intent(in) :: n1,i_IBM_start,i_IBM_end
        logical                 :: onBounds
        real*8, intent(in)      :: h
        real*8, dimension(n1), intent(in):: f
        real*8, dimension(n1), intent(out):: d2f

        integer :: i
        real(kind=8) A(0:1)

        A(1) =  1.d0 / h**2
        A(0) = -2.d0 / h**2

        do i = i_IBM_start,i_IBM_end

            onBounds = (i.le.1).or.(i.ge.n1-1)
        
            if (.not.onBounds) d2f(i) = d2f(i) + A(1)*(f(i+1) + f(i-1)) + A(0)*f(i)

            if (i.eq.1) d2f(1) = d2f(1) + A(1)*(f(2) + f(n1-1)) + A(0)*f(1)

            if (i.eq.n1-1) d2f(n1-1) = d2f(n1-1) + A(1)*(f(1) + f(n1-2)) + A(0)*f(n1-1)

        enddo

        return

    end subroutine D2c_IBM_O2_ACC_3Dx



    subroutine D2c_IBM_O2_MULTACC_3Dy(f,d2f,n2,j_IBM_start,j_IBM_end,h,g)

        implicit none

        integer, intent(in)     :: n2,j_IBM_start,j_IBM_end
        logical                 :: onBounds
        real*8, intent(in)      :: h
        real(kind=8), dimension(n2), intent(in)      :: f
        real(kind=8), dimension(n2), intent(out)     :: d2f
        real(kind=8), dimension(:), intent(in)       :: g

        integer :: j
        real(kind=8) A(0:1)

        A(1) =  1.d0 / h**2
        A(0) = -2.d0*(A(1))

        do j = j_IBM_start,j_IBM_end

            onBounds = (j.le.1).or.(j.ge.n2-1)
        
            if (.not.onBounds) d2f(j) = d2f(j) + (A(1)*(f(j+1) + f(j-1)) + A(0)*f(j))*g(j)

            if (j.eq.n2-1) d2f(n2-1) = d2f(n2-1)+( (-2.0d0*f(n2-1) + (f(n2)+f(n2-2)) ) / h**2)*g(n2-1)

        enddo

        return

    end subroutine D2c_IBM_O2_MULTACC_3Dy



    subroutine D2c_IBM_O2_3Dz(f,d2f,n3,k_IBM_start,k_IBM_end,h)

        implicit none

        integer, intent(in)     :: n3,k_IBM_start,k_IBM_end
        logical                 :: onBounds
        real*8, intent(in)      :: h
        real(kind=8), dimension(n3), intent(in)        :: f
        real(kind=8), dimension(n3), intent(out)       :: d2f

        integer :: k
        real(kind=8) A(0:1)

        A(1) =  1.d0 / h**2
        A(0) = -2.d0*(A(1))

        do k = k_IBM_start,k_IBM_end

            onBounds = (k.le.1).or.(k.ge.n3-1)
        
            if (.not.onBounds) d2f(k) = A(1)*(f(k+1) + f(k-1)) + A(0)*f(k)

            if (k.eq.1) d2f(1) = A(1)*(f(2) + f(n3-1)) + A(0)*f(1)

            if (k.eq.n3-1) d2f(n3-1) = A(1)*(f(1) + f(n3-2)) + A(0)*f(n3-1)

        enddo

        return

    end subroutine D2c_IBM_O2_3Dz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!     D1c Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine D1c_IBM_O2_3Dx(f, d1f, n1,i_IBM_start,i_IBM_end,h)

        implicit none

        integer, intent(in) :: n1,i_IBM_start,i_IBM_end
        logical                 :: onBounds
        real*8, intent(in)      :: h
        real*8, dimension(n1), intent(in):: f
        real*8, dimension(n1), intent(out):: d1f

        integer :: i
        real(kind=8) A

        A =  0.5d0    /h

        do i = i_IBM_start,i_IBM_end

            onBounds = (i.le.1).or.(i.ge.n1-1)

            if (.not.onBounds) d1f(i) = A * (f(i+1) - f(i-1))

            if (i.eq.1) d1f(1) = A * (f(2) - f(n1-1))

            if (i.eq.n1-1) d1f(n1-1) = A * (f(1) - f(n1-2))

        enddo

        return

    end subroutine D1c_IBM_O2_3Dx



    subroutine D1c_IBM_O2_MULT_3Dy(f,d1f,n2,j_IBM_start,j_IBM_end,h,g)

        implicit none

        integer, intent(in)     :: n2,j_IBM_start,j_IBM_end
        logical                 :: onBounds
        real*8, intent(in)      :: h
        real(kind=8), dimension(n2), intent(in)      :: f
        real(kind=8), dimension(n2), intent(out)     :: d1f
        real(kind=8), dimension(:), intent(in)       :: g

        integer :: j
        real(kind=8) A

        A =  0.5d0    /h

        do j = j_IBM_start,j_IBM_end

            onBounds = (j.le.1).or.(j.ge.n2-1)
        
            if (.not.onBounds) d1f(j) = (A*(f(j+1) - f(j-1)))*g(j)

            if (j.eq.n2-1) d1f(n2-1) = A*((f(n2)-f(n2-2)))*g(n2-1)

        enddo

        return

    end subroutine D1c_IBM_O2_MULT_3Dy



    subroutine D1c_IBM_O2_MULTACC_3Dy(f,d1f,n2,j_IBM_start,j_IBM_end,h,g)

        implicit none

        integer, intent(in)     :: n2,j_IBM_start,j_IBM_end
        logical                 :: onBounds
        real*8, intent(in)      :: h
        real(kind=8), dimension(n2), intent(in)      :: f
        real(kind=8), dimension(n2), intent(out)     :: d1f
        real(kind=8), dimension(:), intent(in)       :: g

        integer :: j
        real(kind=8) A

        A =  0.5d0    /h

        do j = j_IBM_start,j_IBM_end

            onBounds = (j.le.1).or.(j.ge.n2-1)
        
            if (.not.onBounds) d1f(j) = d1f(j) + (A*(f(j+1) - f(j-1)))*g(j)

            if (j.eq.n2-1) d1f(n2-1) = d1f(n2-1) + A*((f(n2)-f(n2-2)))*g(n2-1)

        enddo

        return

    end subroutine D1c_IBM_O2_MULTACC_3Dy



    subroutine D1c_IBM_O2_3Dz(f, d1f, n3,k_IBM_start,k_IBM_end,h)

        implicit none

        integer, intent(in) :: n3,k_IBM_start,k_IBM_end
        logical                 :: onBounds
        real*8, intent(in)      :: h
        real*8, dimension(n3), intent(in):: f
        real*8, dimension(n3), intent(out):: d1f

        integer :: k
        real(kind=8) A

        A =  0.5d0    /h

        do k = k_IBM_start,k_IBM_end

            onBounds = (k.le.1).or.(k.ge.n3-1)

            if (.not.onBounds) d1f(k) = A * (f(k+1) - f(k-1))

            if (k.eq.1) d1f(1) = A * (f(2) - f(n3-1))

            if (k.eq.n3-1) d1f(n3-1) = A * (f(1) - f(n3-2))

        enddo

        return

    end subroutine D1c_IBM_O2_3Dz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!     D0ssh Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine D0ssh_IBM_O2_3Dx(f,ff,n1,i_IBM_start,i_IBM_end)
        implicit none

        integer, intent(in) :: n1,i_IBM_start,i_IBM_end
        logical                 :: onBounds
        real*8, dimension(n1), intent(in):: f
        real*8, dimension(n1), intent(out):: ff

        integer :: i

        do i = i_IBM_start,i_IBM_end

            onBounds = (i.le.1).or.(i.ge.n1)
        
            if (.not.onBounds) ff(i)= 0.5d0*(f(i) + f(i-1))

            if (i.eq.1) ff(1)=   0.5d0*(f(1) + f(n1-1))

        enddo

        return

    end subroutine D0ssh_IBM_O2_3Dx



    subroutine D0ssh_IBM_O2_3Dy(f,ff,n2,j_IBM_start,j_IBM_end)

        implicit none

        integer, intent(in)     :: n2,j_IBM_start,j_IBM_end
        logical                 :: onBounds
        real(kind=8), dimension(n2), intent(in)      :: f
        real(kind=8), dimension(n2), intent(out)     :: ff

        integer :: j
        real(kind=8) A

        A =  1.d0     /2.d0

        do j = j_IBM_start,j_IBM_end

            onBounds = (j.le.1).or.(j.ge.n2)
        
            if (.not.onBounds) ff(j)= A*(f(j) + f(j-1))

        enddo

        return

    end subroutine D0ssh_IBM_O2_3Dy



    subroutine D0ssh_IBM_O2_3Dz(f,ff,n3,k_IBM_start,k_IBM_end)

        implicit none

        integer, intent(in)     :: n3,k_IBM_start,k_IBM_end
        logical                 :: onBounds
        real(kind=8), dimension(n3), intent(in)        :: f
        real(kind=8), dimension(n3), intent(out)       :: ff

        integer :: k
        real(kind=8) A

        A =  1.d0     /2.d0

        do k = k_IBM_start,k_IBM_end

            onBounds = (k.le.1).or.(k.ge.n3)
        
            if (.not.onBounds) ff(k)= A*(f(k) + f(k-1))

            if (k.eq.1) ff(1)=   0.5d0*(f(1) + f(n3-1))

        enddo

        return

    end subroutine D0ssh_IBM_O2_3Dz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!     D0s Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine D0s_IBM_O2_3Dx(f,ff,n1,i_IBM_start,i_IBM_end)
        implicit none

        integer, intent(in) :: n1,i_IBM_start,i_IBM_end
        logical                 :: onBounds
        real*8, dimension(n1), intent(in):: f
        real*8, dimension(n1), intent(out):: ff

        integer :: i

        do i = i_IBM_start,i_IBM_end

            onBounds = (i.le.0).or.(i.ge.n1-1)
        
            if (.not.onBounds) ff(i)= 0.5d0*(f(i+1) + f(i))

            if (i.eq.n1-1) ff(n1-1) = 0.5d0*(f(1) + f(n1-1))

        enddo

        return

    end subroutine D0s_IBM_O2_3Dx



    subroutine D0s_IBM_O2_3Dy(f,ff,n2,j_IBM_start,j_IBM_end)

        implicit none

        integer, intent(in)     :: n2,j_IBM_start,j_IBM_end
        logical                 :: onBounds
        real(kind=8), dimension(n2), intent(in)      :: f
        real(kind=8), dimension(n2), intent(out)     :: ff

        integer :: j
        real(kind=8) A

        A =  1.d0     /2.d0

        do j = j_IBM_start,j_IBM_end

            onBounds = (j.le.0).or.(j.ge.n2-1)
        
            if (.not.onBounds) ff(j)= A*(f(j+1) + f(j))

            if (j.eq.n2-1) ff(n2-1) = 0.5d0*(f(n2) + f(n2-1))

        enddo

        return

    end subroutine D0s_IBM_O2_3Dy



    subroutine D0s_IBM_O2_3Dz(f,ff,n3,k_IBM_start,k_IBM_end)

        implicit none

        integer, intent(in)     :: n3,k_IBM_start,k_IBM_end
        logical                 :: onBounds
        real(kind=8), dimension(n3), intent(in)        :: f
        real(kind=8), dimension(n3), intent(out)       :: ff

        integer :: k
        real(kind=8) A

        A =  1.d0     /2.d0

        do k = k_IBM_start,k_IBM_end

            onBounds = (k.le.0).or.(k.ge.n3-1)
        
            if (.not.onBounds) ff(k)= A*(f(k) + f(k+1))

            if (k.eq.n3-1) ff(n3-1)=   0.5d0*(f(1) + f(n3-1))

        enddo

        return

    end subroutine D0s_IBM_O2_3Dz



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!     D1ssh Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine D1ssh_IBM_O2_ACC_3Dx(f,d1f,n1,i_IBM_start,i_IBM_end,h)
        implicit none

        integer, intent(in) :: n1,i_IBM_start,i_IBM_end
        logical                 :: onBounds
        real*8, dimension(n1), intent(in):: f
        real*8, dimension(n1), intent(out):: d1f
        real*8, intent(in)      :: h

        integer :: i
        real(kind=8) A

        A =  1.0d0      /   (h)

        do i = i_IBM_start,i_IBM_end

            onBounds = (i.le.1).or.(i.ge.n1)
        
            if (.not.onBounds) d1f(i) = d1f(i) + A*(f(i) - f(i-1))

            if (i.eq.1) d1f(1) = d1f(1) + A*(f(1) - f(n1-1))

        enddo

        return

    end subroutine D1ssh_IBM_O2_ACC_3Dx



    subroutine D1ssh_IBM_O2_MULTACC_3Dy(f,d1f,n2,j_IBM_start,j_IBM_end,h,g)
        implicit none

        integer, intent(in) :: n2,j_IBM_start,j_IBM_end
        logical                 :: onBounds
        real*8, dimension(n2), intent(in):: f
        real*8, dimension(n2), intent(out):: d1f
        real*8, intent(in)      :: h
        real(kind=8), dimension(:), intent(in)              :: g

        integer :: j
        real(kind=8) A

        A =  1.0d0      /   (h)

        do j = j_IBM_start,j_IBM_end

            onBounds = (j.le.1).or.(j.ge.n2)
        
            if (.not.onBounds) d1f(j) = d1f(j) + (A*(f(j) - f(j-1)))*g(j)

        enddo

        return

    end subroutine D1ssh_IBM_O2_MULTACC_3Dy



    subroutine D1ssh_IBM_O2_ACC_3Dz(f,d1f,n3,k_IBM_start,k_IBM_end,h)
        implicit none

        integer, intent(in) :: n3,k_IBM_start,k_IBM_end
        logical                 :: onBounds
        real*8, dimension(n3), intent(in):: f
        real*8, dimension(n3), intent(out):: d1f
        real*8, intent(in)      :: h

        integer :: k
        real(kind=8) A

        A =  1.0d0      /   (h)

        do k = k_IBM_start,k_IBM_end

            onBounds = (k.le.1).or.(k.ge.n3)
        
            if (.not.onBounds) d1f(k) = d1f(k) + A*(f(k) - f(k-1))

            if (k.eq.1) d1f(1) = d1f(1) + A*(f(1) - f(n3-1))

        enddo

        return

    end subroutine D1ssh_IBM_O2_ACC_3Dz



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!     D1s Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine D1s_IBM_O2_ACC_3Dx(f,d1f,n1,i_IBM_start,i_IBM_end,h)
        implicit none

        integer, intent(in) :: n1,i_IBM_start,i_IBM_end
        logical                 :: onBounds
        real*8, dimension(n1), intent(in):: f
        real*8, dimension(n1), intent(out):: d1f
        real*8, intent(in)      :: h

        integer :: i
        real(kind=8) A

        A =  1.0d0      /   (h)

        do i = i_IBM_start,i_IBM_end

            onBounds = (i.le.0).or.(i.ge.n1-1)
        
            if (.not.onBounds) d1f(i) = d1f(i) + A*(f(i+1) - f(i))

            if (i.eq.n1-1) d1f(n1-1) = d1f(n1-1) + A*(f(1) - f(n1-1))

        enddo

        return

    end subroutine D1s_IBM_O2_ACC_3Dx



    subroutine D1s_IBM_O2_MULTACC_3Dy(f,d1f,n2,j_IBM_start,j_IBM_end,h,g)
        implicit none

        integer, intent(in) :: n2,j_IBM_start,j_IBM_end
        logical                 :: onBounds
        real*8, dimension(n2), intent(in):: f
        real*8, dimension(n2), intent(out):: d1f
        real*8, intent(in)      :: h
        real(kind=8), dimension(:), intent(in)              :: g

        integer :: j
        real(kind=8) A

        A =  1.0d0      /   (h)

        do j = j_IBM_start,j_IBM_end

            onBounds = (j.le.0).or.(j.ge.n2-1)
        
            if (.not.onBounds) d1f(j) = d1f(j) + (A*(f(j+1) - f(j)))*g(j)

            if (j.eq.n2-1) d1f(n2-1) = d1f(n2-1) + (A*(f(n2)-f(n2-1)))*g(n2-1)

        enddo

        return

    end subroutine D1s_IBM_O2_MULTACC_3Dy



    subroutine D1s_IBM_O2_ACC_3Dz(f,d1f,n3,k_IBM_start,k_IBM_end,h)
        implicit none

        integer, intent(in) :: n3,k_IBM_start,k_IBM_end
        logical                 :: onBounds
        real*8, dimension(n3), intent(in):: f
        real*8, dimension(n3), intent(out):: d1f
        real*8, intent(in)      :: h

        integer :: k
        real(kind=8) A

        A =  1.0d0      /   (h)

        do k = k_IBM_start,k_IBM_end

            onBounds = (k.le.0).or.(k.ge.n3-1)
        
            if (.not.onBounds) d1f(k) = d1f(k) + A*(f(k+1) - f(k))

            if (k.eq.1) d1f(n3-1) = d1f(n3-1) + A*(f(1) - f(n3-1))

        enddo

        return

    end subroutine D1s_IBM_O2_ACC_3Dz

end module DRP_IBM