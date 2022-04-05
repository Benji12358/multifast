module DRP_IBM_fine

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
    
    subroutine D2c_IBM_fine_O2_3Dx(f,d2f,n1,i_IBM_start,i_IBM_end,h,f_bef,f_aft)

        implicit none

        integer, intent(in)                                     :: n1,i_IBM_start,i_IBM_end
        real*8, intent(in)                                      :: h
        real*8, dimension(i_IBM_start:i_IBM_end), intent(in)    :: f
        real*8, dimension(i_IBM_start:i_IBM_end), intent(out)   :: d2f
        real*8                                                 :: f_bef, f_aft

        integer :: i
        real*8 A(0:1)

        A(1) =  1.d0 / h**2
        A(0) = -2.d0 / h**2

        do i = i_IBM_start+1,i_IBM_end-1

            d2f(i) = A(1)*(f(i+1) + f(i-1)) + A(0)*f(i)

        enddo

        d2f(i_IBM_start) = A(1)*(f(i_IBM_start+1) + f_bef) + A(0)*f(i_IBM_start)

        d2f(i_IBM_end) = A(1)*(f_aft + f(i_IBM_end-1)) + A(0)*f(i_IBM_end)

        return

    end subroutine D2c_IBM_fine_O2_3Dx



    subroutine D2c_IBM_fine_O2_ACC_3Dx(f,d2f,n1,i_IBM_start,i_IBM_end,h,f_bef,f_aft)

        implicit none

        integer, intent(in)                                     :: n1,i_IBM_start,i_IBM_end
        real*8, intent(in)                                      :: h
        real*8, dimension(i_IBM_start:i_IBM_end), intent(in)    :: f
        real*8, dimension(i_IBM_start:i_IBM_end), intent(out)   :: d2f
        real*8                                                 :: f_bef, f_aft

        integer :: i
        real*8 A(0:1)

        A(1) =  1.d0 / h**2
        A(0) = -2.d0 / h**2

        do i = i_IBM_start+1,i_IBM_end-1

            d2f(i) = d2f(i) + A(1)*(f(i+1) + f(i-1)) + A(0)*f(i)

        enddo

        d2f(i_IBM_start) = d2f(i_IBM_start) + A(1)*(f(i_IBM_start+1) + f_bef) + A(0)*f(i_IBM_start)

        d2f(i_IBM_end) = d2f(i_IBM_end) + A(1)*(f_aft + f(i_IBM_end-1)) + A(0)*f(i_IBM_end)

        return

    end subroutine D2c_IBM_fine_O2_ACC_3Dx



    subroutine D2c_IBM_fine_O2_MULTACC_3Dy(f,d2f,n2,j_IBM_start,j_IBM_end,h,g,f_bef,f_aft)

        implicit none

        integer, intent(in)                                     :: n2,j_IBM_start,j_IBM_end
        real*8, intent(in)                                      :: h
        real*8, dimension(j_IBM_start:j_IBM_end), intent(in)    :: f
        real*8, dimension(j_IBM_start:j_IBM_end), intent(out)   :: d2f
        real*8, dimension(:), intent(in)                        :: g
        real*8                                                 :: f_bef, f_aft

        integer :: j
        real*8 A(0:1)

        A(1) =  1.d0 / h**2
        A(0) = -2.d0*(A(1))

        do j = j_IBM_start+1,j_IBM_end-1

            d2f(j) = d2f(j) + (A(1)*(f(j+1) + f(j-1)) + A(0)*f(j))*g(j)

        enddo

        d2f(j_IBM_start) = d2f(j_IBM_start) + (A(1)*(f(j_IBM_start+1) + f_bef) + A(0)*f(j_IBM_start))*g(j_IBM_start)
        
        d2f(j_IBM_end) = d2f(j_IBM_end)+( (-2.0d0*f(j_IBM_end) + (f_aft+f(j_IBM_end-1)) ) / h**2)*g(j_IBM_end)

        return

    end subroutine D2c_IBM_fine_O2_MULTACC_3Dy



    subroutine D2c_IBM_fine_O2_3Dz(f,d2f,n3,k_IBM_start,k_IBM_end,h,f_bef,f_aft)

        implicit none

        integer, intent(in)                                     :: n3,k_IBM_start,k_IBM_end
        real*8, intent(in)                                      :: h
        real*8, dimension(k_IBM_start:k_IBM_end), intent(in)    :: f
        real*8, dimension(k_IBM_start:k_IBM_end), intent(out)   :: d2f
        real*8                                                 :: f_bef, f_aft

        integer :: k
        real*8 A(0:1)

        A(1) =  1.d0 / h**2
        A(0) = -2.d0 / h**2

        do k = k_IBM_start+1,k_IBM_end-1

            d2f(k) = A(1)*(f(k+1) + f(k-1)) + A(0)*f(k)

        enddo

        d2f(k_IBM_start) = A(1)*(f(k_IBM_start+1) + f_bef) + A(0)*f(k_IBM_start)

        d2f(k_IBM_end) = A(1)*(f_aft + f(k_IBM_end-1)) + A(0)*f(k_IBM_end)

        return

    end subroutine D2c_IBM_fine_O2_3Dz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!     D1c Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine D1c_IBM_fine_O2_3Dx(f, d1f, n1,i_IBM_start,i_IBM_end,h,f_bef,f_aft)

        implicit none

        integer, intent(in)                                     :: n1,i_IBM_start,i_IBM_end
        real*8, intent(in)                                      :: h
        real*8, dimension(i_IBM_start:i_IBM_end), intent(in)    :: f
        real*8, dimension(i_IBM_start:i_IBM_end), intent(out)   :: d1f
        real*8                                                 :: f_bef, f_aft

        integer :: i
        real*8 A

        A =  0.5d0    /h

        do i = i_IBM_start+1,i_IBM_end-1

            d1f(i) = A * (f(i+1) - f(i-1))

        enddo

        d1f(i_IBM_start) = A * (f(i_IBM_start+1) - f_bef)

        d1f(i_IBM_end) = A * (f_aft - f(i_IBM_end-1))

        return

    end subroutine D1c_IBM_fine_O2_3Dx



    subroutine D1c_IBM_fine_O2_MULT_3Dy(f,d1f,n2,j_IBM_start,j_IBM_end,h,g,f_bef,f_aft)

        implicit none

        integer, intent(in)                                     :: n2,j_IBM_start,j_IBM_end
        real*8, intent(in)                                      :: h
        real*8, dimension(j_IBM_start:j_IBM_end), intent(in)    :: f
        real*8, dimension(j_IBM_start:j_IBM_end), intent(out)   :: d1f
        real*8, dimension(:), intent(in)                        :: g
        real*8                                                 :: f_bef, f_aft

        integer :: j
        real*8 A

        A =  0.5d0    /h

        do j = j_IBM_start+1,j_IBM_end-1

            d1f(j) = (A*(f(j+1) - f(j-1)))*g(j)

        enddo

        d1f(j_IBM_start) = A * (f_bef - f(j_IBM_start+1))*g(j_IBM_start)

        d1f(j_IBM_end) = A * (f_aft- f(j_IBM_end-1))*g(j_IBM_end)

        return

    end subroutine D1c_IBM_fine_O2_MULT_3Dy



    subroutine D1c_IBM_fine_O2_MULTACC_3Dy(f,d1f,n2,j_IBM_start,j_IBM_end,h,g,f_bef,f_aft)

        implicit none

        integer, intent(in)                                     :: n2,j_IBM_start,j_IBM_end
        real*8, intent(in)                                      :: h
        real*8, dimension(j_IBM_start:j_IBM_end), intent(in)    :: f
        real*8, dimension(j_IBM_start:j_IBM_end), intent(out)   :: d1f
        real*8, dimension(:), intent(in)                        :: g
        real*8                                                 :: f_bef, f_aft

        integer :: j
        real*8 A

        A =  0.5d0    /h

        do j = j_IBM_start+1,j_IBM_end-1

            d1f(j) = d1f(j) + (A*(f(j+1) - f(j-1)))*g(j)

        enddo

        d1f(j_IBM_start) = d1f(j_IBM_start) + A * (f_bef - f(j_IBM_start+1))*g(j_IBM_start)

        d1f(j_IBM_end) = d1f(j_IBM_end) + A * (f_aft - f(j_IBM_end-1))*g(j_IBM_end)

        return

    end subroutine D1c_IBM_fine_O2_MULTACC_3Dy



    subroutine D1c_IBM_fine_O2_3Dz(f, d1f, n3,k_IBM_start,k_IBM_end,h,f_bef,f_aft)

        implicit none

        integer, intent(in)                                     :: n3,k_IBM_start,k_IBM_end
        real*8, intent(in)                                      :: h
        real*8, dimension(k_IBM_start:k_IBM_end), intent(in)    :: f
        real*8, dimension(k_IBM_start:k_IBM_end), intent(out)   :: d1f
        real*8                                                 :: f_bef, f_aft

        integer :: k
        real*8 A

        A =  0.5d0    /h

        do k = k_IBM_start+1,k_IBM_end-1

            d1f(k) = A * (f(k+1) - f(k-1))

        enddo

        d1f(k_IBM_start) = A * (f(k_IBM_start+1) - f_bef)

        d1f(k_IBM_end) = A * (f_aft - f(k_IBM_end-1))

        return

    end subroutine D1c_IBM_fine_O2_3Dz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!     D0ssh Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine D0ssh_IBM_fine_O2_3Dx(f,ff,n1,i_IBM_start,i_IBM_end,f_bef)
        implicit none

        integer, intent(in)                                     :: n1,i_IBM_start,i_IBM_end
        real*8, dimension(i_IBM_start:i_IBM_end), intent(in)    :: f
        real*8, dimension(i_IBM_start:i_IBM_end), intent(out)   :: ff
        real*8                                                 :: f_bef

        integer :: i

        do i = i_IBM_start+1,i_IBM_end

            ff(i) = 0.5d0*(f(i) + f(i-1))

        enddo

        ff(i_IBM_start) = 0.5d0*(f(i_IBM_start) + f_bef)

        return

    end subroutine D0ssh_IBM_fine_O2_3Dx



    subroutine D0ssh_IBM_fine_O2_3Dy(f,ff,n2,j_IBM_start,j_IBM_end,f_bef)
        implicit none

        integer, intent(in)                                     :: n2,j_IBM_start,j_IBM_end
        real*8, dimension(j_IBM_start:j_IBM_end), intent(in)    :: f
        real*8, dimension(j_IBM_start:j_IBM_end), intent(out)   :: ff
        real*8                                                 :: f_bef

        integer :: j

        do j = j_IBM_start+1,j_IBM_end

            ff(j) = 0.5d0*(f(j) + f(j-1))

        enddo

        ff(j_IBM_start) = 0.5d0*(f(j_IBM_start) + f_bef)

        return

    end subroutine D0ssh_IBM_fine_O2_3Dy



    subroutine D0ssh_IBM_fine_O2_3Dz(f,ff,n3,k_IBM_start,k_IBM_end,f_bef)
        implicit none

        integer, intent(in)                                     :: n3,k_IBM_start,k_IBM_end
        real*8, dimension(k_IBM_start:k_IBM_end), intent(in)    :: f
        real*8, dimension(k_IBM_start:k_IBM_end), intent(out)   :: ff
        real*8                                                 :: f_bef

        integer :: k

        do k = k_IBM_start+1,k_IBM_end

            ff(k) = 0.5d0*(f(k) + f(k-1))

        enddo

        ff(k_IBM_start) = 0.5d0*(f(k_IBM_start) + f_bef)

        return

    end subroutine D0ssh_IBM_fine_O2_3Dz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!     D0s Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine D0s_IBM_fine_O2_3Dx(f,ff,n1,i_IBM_start,i_IBM_end,f_aft)
        implicit none

        integer, intent(in)                                     :: n1,i_IBM_start,i_IBM_end
        real*8, dimension(i_IBM_start:i_IBM_end), intent(in)    :: f
        real*8, dimension(i_IBM_start:i_IBM_end), intent(out)   :: ff
        real*8                                                 :: f_aft

        integer :: i

        do i = i_IBM_start,i_IBM_end-1

            ff(i) = 0.5d0*(f(i+1) + f(i))

        enddo

        ff(i_IBM_end) =  0.5d0*(f_aft + f(i_IBM_end))

        return

    end subroutine D0s_IBM_fine_O2_3Dx



    subroutine D0s_IBM_fine_O2_3Dy(f,ff,n2,j_IBM_start,j_IBM_end,f_aft)
        implicit none

        integer, intent(in)                                     :: n2,j_IBM_start,j_IBM_end
        real*8, dimension(j_IBM_start:j_IBM_end), intent(in)    :: f
        real*8, dimension(j_IBM_start:j_IBM_end), intent(out)   :: ff
        real*8                                                 :: f_aft

        integer :: j

        do j = j_IBM_start,j_IBM_end-1

            ff(j) = 0.5d0*(f(j+1) + f(j))

        enddo

        ff(j_IBM_end) = 0.5d0*(f_aft + f(j_IBM_end))

        return

    end subroutine D0s_IBM_fine_O2_3Dy




    subroutine D0s_IBM_fine_O2_3Dz(f,ff,n3,k_IBM_start,k_IBM_end,f_aft)
        implicit none

        integer, intent(in)                                     :: n3,k_IBM_start,k_IBM_end
        real*8, dimension(k_IBM_start:k_IBM_end), intent(in)    :: f
        real*8, dimension(k_IBM_start:k_IBM_end), intent(out)   :: ff
        real*8                                                 :: f_aft

        integer :: k

        do k = k_IBM_start,k_IBM_end-1

            ff(k) = 0.5d0*(f(k+1) + f(k))

        enddo

        ff(k_IBM_end) =  0.5d0*(f_aft + f(k_IBM_end))

        return

    end subroutine D0s_IBM_fine_O2_3Dz



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!     D1ssh Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine D1ssh_IBM_fine_O2_ACC_3Dx(f,d1f,n1,i_IBM_start,i_IBM_end,h,f_bef)
        implicit none

        integer, intent(in)                                     :: n1,i_IBM_start,i_IBM_end
        real*8, intent(in)                                      :: h
        real*8, dimension(i_IBM_start:i_IBM_end), intent(in)    :: f
        real*8, dimension(i_IBM_start:i_IBM_end), intent(out)   :: d1f
        real*8                                                 :: f_bef

        integer :: i
        real(kind=8) A

        A =  1.0d0      /   (h)

        do i = i_IBM_start+1,i_IBM_end

            d1f(i) = d1f(i) + A*(f(i) - f(i-1))

        enddo

        d1f(i_IBM_start) = d1f(i_IBM_start) + A*(f(i_IBM_start) - f_bef)

        return

    end subroutine D1ssh_IBM_fine_O2_ACC_3Dx



    subroutine D1ssh_IBM_fine_O2_MULTACC_3Dy(f,d1f,n2,j_IBM_start,j_IBM_end,h,g,f_bef)
        implicit none

        integer, intent(in)                                     :: n2,j_IBM_start,j_IBM_end
        real*8, intent(in)                                      :: h
        real*8, dimension(j_IBM_start:j_IBM_end), intent(in)    :: f
        real*8, dimension(j_IBM_start:j_IBM_end), intent(out)   :: d1f
        real*8, dimension(:), intent(in)                        :: g
        real*8                                                 :: f_bef

        integer :: j
        real(kind=8) A

        A =  1.0d0      /   (h)

        do j = j_IBM_start+1,j_IBM_end

            d1f(j) = d1f(j) + A*(f(j) - f(j-1))

        enddo

        d1f(j_IBM_start) = d1f(j_IBM_start) + A*(f(j_IBM_start) - f_bef)

        return

    end subroutine D1ssh_IBM_fine_O2_MULTACC_3Dy



    subroutine D1ssh_IBM_fine_O2_ACC_3Dz(f,d1f,n3,k_IBM_start,k_IBM_end,h,f_bef)
        implicit none

        integer, intent(in)                                     :: n3,k_IBM_start,k_IBM_end
        real*8, intent(in)                                      :: h
        real*8, dimension(k_IBM_start:k_IBM_end), intent(in)    :: f
        real*8, dimension(k_IBM_start:k_IBM_end), intent(out)   :: d1f
        real*8                                                  :: f_bef

        integer :: k
        real(kind=8) A

        A =  1.0d0      /   (h)

        do k = k_IBM_start+1,k_IBM_end

            d1f(k) = d1f(k) + A*(f(k) - f(k-1))

        enddo

        d1f(k_IBM_start) = d1f(k_IBM_start) + A*(f(k_IBM_start) - f_bef)

        return

    end subroutine D1ssh_IBM_fine_O2_ACC_3Dz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!     D1s Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine D1s_IBM_fine_O2_ACC_3Dx(f,d1f,n1,i_IBM_start,i_IBM_end,h,f_aft)
        implicit none

        integer, intent(in)                                     :: n1,i_IBM_start,i_IBM_end
        real*8, dimension(i_IBM_start:i_IBM_end), intent(in)    :: f
        real*8, dimension(i_IBM_start:i_IBM_end), intent(out)   :: d1f
        real*8, intent(in)                                      :: h
        real*8                                                 :: f_aft

        integer :: i
        real*8 A

        A =  1.0d0      /   (h)
        
        do i = i_IBM_start,i_IBM_end-1

            d1f(i) = d1f(i) + A*(f(i+1) - f(i))

        enddo
        
        d1f(i_IBM_end) = d1f(i_IBM_end) + A*(f_aft - f(i_IBM_end))

        return

    end subroutine D1s_IBM_fine_O2_ACC_3Dx



    subroutine D1s_IBM_fine_O2_MULTACC_3Dy(f,d1f,n2,j_IBM_start,j_IBM_end,h,g,f_aft)
        implicit none

        integer, intent(in)                                     :: n2,j_IBM_start,j_IBM_end
        real*8, dimension(j_IBM_start:j_IBM_end), intent(in)    :: f
        real*8, dimension(j_IBM_start:j_IBM_end), intent(out)   :: d1f
        real*8, intent(in)                                      :: h
        real*8, dimension(:), intent(in)                        :: g
        real*8                                                 :: f_aft

        integer :: j
        real*8 A

        A =  1.0d0      /   (h)

        do j = j_IBM_start,j_IBM_end-1

            d1f(j) = d1f(j) + (A*(f(j+1) - f(j)))*g(j)

        enddo

        d1f(j_IBM_end) = d1f(j_IBM_end) + A*(f_aft - f(j_IBM_end))*g(j_IBM_end)

        return

    end subroutine D1s_IBM_fine_O2_MULTACC_3Dy



    subroutine D1s_IBM_fine_O2_ACC_3Dz(f,d1f,n3,k_IBM_start,k_IBM_end,h,f_aft)
        implicit none

        integer, intent(in)                                     :: n3,k_IBM_start,k_IBM_end
        real*8, dimension(k_IBM_start:k_IBM_end), intent(in)    :: f
        real*8, dimension(k_IBM_start:k_IBM_end), intent(out)   :: d1f
        real*8, intent(in)                                      :: h
        real*8                                                 :: f_aft

        integer :: k
        real*8 A

        A =  1.0d0      /   (h)

        do k = k_IBM_start,k_IBM_end-1

            d1f(k) = d1f(k) + A*(f(k+1) - f(k))

        enddo
        
        d1f(k_IBM_end) = d1f(k_IBM_end) + A*(f_aft - f(k_IBM_end))

        return

    end subroutine D1s_IBM_fine_O2_ACC_3Dz


!     subroutine D2c_IBM_fine_O2_3Dx(f,d2f,h,f_x_coarse,f_y_coarse,f_z_coarse)

!         implicit none

!         real*8, intent(in)                                                              :: h
!         real*8, dimension(max(i_start_ibm_2nd_f,xstart_ibm(1)):min(i_end_ibm_2nd_f,xend_ibm(1)), &
!             max(j_start_ibm_2nd_f,xstart_ibm(2)):min(j_end_ibm_2nd_f,xend_ibm(2)), &
!             max(k_start_ibm_2nd_f,xstart_ibm(3)):min(k_end_ibm_2nd_f,xend_ibm(3))), intent(in)                                      :: f
!         real*8, dimension(max(i_start_ibm_2nd_f,xstart_ibm(1)):min(i_end_ibm_2nd_f,xend_ibm(1)), &
!             max(j_start_ibm_2nd_f,xstart_ibm(2)):min(j_end_ibm_2nd_f,xend_ibm(2)), &
!             max(k_start_ibm_2nd_f,xstart_ibm(3)):min(k_end_ibm_2nd_f,xend_ibm(3))), intent(out)                           :: d2f
!         real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))        :: f_x_coarse
!         real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))        :: f_y_coarse
!         real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))        :: f_z_coarse
!         real*8                                                                          :: f_y, f_z, f_coarse

!         integer :: i,j,k
!         integer :: i_coarse, j_coarse, k_coarse
!         real*8 A(0:1)

!         A(1) =  1.d0 / h**2
!         A(0) = -2.d0 / h**2

!         do j = max(j_start_ibm_2nd_f,xstart_ibm(2)),min(j_end_ibm_2nd_f,xend_ibm(2))
!             do k = max(k_start_ibm_2nd_f,xstart_ibm(3)),min(k_end_ibm_2nd_f,xend_ibm(3))

!                 do i = i_start_ibm_2nd_f+1,i_end_ibm_2nd_f-1

!                     d2f(i,j,k) = A(1)*(f(i+1,j,k) + f(i-1,j,k)) + A(0)*f(i,j,k)

!                 enddo

!                 ! Compute derivative using a backward method

!                 if (i_start_ibm_2nd_f.eq.1) then

!                     i_coarse = n1-1

!                 else

!                     i_coarse = (i_start_ibm_2nd_f-1)/2 + 1

!                 endif

!                 j_coarse = (j-1)/2 + 1
!                 k_coarse = (k-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )
                    
!                 endif

!                 d2f(i_start_ibm_2nd_f,j,k) = ( 2.d0*f(i_start_ibm_2nd_f+1,j,k) + f_coarse - 3.d0*f(i_start_ibm_2nd_f,j,k) ) / 6.d0*h**2


!                 ! Compute derivative using a forward method

!                 if (i_end_ibm_2nd_f.eq.(n1-1)) then

!                     i_coarse = 1

!                 else

!                     i_coarse = (i_end_ibm_2nd_f-1)/2 + 1

!                 endif

!                 j_coarse = (j-1)/2 + 1
!                 k_coarse = (k-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )
                    
!                 endif

!                 d2f(i_end_ibm_2nd_f,j,k) = ( 2.d0*f(i_end_ibm_2nd_f-1,j,k) + f_coarse - 3.d0*f(i_end_ibm_2nd_f,j,k) ) / 6.d0*h**2

!             enddo
!         enddo

!         return

!     end subroutine D2c_IBM_fine_O2_3Dx
    
!     subroutine D2c_IBM_fine_O2_ACC_3Dx(f,d2f,h,f_x_coarse,f_y_coarse,f_z_coarse)

!         implicit none

!         real*8, intent(in)                                                              :: h
!         real*8, dimension(max(i_start_ibm_2nd_f,xstart_ibm(1)):min(i_end_ibm_2nd_f,xend_ibm(1)), &
!             max(j_start_ibm_2nd_f,xstart_ibm(2)):min(j_end_ibm_2nd_f,xend_ibm(2)), &
!             max(k_start_ibm_2nd_f,xstart_ibm(3)):min(k_end_ibm_2nd_f,xend_ibm(3))), intent(in)                                      :: f
!         real*8, dimension(max(i_start_ibm_2nd_f,xstart_ibm(1)):min(i_end_ibm_2nd_f,xend_ibm(1)), &
!             max(j_start_ibm_2nd_f,xstart_ibm(2)):min(j_end_ibm_2nd_f,xend_ibm(2)), &
!             max(k_start_ibm_2nd_f,xstart_ibm(3)):min(k_end_ibm_2nd_f,xend_ibm(3))), intent(out)                           :: d2f
!         real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))        :: f_x_coarse
!         real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))        :: f_y_coarse
!         real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))        :: f_z_coarse
!         real*8                                                                          :: f_y, f_z, f_coarse

!         integer :: i,j,k
!         integer :: i_coarse, j_coarse, k_coarse
!         real*8 A(0:1)

!         A(1) =  1.d0 / h**2
!         A(0) = -2.d0 / h**2

!         do j = max(j_start_ibm_2nd_f,xstart_ibm(2)),min(j_end_ibm_2nd_f,xend_ibm(2))
!             do k = max(k_start_ibm_2nd_f,xstart_ibm(3)),min(k_end_ibm_2nd_f,xend_ibm(3))

!                 do i = i_start_ibm_2nd_f+1,i_end_ibm_2nd_f-1

!                     d2f(i,j,k) = d2f(i,j,k) + A(1)*(f(i+1,j,k) + f(i-1,j,k)) + A(0)*f(i,j,k)

!                 enddo

!                 ! Compute derivative using a backward method

!                 if (i_start_ibm_2nd_f.eq.1) then

!                     i_coarse = n1-1

!                 else

!                     i_coarse = (i_start_ibm_2nd_f-1)/2 + 1

!                 endif

!                 j_coarse = (j-1)/2 + 1
!                 k_coarse = (k-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )
                    
!                 endif

!                 d2f(i_start_ibm_2nd_f,j,k) = d2f(i_start_ibm_2nd_f,j,k) + ( 2.d0*f(i_start_ibm_2nd_f+1,j,k) + f_coarse - 3.d0*f(i_start_ibm_2nd_f,j,k) ) / 6.d0*h**2


!                 ! Compute derivative using a forward method

!                 if (i_end_ibm_2nd_f.eq.(n1-1)) then

!                     i_coarse = 1

!                 else

!                     i_coarse = (i_end_ibm_2nd_f-1)/2 + 1

!                 endif

!                 j_coarse = (j-1)/2 + 1
!                 k_coarse = (k-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )
                    
!                 endif

!                 d2f(i_end_ibm_2nd_f,j,k) = d2f(i_end_ibm_2nd_f,j,k) + ( 2.d0*f(i_end_ibm_2nd_f-1,j,k) + f_coarse - 3.d0*f(i_end_ibm_2nd_f,j,k) ) / 6.d0*h**2

!             enddo
!         enddo

!         return

!     end subroutine D2c_IBM_fine_O2_ACC_3Dx
    
!     ! subroutine D2c_IBM_fine_O2_MULTACC_3Dy(f,d2f,h,f_x_coarse,f_y_coarse,f_z_coarse)

!     !     implicit none

!     !     real*8, intent(in)                                                              :: h
!     !     real*8, dimension(max(i_start_ibm_2nd_f,xstart_ibm(1)):min(i_end_ibm_2nd_f,xend_ibm(1)), &
!     !         max(j_start_ibm_2nd_f,xstart_ibm(2)):min(j_end_ibm_2nd_f,xend_ibm(2)), &
!     !         max(k_start_ibm_2nd_f,xstart_ibm(3)):min(k_end_ibm_2nd_f,xend_ibm(3))), intent(in)                                      :: f
!     !     real*8, dimension(max(i_start_ibm_2nd_f,xstart_ibm(1)):min(i_end_ibm_2nd_f,xend_ibm(1)), &
!     !         max(j_start_ibm_2nd_f,xstart_ibm(2)):min(j_end_ibm_2nd_f,xend_ibm(2)), &
!     !         max(k_start_ibm_2nd_f,xstart_ibm(3)):min(k_end_ibm_2nd_f,xend_ibm(3))), intent(out)                           :: d2f
!     !     real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))        :: f_x_coarse
!     !     real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))        :: f_y_coarse
!     !     real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))        :: f_z_coarse
!     !     real*8                                                                          :: f_y, f_z, f_coarse

!     !     integer :: i,j,k
!     !     integer :: i_coarse, j_coarse, k_coarse
!     !     real*8 A(0:1)

!     !     A(1) =  1.d0 / h**2
!     !     A(0) = -2.d0 / h**2

!     !     do j = max(j_start_ibm_2nd_f,xstart_ibm(2)),min(j_end_ibm_2nd_f,xend_ibm(2))
!     !         do k = max(k_start_ibm_2nd_f,xstart_ibm(3)),min(k_end_ibm_2nd_f,xend_ibm(3))

!     !             do i = i_start_ibm_2nd_f+1,i_end_ibm_2nd_f-1

!     !                 d2f(i,j,k) = d2f(i,j,k) + A(1)*(f(i+1,j,k) + f(i-1,j,k)) + A(0)*f(i,j,k)

!     !             enddo

!     !             ! Compute derivative using a backward method

!     !             if (i_start_ibm_2nd_f.eq.1) then

!     !                 i_coarse = n1-1

!     !             else

!     !                 i_coarse = (i_start_ibm_2nd_f-1)/2 + 1

!     !             endif

!     !             j_coarse = (j-1)/2 + 1
!     !             k_coarse = (k-1)/2 + 1

!     !             ! case lower left
!     !             if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!     !                 ! if ... else ... but condensed
!     !                 if ((j_coarse-1).le.0) f_y = 0.d0
!     !                 if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!     !                 ! if ... else ... but condensed
!     !                 if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!     !                 if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!     !                 f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!     !             ! case upper left
!     !             elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!     !                 ! if ... else ... but condensed
!     !                 if ((j_coarse+1).ge.n2) f_y = 0.d0
!     !                 if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!     !                 ! if ... else ... but condensed
!     !                 if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!     !                 if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!     !                 f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!     !             ! case lower right
!     !             elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!     !                 ! if ... else ... but condensed
!     !                 if ((j_coarse-1).le.0) f_y = 0.d0
!     !                 if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!     !                 ! if ... else ... but condensed
!     !                 if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!     !                 if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!     !                 f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!     !             ! case upper right
!     !             elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!     !                 ! if ... else ... but condensed
!     !                 if ((j_coarse+1).ge.n2) f_y = 0.d0
!     !                 if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!     !                 ! if ... else ... but condensed
!     !                 if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!     !                 if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!     !                 f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )
                    
!     !             endif

!     !             d2f(i_start_ibm_2nd_f,j,k) = d2f(i_start_ibm_2nd_f,j,k) + ( 2.d0*f(i_start_ibm_2nd_f+1,j,k) + f_coarse - 3.d0*f(i_start_ibm_2nd_f,j,k) ) / 6.d0*h**2


!     !             ! Compute derivative using a forward method

!     !             if (i_end_ibm_2nd_f.eq.(n1-1)) then

!     !                 i_coarse = 1

!     !             else

!     !                 i_coarse = (i_start_ibm_2nd_f-1)/2 + 1

!     !             endif

!     !             j_coarse = (j-1)/2 + 1
!     !             k_coarse = (k-1)/2 + 1

!     !             ! case lower left
!     !             if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!     !                 ! if ... else ... but condensed
!     !                 if ((j_coarse-1).le.0) f_y = 0.d0
!     !                 if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!     !                 ! if ... else ... but condensed
!     !                 if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!     !                 if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!     !                 f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!     !             ! case upper left
!     !             elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!     !                 ! if ... else ... but condensed
!     !                 if ((j_coarse+1).ge.n2) f_y = 0.d0
!     !                 if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!     !                 ! if ... else ... but condensed
!     !                 if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!     !                 if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!     !                 f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!     !             ! case lower right
!     !             elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!     !                 ! if ... else ... but condensed
!     !                 if ((j_coarse-1).le.0) f_y = 0.d0
!     !                 if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!     !                 ! if ... else ... but condensed
!     !                 if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!     !                 if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!     !                 f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!     !             ! case upper right
!     !             elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!     !                 ! if ... else ... but condensed
!     !                 if ((j_coarse+1).ge.n2) f_y = 0.d0
!     !                 if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!     !                 ! if ... else ... but condensed
!     !                 if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!     !                 if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!     !                 f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )
                    
!     !             endif

!     !             d2f(i_end_ibm_2nd_f,j,k) = d2f(i_end_ibm_2nd_f,j,k) + ( 2.d0*f(i_end_ibm_2nd_f-1,j,k) + f_coarse - 3.d0*f(i_end_ibm_2nd_f,j,k) ) / 6.d0*h**2

!     !         enddo
!     !     enddo

!     !     return

!     ! end subroutine D2c_IBM_fine_O2_MULTACC_3Dy
    
!     subroutine D2c_IBM_fine_O2_3Dz(f,d2f,h,f_x_coarse,f_y_coarse,f_z_coarse)

!         implicit none

!         real*8, intent(in)                                                              :: h
!         real*8, dimension(max(i_start_ibm_2nd_f,zstart_ibm(1)):min(i_end_ibm_2nd_f,zend_ibm(1)), &
!             max(j_start_ibm_2nd_f,zstart_ibm(2)):min(j_end_ibm_2nd_f,zend_ibm(2)), &
!             max(k_start_ibm_2nd_f,zstart_ibm(3)):min(k_end_ibm_2nd_f,zend_ibm(3))), intent(in)                                      :: f
!         real*8, dimension(max(i_start_ibm_2nd_f,zstart_ibm(1)):min(i_end_ibm_2nd_f,zend_ibm(1)), &
!             max(j_start_ibm_2nd_f,zstart_ibm(2)):min(j_end_ibm_2nd_f,zend_ibm(2)), &
!             max(k_start_ibm_2nd_f,zstart_ibm(3)):min(k_end_ibm_2nd_f,zend_ibm(3))), intent(out)                           :: d2f
!         real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))        :: f_x_coarse
!         real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))        :: f_y_coarse
!         real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))        :: f_z_coarse
!         real*8                                                                          :: f_y, f_x, f_coarse

!         integer :: i,j,k
!         integer :: i_coarse, j_coarse, k_coarse
!         real*8 A(0:1)

!         A(1) =  1.d0 / h**2
!         A(0) = -2.d0 / h**2

!         do i = max(i_start_ibm_2nd_f,zstart_ibm(1)),min(i_end_ibm_2nd_f,zend_ibm(1))
!             do j = max(j_start_ibm_2nd_f,zstart_ibm(2)),min(j_end_ibm_2nd_f,zend_ibm(2))

!                 do k = k_start_ibm_2nd_f+1,k_end_ibm_2nd_f-1

!                     d2f(i,j,k) = A(1)*(f(i,j,k+1) + f(i,j,k-1)) + A(0)*f(i,j,k)

!                 enddo

!                 ! Compute derivative using a backward method

!                 if (k_start_ibm_2nd_f.eq.1) then

!                     k_coarse = n3-1

!                 else

!                     k_coarse = (k_start_ibm_2nd_f-1)/2 + 1

!                 endif

!                 i_coarse = (i-1)/2 + 1
!                 j_coarse = (j-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse-1).le.0) f_x = f_x_coarse(n1-1,j_coarse,k_coarse)
!                     if ((i_coarse-1).ge.1) f_x = f_x_coarse(i_coarse-1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse-1).le.0) f_x = f_x_coarse(n1-1,j_coarse,k_coarse)
!                     if ((i_coarse-1).ge.1) f_x = f_x_coarse(i_coarse-1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse+1).ge.n1) f_x = f_x_coarse(1,j_coarse,k_coarse)
!                     if ((i_coarse+1).le.(n1-1)) f_x = f_x_coarse(i_coarse+1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse+1).ge.n1) f_x = f_x_coarse(1,j_coarse,k_coarse)
!                     if ((i_coarse+1).le.(n1-1)) f_x = f_x_coarse(i_coarse+1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )
                    
!                 endif

!                 d2f(i,j,k_start_ibm_2nd_f) = ( 2.d0*f(i,j,k_start_ibm_2nd_f+1) + f_coarse - 3.d0*f(i,j,k_start_ibm_2nd_f) ) / 6.d0*h**2


!                 ! Compute derivative using a forward method

!                 if (k_end_ibm_2nd_f.eq.1) then

!                     k_coarse = n3-1

!                 else

!                     k_coarse = (k_end_ibm_2nd_f-1)/2 + 1

!                 endif

!                 i_coarse = (i-1)/2 + 1
!                 j_coarse = (j-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse-1).le.0) f_x = f_x_coarse(n1-1,j_coarse,k_coarse)
!                     if ((i_coarse-1).ge.1) f_x = f_x_coarse(i_coarse-1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse-1).le.0) f_x = f_x_coarse(n1-1,j_coarse,k_coarse)
!                     if ((i_coarse-1).ge.1) f_x = f_x_coarse(i_coarse-1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse+1).ge.n1) f_x = f_x_coarse(1,j_coarse,k_coarse)
!                     if ((i_coarse+1).le.(n1-1)) f_x = f_x_coarse(i_coarse+1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse+1).ge.n1) f_x = f_x_coarse(1,j_coarse,k_coarse)
!                     if ((i_coarse+1).le.(n1-1)) f_x = f_x_coarse(i_coarse+1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )
                    
!                 endif

!                 d2f(i,j,k_start_ibm_2nd_f) = ( 2.d0*f(i,j,k_start_ibm_2nd_f+1) + f_coarse - 3.d0*f(i,j,k_start_ibm_2nd_f) ) / 6.d0*h**2

!             enddo
!         enddo

!         return

!     end subroutine D2c_IBM_fine_O2_3Dz

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!     D1c Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     subroutine D1c_IBM_fine_O2_3Dx(f,d1f,h,f_x_coarse,f_y_coarse,f_z_coarse)

!         implicit none

!         real*8, intent(in)                                                              :: h
!         real*8, dimension(max(i_start_ibm_2nd_f,xstart_ibm(1)):min(i_end_ibm_2nd_f,xend_ibm(1)), &
!             max(j_start_ibm_2nd_f,xstart_ibm(2)):min(j_end_ibm_2nd_f,xend_ibm(2)), &
!             max(k_start_ibm_2nd_f,xstart_ibm(3)):min(k_end_ibm_2nd_f,xend_ibm(3))), intent(in)                                      :: f
!         real*8, dimension(max(i_start_ibm_2nd_f,xstart_ibm(1)):min(i_end_ibm_2nd_f,xend_ibm(1)), &
!             max(j_start_ibm_2nd_f,xstart_ibm(2)):min(j_end_ibm_2nd_f,xend_ibm(2)), &
!             max(k_start_ibm_2nd_f,xstart_ibm(3)):min(k_end_ibm_2nd_f,xend_ibm(3))), intent(out)                           :: d1f
!         real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))        :: f_x_coarse
!         real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))        :: f_y_coarse
!         real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))        :: f_z_coarse
!         real*8                                                                          :: f_y, f_z, f_coarse

!         integer :: i,j,k
!         integer :: i_coarse, j_coarse, k_coarse
!         real*8 A

!         A =  0.5d0/h

!         do j = max(j_start_ibm_2nd_f,xstart_ibm(2)),min(j_end_ibm_2nd_f,xend_ibm(2))
!             do k = max(k_start_ibm_2nd_f,xstart_ibm(3)),min(k_end_ibm_2nd_f,xend_ibm(3))

!                 do i = i_start_ibm_2nd_f+1,i_end_ibm_2nd_f-1

!                     d1f(i,j,k) = A*(f(i+1,j,k) + f(i-1,j,k))

!                 enddo

!                 ! Compute derivative using a backward method

!                 if (i_start_ibm_2nd_f.eq.1) then

!                     i_coarse = n1-1

!                 else

!                     i_coarse = (i_start_ibm_2nd_f-1)/2 + 1

!                 endif

!                 j_coarse = (j-1)/2 + 1
!                 k_coarse = (k-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )
                    
!                 endif

!                 d1f(i_start_ibm_2nd_f,j,k) = ( f(i_start_ibm_2nd_f+1,j,k) - f_coarse ) / 3.d0*h


!                 ! Compute derivative using a forward method

!                 if (i_end_ibm_2nd_f.eq.(n1-1)) then

!                     i_coarse = 1

!                 else

!                     i_coarse = (i_end_ibm_2nd_f-1)/2 + 1

!                 endif

!                 j_coarse = (j-1)/2 + 1
!                 k_coarse = (k-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )
                    
!                 endif

!                 d1f(i_end_ibm_2nd_f,j,k) = ( f_coarse - f(i_end_ibm_2nd_f-1,j,k) ) / 3.d0*h

!             enddo
!         enddo

!         return

!     end subroutine D1c_IBM_fine_O2_3Dx















!     subroutine D1c_IBM_fine_O2_3Dz(f,d1f,h,f_x_coarse,f_y_coarse,f_z_coarse)

!         implicit none

!         real*8, intent(in)                                                              :: h
!         real*8, dimension(max(i_start_ibm_2nd_f,zstart_ibm(1)):min(i_end_ibm_2nd_f,zend_ibm(1)), &
!             max(j_start_ibm_2nd_f,zstart_ibm(2)):min(j_end_ibm_2nd_f,zend_ibm(2)), &
!             max(k_start_ibm_2nd_f,zstart_ibm(3)):min(k_end_ibm_2nd_f,zend_ibm(3))), intent(in)                                      :: f
!         real*8, dimension(max(i_start_ibm_2nd_f,zstart_ibm(1)):min(i_end_ibm_2nd_f,zend_ibm(1)), &
!             max(j_start_ibm_2nd_f,zstart_ibm(2)):min(j_end_ibm_2nd_f,zend_ibm(2)), &
!             max(k_start_ibm_2nd_f,zstart_ibm(3)):min(k_end_ibm_2nd_f,zend_ibm(3))), intent(out)                           :: d1f
!         real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))        :: f_x_coarse
!         real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))        :: f_y_coarse
!         real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))        :: f_z_coarse
!         real*8                                                                          :: f_y, f_x, f_coarse

!         integer :: i,j,k
!         integer :: i_coarse, j_coarse, k_coarse
!         real*8 A

!         A =  0.5d0/h

!         do i = max(i_start_ibm_2nd_f,zstart_ibm(1)),min(i_end_ibm_2nd_f,zend_ibm(1))
!             do j = max(j_start_ibm_2nd_f,zstart_ibm(2)),min(j_end_ibm_2nd_f,zend_ibm(2))

!                 do k = k_start_ibm_2nd_f+1,k_end_ibm_2nd_f-1

!                     d1f(i,j,k) = A*(f(i,j,k+1) + f(i,j,k-1))

!                 enddo

!                 ! Compute derivative using a backward method

!                 if (k_start_ibm_2nd_f.eq.1) then

!                     k_coarse = n3-1

!                 else

!                     k_coarse = (k_start_ibm_2nd_f-1)/2 + 1

!                 endif

!                 i_coarse = (i-1)/2 + 1
!                 j_coarse = (j-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse-1).le.0) f_x = f_x_coarse(n1-1,j_coarse,k_coarse)
!                     if ((i_coarse-1).ge.1) f_x = f_x_coarse(i_coarse-1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse-1).le.0) f_x = f_x_coarse(n1-1,j_coarse,k_coarse)
!                     if ((i_coarse-1).ge.1) f_x = f_x_coarse(i_coarse-1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse+1).ge.n1) f_x = f_x_coarse(1,j_coarse,k_coarse)
!                     if ((i_coarse+1).le.(n1-1)) f_x = f_x_coarse(i_coarse+1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse+1).ge.n1) f_x = f_x_coarse(1,j_coarse,k_coarse)
!                     if ((i_coarse+1).le.(n1-1)) f_x = f_x_coarse(i_coarse+1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )
                    
!                 endif

!                 d1f(i,j,i_start_ibm_2nd_f) = ( f(i,j,k_start_ibm_2nd_f+1) - f_coarse ) / 3.d0*h


!                 ! Compute derivative using a forward method

!                 if (k_end_ibm_2nd_f.eq.(n3-1)) then

!                     k_coarse = 1

!                 else

!                     k_coarse = (k_end_ibm_2nd_f-1)/2 + 1

!                 endif

!                 i_coarse = (i-1)/2 + 1
!                 j_coarse = (j-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse-1).le.0) f_x = f_x_coarse(n1-1,j_coarse,k_coarse)
!                     if ((i_coarse-1).ge.1) f_x = f_x_coarse(i_coarse-1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse-1).le.0) f_x = f_x_coarse(n1-1,j_coarse,k_coarse)
!                     if ((i_coarse-1).ge.1) f_x = f_x_coarse(i_coarse-1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse+1).ge.n1) f_x = f_x_coarse(1,j_coarse,k_coarse)
!                     if ((i_coarse+1).le.(n1-1)) f_x = f_x_coarse(i_coarse+1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse+1).ge.n1) f_x = f_x_coarse(1,j_coarse,k_coarse)
!                     if ((i_coarse+1).le.(n1-1)) f_x = f_x_coarse(i_coarse+1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )
                    
!                 endif

!                 d1f(i,j,k_end_ibm_2nd_f) = ( f_coarse - f(i,j,k_end_ibm_2nd_f-1) ) / 3.d0*h

!             enddo
!         enddo

!         return

!     end subroutine D1c_IBM_fine_O2_3Dz
















! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!     D0ssh Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     subroutine D0ssh_IBM_fine_O2_3Dx(f,ff,f_x_coarse,f_y_coarse,f_z_coarse)

!         implicit none

!         real*8, dimension(max(i_start_ibm_2nd_f,xstart_ibm(1)):min(i_end_ibm_2nd_f,xend_ibm(1)), &
!             max(j_start_ibm_2nd_f,xstart_ibm(2)):min(j_end_ibm_2nd_f,xend_ibm(2)), &
!             max(k_start_ibm_2nd_f,xstart_ibm(3)):min(k_end_ibm_2nd_f,xend_ibm(3))), intent(in)                                      :: f
!         real*8, dimension(max(i_start_ibm_2nd_f,xstart_ibm(1)):min(i_end_ibm_2nd_f,xend_ibm(1)), &
!             max(j_start_ibm_2nd_f,xstart_ibm(2)):min(j_end_ibm_2nd_f,xend_ibm(2)), &
!             max(k_start_ibm_2nd_f,xstart_ibm(3)):min(k_end_ibm_2nd_f,xend_ibm(3))), intent(out)                           :: ff
!         real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))        :: f_x_coarse
!         real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))        :: f_y_coarse
!         real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))        :: f_z_coarse
!         real*8                                                                          :: f_y, f_z, f_coarse

!         integer :: i,j,k
!         integer :: i_coarse, j_coarse, k_coarse

!         do j = max(j_start_ibm_2nd_f,xstart_ibm(2)),min(j_end_ibm_2nd_f,xend_ibm(2))
!             do k = max(k_start_ibm_2nd_f,xstart_ibm(3)),min(k_end_ibm_2nd_f,xend_ibm(3))

!                 do i = i_start_ibm_2nd_f+1,i_end_ibm_2nd_f

!                     ff(i,j,k) = 0.5d0*(f(i,j,k) + f(i-1,j,k))

!                 enddo

!                 ! Compute derivative using a backward method

!                 if (i_start_ibm_2nd_f.eq.1) then

!                     i_coarse = n1-1

!                 else

!                     i_coarse = (i_start_ibm_2nd_f-1)/2 + 1

!                 endif

!                 j_coarse = (j-1)/2 + 1
!                 k_coarse = (k-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )
                    
!                 endif

!                 ff(i_start_ibm_2nd_f,j,k) = 0.5d0 * ( f(i_start_ibm_2nd_f,j,k) + f_coarse )

!             enddo
!         enddo

!         return

!     end subroutine D0ssh_IBM_fine_O2_3Dx






!     subroutine D0ssh_IBM_fine_O2_3Dz(f,ff,f_x_coarse,f_y_coarse,f_z_coarse)

!         implicit none

!         real*8, dimension(max(i_start_ibm_2nd_f,zstart_ibm(1)):min(i_end_ibm_2nd_f,zend_ibm(1)), &
!             max(j_start_ibm_2nd_f,zstart_ibm(2)):min(j_end_ibm_2nd_f,zend_ibm(2)), &
!             max(k_start_ibm_2nd_f,zstart_ibm(3)):min(k_end_ibm_2nd_f,zend_ibm(3))), intent(in)                                      :: f
!         real*8, dimension(max(i_start_ibm_2nd_f,zstart_ibm(1)):min(i_end_ibm_2nd_f,zend_ibm(1)), &
!             max(j_start_ibm_2nd_f,zstart_ibm(2)):min(j_end_ibm_2nd_f,zend_ibm(2)), &
!             max(k_start_ibm_2nd_f,zstart_ibm(3)):min(k_end_ibm_2nd_f,zend_ibm(3))), intent(out)                           :: ff
!         real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))        :: f_x_coarse
!         real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))        :: f_y_coarse
!         real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))        :: f_z_coarse
!         real*8                                                                          :: f_y, f_x, f_coarse

!         integer :: i,j,k
!         integer :: i_coarse, j_coarse, k_coarse

!         do i = max(i_start_ibm_2nd_f,zstart_ibm(1)),min(i_end_ibm_2nd_f,zend_ibm(1))
!             do j = max(j_start_ibm_2nd_f,zstart_ibm(2)),min(j_end_ibm_2nd_f,zend_ibm(2))

!                 do k = k_start_ibm_2nd_f+1,k_end_ibm_2nd_f

!                     ff(i,j,k) = 0.5d0*(f(i,j,k) + f(i,j,k-1))

!                 enddo

!                 ! Compute derivative using a backward method

!                 if (k_start_ibm_2nd_f.eq.1) then

!                     k_coarse = n3-1

!                 else

!                     k_coarse = (k_start_ibm_2nd_f-1)/2 + 1

!                 endif

!                 i_coarse = (i-1)/2 + 1
!                 j_coarse = (j-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse-1).le.0) f_x = f_x_coarse(n1-1,j_coarse,k_coarse)
!                     if ((i_coarse-1).ge.1) f_x = f_x_coarse(i_coarse-1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse-1).le.0) f_x = f_x_coarse(n1-1,j_coarse,k_coarse)
!                     if ((i_coarse-1).ge.1) f_x = f_x_coarse(i_coarse-1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse+1).ge.n1) f_x = f_x_coarse(1,j_coarse,k_coarse)
!                     if ((i_coarse+1).le.(n1-1)) f_x = f_x_coarse(i_coarse+1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse+1).ge.n1) f_x = f_x_coarse(1,j_coarse,k_coarse)
!                     if ((i_coarse+1).le.(n1-1)) f_x = f_x_coarse(i_coarse+1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )
                    
!                 endif

!                 ff(i,j,k_start_ibm_2nd_f) = 0.5d0 * ( f(i,j,k_start_ibm_2nd_f) + f_coarse )

!             enddo
!         enddo

!         return

!     end subroutine D0ssh_IBM_fine_O2_3Dz










! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!     D0s Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     subroutine D0s_IBM_fine_O2_3Dx(f,ff,f_x_coarse,f_y_coarse,f_z_coarse)

!         implicit none

!         real*8, dimension(max(i_start_ibm_2nd_f,xstart_ibm(1)):min(i_end_ibm_2nd_f,xend_ibm(1)), &
!             max(j_start_ibm_2nd_f,xstart_ibm(2)):min(j_end_ibm_2nd_f,xend_ibm(2)), &
!             max(k_start_ibm_2nd_f,xstart_ibm(3)):min(k_end_ibm_2nd_f,xend_ibm(3))), intent(in)                                      :: f
!         real*8, dimension(max(i_start_ibm_2nd_f,xstart_ibm(1)):min(i_end_ibm_2nd_f,xend_ibm(1)), &
!             max(j_start_ibm_2nd_f,xstart_ibm(2)):min(j_end_ibm_2nd_f,xend_ibm(2)), &
!             max(k_start_ibm_2nd_f,xstart_ibm(3)):min(k_end_ibm_2nd_f,xend_ibm(3))), intent(out)                           :: ff
!         real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))        :: f_x_coarse
!         real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))        :: f_y_coarse
!         real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))        :: f_z_coarse
!         real*8                                                                          :: f_y, f_z, f_coarse

!         integer :: i,j,k
!         integer :: i_coarse, j_coarse, k_coarse

!         do j = max(j_start_ibm_2nd_f,xstart_ibm(2)),min(j_end_ibm_2nd_f,xend_ibm(2))
!             do k = max(k_start_ibm_2nd_f,xstart_ibm(3)),min(k_end_ibm_2nd_f,xend_ibm(3))

!                 do i = i_start_ibm_2nd_f,i_end_ibm_2nd_f-1

!                     ff(i,j,k) = 0.5d0*(f(i+1,j,k) + f(i,j,k))

!                 enddo

!                 ! Compute derivative using a forward method

!                 if (i_end_ibm_2nd_f.eq.(n1-1)) then

!                     i_coarse = 1

!                 else

!                     i_coarse = (i_end_ibm_2nd_f-1)/2 + 1

!                 endif

!                 j_coarse = (j-1)/2 + 1
!                 k_coarse = (k-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )
                    
!                 endif

!                 ff(i_end_ibm_2nd_f,j,k) = 0.5d0*( f_coarse + f(i_end_ibm_2nd_f,j,k))

!             enddo
!         enddo

!         return

!     end subroutine D0s_IBM_fine_O2_3Dx









!     subroutine D0s_IBM_fine_O2_3Dz(f,ff,f_x_coarse,f_y_coarse,f_z_coarse)

!         implicit none

!         real*8, dimension(max(i_start_ibm_2nd_f,zstart_ibm(1)):min(i_end_ibm_2nd_f,zend_ibm(1)), &
!             max(j_start_ibm_2nd_f,zstart_ibm(2)):min(j_end_ibm_2nd_f,zend_ibm(2)), &
!             max(k_start_ibm_2nd_f,zstart_ibm(3)):min(k_end_ibm_2nd_f,zend_ibm(3))), intent(in)                                      :: f
!         real*8, dimension(max(i_start_ibm_2nd_f,zstart_ibm(1)):min(i_end_ibm_2nd_f,zend_ibm(1)), &
!             max(j_start_ibm_2nd_f,zstart_ibm(2)):min(j_end_ibm_2nd_f,zend_ibm(2)), &
!             max(k_start_ibm_2nd_f,zstart_ibm(3)):min(k_end_ibm_2nd_f,zend_ibm(3))), intent(out)                           :: ff
!         real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))        :: f_x_coarse
!         real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))        :: f_y_coarse
!         real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))        :: f_z_coarse
!         real*8                                                                          :: f_y, f_x, f_coarse

!         integer :: i,j,k
!         integer :: i_coarse, j_coarse, k_coarse

!         do i = max(i_start_ibm_2nd_f,zstart_ibm(1)),min(i_end_ibm_2nd_f,zend_ibm(1))
!             do j = max(j_start_ibm_2nd_f,zstart_ibm(2)),min(j_end_ibm_2nd_f,zend_ibm(2))

!                 do k = k_start_ibm_2nd_f,k_end_ibm_2nd_f-1

!                     ff(i,j,k) = 0.5d0*(f(i,j,k+1) + f(i,j,k))

!                 enddo

!                 ! Compute derivative using a forward method

!                 if (k_end_ibm_2nd_f.eq.(n3-1)) then

!                     k_coarse = 1

!                 else

!                     k_coarse = (k_end_ibm_2nd_f-1)/2 + 1

!                 endif

!                 i_coarse = (i-1)/2 + 1
!                 j_coarse = (j-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse-1).le.0) f_x = f_x_coarse(n1-1,j_coarse,k_coarse)
!                     if ((i_coarse-1).ge.1) f_x = f_x_coarse(i_coarse-1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse-1).le.0) f_x = f_x_coarse(n1-1,j_coarse,k_coarse)
!                     if ((i_coarse-1).ge.1) f_x = f_x_coarse(i_coarse-1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse+1).ge.n1) f_x = f_x_coarse(1,j_coarse,k_coarse)
!                     if ((i_coarse+1).le.(n1-1)) f_x = f_x_coarse(i_coarse+1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse+1).ge.n1) f_x = f_x_coarse(1,j_coarse,k_coarse)
!                     if ((i_coarse+1).le.(n1-1)) f_x = f_x_coarse(i_coarse+1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )
                    
!                 endif

!                 ff(i,j,k_end_ibm_2nd_f) = 0.5d0*( f_coarse + f(i,j,k_end_ibm_2nd_f) )

!             enddo
!         enddo

!         return

!     end subroutine D0s_IBM_fine_O2_3Dz






! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!     D1ssh Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     subroutine D1ssh_IBM_fine_O2_ACC_3Dx(f,d1f,h,f_x_coarse,f_y_coarse,f_z_coarse)

!         implicit none

!         real*8, intent(in)                                                              :: h
!         real*8, dimension(max(i_start_ibm_2nd_f,xstart_ibm(1)):min(i_end_ibm_2nd_f,xend_ibm(1)), &
!             max(j_start_ibm_2nd_f,xstart_ibm(2)):min(j_end_ibm_2nd_f,xend_ibm(2)), &
!             max(k_start_ibm_2nd_f,xstart_ibm(3)):min(k_end_ibm_2nd_f,xend_ibm(3))), intent(in)                                      :: f
!         real*8, dimension(max(i_start_ibm_2nd_f,xstart_ibm(1)):min(i_end_ibm_2nd_f,xend_ibm(1)), &
!             max(j_start_ibm_2nd_f,xstart_ibm(2)):min(j_end_ibm_2nd_f,xend_ibm(2)), &
!             max(k_start_ibm_2nd_f,xstart_ibm(3)):min(k_end_ibm_2nd_f,xend_ibm(3))), intent(out)                           :: d1f
!         real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))        :: f_x_coarse
!         real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))        :: f_y_coarse
!         real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))        :: f_z_coarse
!         real*8                                                                          :: f_y, f_z, f_coarse

!         integer :: i,j,k
!         integer :: i_coarse, j_coarse, k_coarse
!         real(kind=8) A

!         A =  1.0d0      /   (h)

!         do j = max(j_start_ibm_2nd_f,xstart_ibm(2)),min(j_end_ibm_2nd_f,xend_ibm(2))
!             do k = max(k_start_ibm_2nd_f,xstart_ibm(3)),min(k_end_ibm_2nd_f,xend_ibm(3))

!                 do i = i_start_ibm_2nd_f+1,i_end_ibm_2nd_f

!                     d1f(i,j,k) = d1f(i,j,k) + A*(f(i,j,k) - f(i-1,j,k))

!                 enddo

!                 ! Compute derivative using a backward method

!                 if (i_start_ibm_2nd_f.eq.1) then

!                     i_coarse = n1-1

!                 else

!                     i_coarse = (i_start_ibm_2nd_f-1)/2 + 1

!                 endif

!                 j_coarse = (j-1)/2 + 1
!                 k_coarse = (k-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )
                    
!                 endif

!                 d1f(i_start_ibm_2nd_f,j,k) = d1f(i_start_ibm_2nd_f,j,k) + A*(f(i_start_ibm_2nd_f,j,k) - f_coarse)

!             enddo
!         enddo

!         return

!     end subroutine D1ssh_IBM_fine_O2_ACC_3Dx








!     subroutine D1ssh_IBM_fine_O2_ACC_3Dz(f,d1f,h,f_x_coarse,f_y_coarse,f_z_coarse)

!         implicit none

!         real*8, intent(in)                                                              :: h
!         real*8, dimension(max(i_start_ibm_2nd_f,zstart_ibm(1)):min(i_end_ibm_2nd_f,zend_ibm(1)), &
!             max(j_start_ibm_2nd_f,zstart_ibm(2)):min(j_end_ibm_2nd_f,zend_ibm(2)), &
!             max(k_start_ibm_2nd_f,zstart_ibm(3)):min(k_end_ibm_2nd_f,zend_ibm(3))), intent(in)                                      :: f
!         real*8, dimension(max(i_start_ibm_2nd_f,zstart_ibm(1)):min(i_end_ibm_2nd_f,zend_ibm(1)), &
!             max(j_start_ibm_2nd_f,zstart_ibm(2)):min(j_end_ibm_2nd_f,zend_ibm(2)), &
!             max(k_start_ibm_2nd_f,zstart_ibm(3)):min(k_end_ibm_2nd_f,zend_ibm(3))), intent(out)                           :: d1f
!         real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))        :: f_x_coarse
!         real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))        :: f_y_coarse
!         real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))        :: f_z_coarse
!         real*8                                                                          :: f_y, f_x, f_coarse

!         integer :: i,j,k
!         integer :: i_coarse, j_coarse, k_coarse
!         real(kind=8) A

!         A =  1.0d0      /   (h)

!         do i = max(i_start_ibm_2nd_f,zstart_ibm(1)),min(i_end_ibm_2nd_f,zend_ibm(1))
!             do j = max(j_start_ibm_2nd_f,zstart_ibm(2)),min(j_end_ibm_2nd_f,zend_ibm(2))

!                 do k = k_start_ibm_2nd_f+1,k_end_ibm_2nd_f

!                     d1f(i,j,k) = d1f(i,j,k) + A*(f(i,j,k) - f(i,j,k-1))

!                 enddo

!                 ! Compute derivative using a backward method

!                 if (k_start_ibm_2nd_f.eq.1) then

!                     k_coarse = n3-1

!                 else

!                     k_coarse = (k_start_ibm_2nd_f-1)/2 + 1

!                 endif

!                 i_coarse = (i-1)/2 + 1
!                 j_coarse = (j-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse-1).le.0) f_x = f_x_coarse(n1-1,j_coarse,k_coarse)
!                     if ((i_coarse-1).ge.1) f_x = f_x_coarse(i_coarse-1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse-1).le.0) f_x = f_x_coarse(n1-1,j_coarse,k_coarse)
!                     if ((i_coarse-1).ge.1) f_x = f_x_coarse(i_coarse-1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse+1).ge.n1) f_x = f_x_coarse(1,j_coarse,k_coarse)
!                     if ((i_coarse+1).le.(n1-1)) f_x = f_x_coarse(i_coarse+1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse+1).ge.n1) f_x = f_x_coarse(1,j_coarse,k_coarse)
!                     if ((i_coarse+1).le.(n1-1)) f_x = f_x_coarse(i_coarse+1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )
                    
!                 endif

!                 d1f(i,j,k_start_ibm_2nd_f) = d1f(i,j,k_start_ibm_2nd_f) + A*(f(i,j,k_start_ibm_2nd_f) - f_coarse)

!             enddo
!         enddo

!         return

!     end subroutine D1ssh_IBM_fine_O2_ACC_3Dz





! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!     D1s Scheme    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     subroutine D1s_IBM_fine_O2_ACC_3Dx(f,d1f,h,f_x_coarse,f_y_coarse,f_z_coarse)

!         implicit none

!         real*8, intent(in)                                                              :: h
!         real*8, dimension(max(i_start_ibm_2nd_f,xstart_ibm(1)):min(i_end_ibm_2nd_f,xend_ibm(1)), &
!             max(j_start_ibm_2nd_f,xstart_ibm(2)):min(j_end_ibm_2nd_f,xend_ibm(2)), &
!             max(k_start_ibm_2nd_f,xstart_ibm(3)):min(k_end_ibm_2nd_f,xend_ibm(3))), intent(in)                                      :: f
!         real*8, dimension(max(i_start_ibm_2nd_f,xstart_ibm(1)):min(i_end_ibm_2nd_f,xend_ibm(1)), &
!             max(j_start_ibm_2nd_f,xstart_ibm(2)):min(j_end_ibm_2nd_f,xend_ibm(2)), &
!             max(k_start_ibm_2nd_f,xstart_ibm(3)):min(k_end_ibm_2nd_f,xend_ibm(3))), intent(out)                           :: d1f
!         real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))        :: f_x_coarse
!         real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))        :: f_y_coarse
!         real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))        :: f_z_coarse
!         real*8                                                                          :: f_y, f_z, f_coarse

!         integer :: i,j,k
!         integer :: i_coarse, j_coarse, k_coarse
!         real*8 A

!         A =  1.0d0      /   (h)

!         do j = max(j_start_ibm_2nd_f,xstart_ibm(2)),min(j_end_ibm_2nd_f,xend_ibm(2))
!             do k = max(k_start_ibm_2nd_f,xstart_ibm(3)),min(k_end_ibm_2nd_f,xend_ibm(3))

!                 do i = i_start_ibm_2nd_f,i_end_ibm_2nd_f-1

!                     d1f(i,j,k) = d1f(i,j,k) + A*(f(i+1,j,k) - f(i,j,k))

!                 enddo

!                 ! Compute derivative using a forward method

!                 if (i_end_ibm_2nd_f.eq.(n1-1)) then

!                     i_coarse = 1

!                 else

!                     i_coarse = (i_end_ibm_2nd_f-1)/2 + 1

!                 endif

!                 j_coarse = (j-1)/2 + 1
!                 k_coarse = (k-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse-1).le.0) f_z = f_z_coarse(i_coarse,j_coarse,n3-1)
!                     if ((k_coarse-1).ge.1) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse-1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((k_coarse+1).ge.n3) f_z = f_z_coarse(i_coarse,j_coarse,1)
!                     if ((k_coarse+1).le.(n3-1)) f_z = f_z_coarse(i_coarse,j_coarse,k_coarse+1)

!                     f_coarse = 0.5d0 * f_x_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_z )
                    
!                 endif

!                 d1f(i_start_ibm_2nd_f,j,k) = d1f(i_start_ibm_2nd_f,j,k) + A*(f(i_start_ibm_2nd_f,j,k) - f_coarse)

!             enddo
!         enddo

!         return

!     end subroutine D1s_IBM_fine_O2_ACC_3Dx




!     subroutine D1s_IBM_fine_O2_ACC_3Dz(f,d1f,h,f_x_coarse,f_y_coarse,f_z_coarse)

!         implicit none

!         real*8, intent(in)                                                              :: h
!         real*8, dimension(max(i_start_ibm_2nd_f,zstart_ibm(1)):min(i_end_ibm_2nd_f,zend_ibm(1)), &
!             max(j_start_ibm_2nd_f,zstart_ibm(2)):min(j_end_ibm_2nd_f,zend_ibm(2)), &
!             max(k_start_ibm_2nd_f,zstart_ibm(3)):min(k_end_ibm_2nd_f,zend_ibm(3))), intent(in)                                      :: f
!         real*8, dimension(max(i_start_ibm_2nd_f,zstart_ibm(1)):min(i_end_ibm_2nd_f,zend_ibm(1)), &
!             max(j_start_ibm_2nd_f,zstart_ibm(2)):min(j_end_ibm_2nd_f,zend_ibm(2)), &
!             max(k_start_ibm_2nd_f,zstart_ibm(3)):min(k_end_ibm_2nd_f,zend_ibm(3))), intent(out)                           :: d1f
!         real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))        :: f_x_coarse
!         real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))        :: f_y_coarse
!         real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))        :: f_z_coarse
!         real*8                                                                          :: f_y, f_x, f_coarse

!         integer :: i,j,k
!         integer :: i_coarse, j_coarse, k_coarse
!         real*8 A

!         A =  1.0d0      /   (h)

!         do i = max(i_start_ibm_2nd_f,zstart_ibm(1)),min(i_end_ibm_2nd_f,zend_ibm(1))
!             do j = max(j_start_ibm_2nd_f,zstart_ibm(2)),min(j_end_ibm_2nd_f,zend_ibm(2))

!                 do k = k_start_ibm_2nd_f,k_end_ibm_2nd_f-1

!                     d1f(i,j,k) = d1f(i,j,k) + A*(f(i,j,k+1) - f(i,j,k))

!                 enddo

!                 ! Compute derivative using a forward method

!                 if (k_end_ibm_2nd_f.eq.(n3-1)) then

!                     k_coarse = 1

!                 else

!                     k_coarse = (k_end_ibm_2nd_f-1)/2 + 1

!                 endif

!                 i_coarse = (i-1)/2 + 1
!                 j_coarse = (j-1)/2 + 1

!                 ! case lower left
!                 if (((2*j_coarse-1).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse-1).le.0) f_x = f_x_coarse(n1-1,j_coarse,k_coarse)
!                     if ((i_coarse-1).ge.1) f_x = f_x_coarse(i_coarse-1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case upper left
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse-1).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse-1).le.0) f_x = f_x_coarse(n1-1,j_coarse,k_coarse)
!                     if ((i_coarse-1).ge.1) f_x = f_x_coarse(i_coarse-1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case lower right
!                 elseif (((2*j_coarse-1).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse+1).ge.n1) f_x = f_x_coarse(1,j_coarse,k_coarse)
!                     if ((i_coarse+1).le.(n1-1)) f_x = f_x_coarse(i_coarse+1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse-1).le.0) f_y = 0.d0
!                     if ((j_coarse-1).ge.1) f_y = f_y_coarse(i_coarse,j_coarse-1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )

!                 ! case upper right
!                 elseif (((2*j_coarse).eq.j).and.((2*k_coarse).eq.k)) then

!                     ! if ... else ... but condensed
!                     if ((i_coarse+1).ge.n1) f_x = f_x_coarse(1,j_coarse,k_coarse)
!                     if ((i_coarse+1).le.(n1-1)) f_x = f_x_coarse(i_coarse+1,j_coarse,k_coarse)

!                     ! if ... else ... but condensed
!                     if ((j_coarse+1).ge.n2) f_y = 0.d0
!                     if ((j_coarse+1).le.(n2-1)) f_y = f_y_coarse(i_coarse,j_coarse+1,k_coarse)

!                     f_coarse = 0.5d0 * f_z_coarse(i_coarse,j_coarse,k_coarse) + 0.25d0 * ( f_y + f_x )
                    
!                 endif

!                 d1f(i,j,k_start_ibm_2nd_f) = d1f(i,j,k_start_ibm_2nd_f) + A*(f(i,j,k_start_ibm_2nd_f) - f_coarse)

!             enddo
!         enddo

!         return

!     end subroutine D1s_IBM_fine_O2_ACC_3Dz




    
end module DRP_IBM_fine