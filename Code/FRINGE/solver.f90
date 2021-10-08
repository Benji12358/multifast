module FRINGE_solver

    use decomp_2d
    use mesh
    use boundaries
    use schemes3D_interface
    use irregular_derivative_coefficients
    use DNS_settings
    use FRINGE_data
    use physical_fields
    use mpi

    implicit none

contains

    subroutine set_inflow

      implicit  none
      integer               :: i,j,k

        ! the inflow is imposed
        if (streamwise==1) then
          do j=xstart(2),min(xend(2),n2-1)
            do k=xstart(3),min(xend(3),n3-1)
              q1_x(1,j,k) = q1_inflow(j,k)
              q2_x(1,j,k) = q2_inflow(j,k)
              q3_x(1,j,k) = q3_inflow(j,k)
            enddo
          enddo
        endif

        if (streamwise==3) then
          do i=zstart(1),zend(1)
            do j=zstart(2),zend(2)
              q1_z(i,j,1) = q1_inflow(i,j)
              q2_z(i,j,1) = q2_inflow(i,j)
              q3_z(i,j,1) = q3_inflow(i,j)
            enddo
          enddo
        endif

    end subroutine set_inflow

    subroutine compute_fringe_force_x(q1_x, q2_x, q3_x, f1_fringe_x, f2_fringe_x, f3_fringe_x, ntime, last_substep)

      implicit none

      integer, intent(in)   :: ntime
      logical, intent(in)   :: last_substep
      integer :: i,j,k
      integer :: shift_i
      real*8  :: u,v,w
      real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in) :: q1_x, q2_x, q3_x
      real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(inout) :: f1_fringe_x, f2_fringe_x, f3_fringe_x

      lambda_x = 0.d0
      f1_fringe_x = 0.d0
      f2_fringe_x = 0.d0
      f3_fringe_x = 0.d0

      do i=xstart(1),min(xend(1),n1-1)
        lambda_x(i) = max_strength_damping * ( fringe_smooth_step_function( real(i-n_fringe_start,8)/real(n_delta_rise, 8) ) - fringe_smooth_step_function( real(i-n_fringe_end,8)/real(n_delta_fall,8) + 1) )
      enddo

      ! do i=xstart(1),min(xend(1),n1-1)
      !   lambda_x(i) = max_strength_damping * fringe_linear_function( real(i-n_fringe_start,8) )
      ! enddo

      do i=n_fringe_start,min(xend(1),n1-1)
        do j=xstart(2),min(xend(2),n2-1)
          do k=xstart(3),min(xend(3),n3-1)

              ! f1_fringe_x(i,j,k) = lambda_x(i) * ( q1_inflow(j,k) - q1_x(i,j,k) )
              ! f2_fringe_x(i,j,k) = lambda_x(i) * ( q2_inflow(j,k) - q2_x(i,j,k) )
              ! f3_fringe_x(i,j,k) = lambda_x(i) * ( q3_inflow(j,k) - q3_x(i,j,k) )

              shift_i = n_fringe_start - (real(i-n_fringe_start)/real(n_fringe_region))*n_fringe_start

              if (shift_i==0) shift_i=1

              u = ( q1_inflow(j,k) - q1_x(i,j,k) ) * fringe_smooth_step_function( real(i-n_fringe_start,8)/real(n_delta_rise, 8))
              v = ( q2_inflow(j,k) - q2_x(i,j,k) ) * fringe_smooth_step_function( real(i-n_fringe_start,8)/real(n_delta_rise, 8))
              w = ( q3_inflow(j,k) - q3_x(i,j,k) ) * fringe_smooth_step_function( real(i-n_fringe_start,8)/real(n_delta_rise, 8))

              f1_fringe_x(i,j,k) = lambda_x(i) * ( u )
              f2_fringe_x(i,j,k) = lambda_x(i) * ( v )
              f3_fringe_x(i,j,k) = lambda_x(i) * ( w )

          enddo
        enddo
      enddo

    contains

    function fringe_smooth_step_function(x) result(y)

        implicit  none

        real*8, intent(in) :: x
        real*8 :: y

        if (x<=0) then

          y = 0.d0

        elseif (x>=1) then

          y = 1.d0

        else ! 0 < x < 1

          y = 1.d0 / ( 1.d0 + exp(1.d0/(x-1.d0) + 1.d0/x))

        endif

    end function fringe_smooth_step_function


    function fringe_linear_function(x) result(y)

        implicit  none

        real*8, intent(in) :: x
        real*8 :: y

        if (x<=0) then

          y = 0.d0

        elseif (x>=n_fringe_end) then

          y = 1.d0

        else ! 0 < x < 1

          y = x * 1.d0/(n_fringe_region-1)

        endif

    end function fringe_linear_function

  end subroutine compute_fringe_force_x

      subroutine compute_fringe_force_z(q1_z, q2_z, q3_z, f1_fringe_z, f2_fringe_z, f3_fringe_z, ntime, last_substep)

      implicit none

      integer, intent(in)   :: ntime
      logical, intent(in)   :: last_substep
      integer :: i,j,k
      real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in) :: q1_z, q2_z, q3_z
      real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(inout) :: f1_fringe_z, f2_fringe_z, f3_fringe_z

      do k=zstart(3),min(zend(3),n3-1)
        lambda_z(k) = max_strength_damping * ( fringe_smooth_step_function( real(k-n_fringe_start,8)/real(n_delta_rise, 8) ) - fringe_smooth_step_function( real(k-n_fringe_end,8)/real(n_delta_fall,8) + 1) )
      enddo

      do i=zstart(1),min(zend(1),n1-1)
        do j=zstart(2),min(zend(2),n2-1)
          do k=zstart(3),min(zend(3),n3-1)

              f1_fringe_z(i,j,k) = lambda_z(k) * ( q1_inflow(i,j) - q1_z(i,j,k) )
              f2_fringe_z(i,j,k) = lambda_z(k) * ( q2_inflow(i,j) - q2_z(i,j,k) )
              f3_fringe_z(i,j,k) = lambda_z(k) * ( q3_inflow(i,j) - q3_z(i,j,k) )
            
          enddo
        enddo
      enddo

    contains

    function fringe_smooth_step_function(x) result(y)

        implicit  none

        real*8, intent(in) :: x
        real*8 :: y

        if (x<=0) then

          y = 0.d0

        elseif (x>=1) then

          y = 1.d0

        else ! 0 < x < 1

          y = 1 / ( 1 + exp(1/(x-1) + 1/x))

        endif

    end function fringe_smooth_step_function

  end subroutine compute_fringe_force_z

end module FRINGE_solver
