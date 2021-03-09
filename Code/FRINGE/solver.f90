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

    subroutine set_inflow(ntime)

      implicit  none
      integer, intent(in)   :: ntime
      integer               :: i,j

      if (ntime.gt.number_it_periodic_activation) then
        ! the outflow is reinject at the inflow
        write(*,*) 'The outflow is reinjected at the inflow'

        if (streamwise==1) then
          q1_x(1,:,:) = q1_x(n1-1,:,:)
          q2_x(1,:,:) = q2_x(n1-1,:,:)
          q3_x(1,:,:) = q3_x(n1-1,:,:)
        endif

        if (streamwise==3) then
          q1_z(:,:,1) = q1_z(:,:,n3-1)
          q2_z(:,:,1) = q2_z(:,:,n3-1)
          q3_z(:,:,1) = q3_z(:,:,n3-1)
        endif

      else
        ! the inflow is imposed
        if (streamwise==1) then
          q1_x(1,:,:) = q1_inflow(:,:)
          q2_x(1,:,:) = q2_inflow(:,:)
          q3_x(1,:,:) = q3_inflow(:,:)
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

      endif

    end subroutine set_inflow

    subroutine compute_fringe_force_x(q1_x, q2_x, q3_x, f1_fringe_x, f2_fringe_x, f3_fringe_x, ntime, last_substep)

      implicit none

      integer, intent(in)   :: ntime
      logical, intent(in)   :: last_substep
      integer :: i,j,k
      real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in) :: q1_x, q2_x, q3_x
      real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(inout) :: f1_fringe_x, f2_fringe_x, f3_fringe_x

      do i=xstart(1),xend(1)
        lambda_x(i) = max_strength_damping * ( fringe_smooth_step_function( real(i-n_fringe_start,8)/real(n_delta_rise, 8) ) - fringe_smooth_step_function( real(i-n_fringe_end,8)/real(n_delta_fall,8) + 1) )
      enddo

      do i=xstart(1),xend(1)
        do j=xstart(2),xend(2)
          do k=xstart(3),xend(3)

            if (ntime.gt.number_it_periodic_activation) then
              f1_fringe_x(i,j,k) = lambda_x(i) * ( q1_x(1,j,k) - q1_x(i,j,k) )
              f2_fringe_x(i,j,k) = lambda_x(i) * ( q2_x(1,j,k) - q2_x(i,j,k) )
              f3_fringe_x(i,j,k) = lambda_x(i) * ( q3_x(1,j,k) - q3_x(i,j,k) )
            else
              f1_fringe_x(i,j,k) = lambda_x(i) * ( q1_inflow(j,k) - q1_x(i,j,k) )
              f2_fringe_x(i,j,k) = lambda_x(i) * ( q2_inflow(j,k) - q2_x(i,j,k) )
              f3_fringe_x(i,j,k) = lambda_x(i) * ( q3_inflow(j,k) - q3_x(i,j,k) )
            endif

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

  end subroutine compute_fringe_force_x

      subroutine compute_fringe_force_z(q1_z, q2_z, q3_z, f1_fringe_z, f2_fringe_z, f3_fringe_z, ntime, last_substep)

      implicit none

      integer, intent(in)   :: ntime
      logical, intent(in)   :: last_substep
      integer :: i,j,k
      real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in) :: q1_z, q2_z, q3_z
      real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(inout) :: f1_fringe_z, f2_fringe_z, f3_fringe_z

      do k=zstart(3),zend(3)
        lambda_z(k) = max_strength_damping * ( fringe_smooth_step_function( real(k-n_fringe_start,8)/real(n_delta_rise, 8) ) - fringe_smooth_step_function( real(k-n_fringe_end,8)/real(n_delta_fall,8) + 1) )
      enddo

      do i=zstart(1),zend(1)
        do j=zstart(2),zend(2)
          do k=zstart(3),zend(3)

            if (ntime.gt.number_it_periodic_activation) then
              f1_fringe_z(i,j,k) = lambda_z(k) * ( q1_z(i,j,1) - q1_z(i,j,k) )
              f2_fringe_z(i,j,k) = lambda_z(k) * ( q2_z(i,j,1) - q2_z(i,j,k) )
              f3_fringe_z(i,j,k) = lambda_z(k) * ( q3_z(i,j,1) - q3_z(i,j,k) )
            else
              f1_fringe_z(i,j,k) = lambda_z(k) * ( q1_inflow(i,j) - q1_z(i,j,k) )
              f2_fringe_z(i,j,k) = lambda_z(k) * ( q2_inflow(i,j) - q2_z(i,j,k) )
              f3_fringe_z(i,j,k) = lambda_z(k) * ( q3_inflow(i,j) - q3_z(i,j,k) )
            endif

            
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
