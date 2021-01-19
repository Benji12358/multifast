module FRINGE_solver
!
!*********** ASSUMING THAT STREAMWISE = 1 ******************
!
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

    subroutine set_inflow()

      implicit  none

      q1_x(1,2:n2-1,:) = q1_inflow(2:n2-1,:)
      q2_x(1,2:n2-1,:) = q2_inflow(2:n2-1,:)
      q3_x(1,2:n2-1,:) = q3_inflow(2:n2-1,:)

    end subroutine set_inflow

    subroutine compute_fringe_force(q1_x, q2_x, q3_x, f1_fringe_x, f2_fringe_x, f3_fringe_x, ntime, last_substep)

      implicit none

      integer, intent(in)   :: ntime
      logical, intent(in)   :: last_substep
      integer :: i,j,k
      real*8 :: tmp1, tmp2
      real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: q1_x, q2_x, q3_x, f1_fringe_x, f2_fringe_x, f3_fringe_x
      real*8, dimension(xstart(1):xend(1)) :: lambda_x

      do i=xstart(1),xend(1)
        lambda_x(i) = max_strength_damping * ( fringe_smooth_step_function( real((i-n_fringe_start)/n_delta_rise, 8) ) - fringe_smooth_step_function( real((i-n_fringe_end)/n_delta_fall + 1, 8)) )
      enddo

      do i=xstart(1),xend(1)
        do j=ystart(2),yend(2)
          do k=zstart(3),zend(3)
            f1_fringe_x(i,j,k) = lambda_x(i) * ( q1_inflow(j,k) - q1_x(i,j,k) )
            f2_fringe_x(i,j,k) = lambda_x(i) * ( q2_inflow(j,k) - q2_x(i,j,k) )
            f3_fringe_x(i,j,k) = lambda_x(i) * ( q3_inflow(j,k) - q3_x(i,j,k) )
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

  end subroutine compute_fringe_force

end module FRINGE_solver
