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

    use SCALAR_data, only: SCA_state, sca_x, sca_y
    use SCALAR_solver, only: Fextsca

    implicit none

contains

    ! First and naive implementation of the fringe where the inflow was imposed like the open boundaries  
    subroutine set_inflow

      implicit  none
      integer               :: i,j,k

        ! the inflow is imposed
        if (streamwise==1) then
          do j=xstart(2),min(xend(2),n2-1)
            do k=xstart(3),min(xend(3),n3-1)
              q1_x(n_interest_region_start,j,k) = q1_inflow(j,k) 
              q2_x(n_interest_region_start,j,k) = q2_inflow(j,k) 
              q3_x(n_interest_region_start,j,k) = q3_inflow(j,k)
            enddo
          enddo
        endif

        if (SCA_state==1) then
          do j=xstart(2),min(xend(2),n2-1)
            do k=xstart(3),min(xend(3),n3-1)
              sca_x(n_interest_region_start,j,k) = sca_inflow(j,k)
            enddo
          enddo

          call transpose_x_to_y(sca_x,sca_y)
        endif

    end subroutine set_inflow

    subroutine compute_fringe_force_x(q1_x, q2_x, q3_x, f1_fringe_x, f2_fringe_x, f3_fringe_x, ntime, last_substep)

      use DNS_settings, only: dt

      implicit none

      integer, intent(in)   :: ntime
      logical, intent(in)   :: last_substep
      integer :: i,j,k,i0
      integer :: shift_i
      real*8  :: u,v,w,sca
      real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in) :: q1_x, q2_x, q3_x
      real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(inout) :: f1_fringe_x, f2_fringe_x, f3_fringe_x

      real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: sca_fringe_x

      lambda_x = 0.d0
      f1_fringe_x = 0.d0
      f2_fringe_x = 0.d0
      f3_fringe_x = 0.d0
      sca_fringe_x = 0.d0

      do i=xstart(1),min(xend(1),n1-1)
        lambda_x(i) = (max_strength_damping/dt) * ( fringe_smooth_step_function( real(i-n_fringe_start,8)/real(n_delta_rise, 8) ) - fringe_smooth_step_function( real(i-n_fringe_end,8)/real(n_delta_fall,8) + 1) )
      enddo

      ! If another function is used in the fringe region ....
      ! do i=xstart(1),min(xend(1),n1-1)
      !   lambda_x(i) = max_strength_damping * fringe_linear_function( real(i-n_fringe_start,8) )
      ! enddo

      do i=n_fringe_start,n_fringe_end
        do j=xstart(2),min(xend(2),n2-1)
          do k=xstart(3),min(xend(3),n3-1)

              f1_fringe_x(i,j,k) = lambda_x(i) * ( q1_inflow(j,k) - q1_x(i,j,k) )
              f2_fringe_x(i,j,k) = lambda_x(i) * ( q2_inflow(j,k) - q2_x(i,j,k) )
              f3_fringe_x(i,j,k) = lambda_x(i) * ( q3_inflow(j,k) - q3_x(i,j,k) )

          enddo
        enddo
      enddo

      if (SCA_state==1) then

        do i=n_fringe_start,n_fringe_end
          do j=xstart(2),min(xend(2),n2-1)
            do k=xstart(3),min(xend(3),n3-1)

              sca_fringe_x(i,j,k) = lambda_x(i) * ( sca_inflow(j,k) - sca_x(i,j,k) )

            enddo
          enddo
        enddo

        call transpose_x_to_y(sca_fringe_x,Fextsca)

      endif

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

end module FRINGE_solver
