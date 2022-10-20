module embedded_fringe_solver

    use decomp_2d

    use embedded_fringe_data
    use embedded_physical_fields
    use mpi
    use embedded_mesh
    use embedded_data

    use embedded_scalar_data, only: SCA_state, sca_x, sca_y
    use embedded_scalar_solver, only: Fextsca

    implicit none

contains

    subroutine set_inflow

      implicit  none
      integer               :: i,j,k

        ! assume inflow is square, so:
        ! u_inflow = u_bulk
        ! v_inflow = 0
        ! w_inflow = 0

        ! the inflow is imposed
        q1_x(1,:,:) = u_bulk
        q2_x(1,:,:) = 0.d0
        q3_x(1,:,:) = 0.d0

        if (SCA_state.eq.1) then

          do j=xstart_e(2),min(xend_e(2),n2-1)
            do k=xstart_e(3),min(xend_e(3),n3-1)

              sca_x(1,j,k) = sca_inflow(j,k)

            enddo
          enddo

          call transpose_x_to_y(sca_x, sca_y, decomp_embedded)

        endif

    end subroutine set_inflow

    subroutine compute_fringe_force_x(q1_x, q2_x, q3_x, f1_fringe_x, f2_fringe_x, f3_fringe_x)

      use DNS_settings, only: dt

      implicit none

      integer :: i,j,k,i0
      integer :: shift_i
      real*8  :: u,v,w,sca
      real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)), intent(in) :: q1_x, q2_x, q3_x
      real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)), intent(inout) :: f1_fringe_x, f2_fringe_x, f3_fringe_x

      real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3))              :: sca_fringe_x

      lambda_x = 0.d0
      f1_fringe_x = 0.d0
      f2_fringe_x = 0.d0
      f3_fringe_x = 0.d0
      sca_fringe_x = 0.d0

      do i=xstart_e(1),min(xend_e(1),n1-1)
        lambda_x(i) = (max_strength_damping/dt) * ( fringe_smooth_step_function( real(i-n_fringe_start,8)/real(n_delta_rise, 8) ) - fringe_smooth_step_function( real(i-n_fringe_end,8)/real(n_delta_fall,8) + 1) )
      enddo

      do i=n_fringe_start,min(xend_e(1),n1-1)
        do j=xstart_e(2),min(xend_e(2),n2-1)
          do k=xstart_e(3),min(xend_e(3),n3-1)

              ! assume inflow is square, so:
              ! u_inflow = u_bulk
              ! v_inflow = 0
              ! w_inflow = 0
              f1_fringe_x(i,j,k) = lambda_x(i) * ( u_bulk - q1_x(i,j,k) )
              f2_fringe_x(i,j,k) = lambda_x(i) * ( 0.d0 - q2_x(i,j,k) )
              f3_fringe_x(i,j,k) = lambda_x(i) * ( 0.d0 - q3_x(i,j,k) )

          enddo
        enddo
      enddo

      if (SCA_state==1) then

      do i=n_fringe_start,min(xend_e(1),n1-1)
        do j=xstart_e(2),min(xend_e(2),n2-1)
          do k=xstart_e(3),min(xend_e(3),n3-1)

                sca_fringe_x(i,j,k) = lambda_x(i) * ( sca_inflow(j,k) - sca_x(i,j,k) )

            enddo
          enddo
        enddo

        call transpose_x_to_y(sca_fringe_x,Fextsca, decomp_embedded)

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

  end subroutine compute_fringe_force_x

end module embedded_fringe_solver
