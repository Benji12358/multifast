module multiphysics
!
! MASTER MODULE
!
contains

    subroutine perform_multiphysics(ntime)
        use Bubble_generator

        integer, intent(in) :: ntime

        call move_bubble(ntime,sub_step)

    end subroutine perform_multiphysics



    subroutine perform_multiphysics_without_bubbles(ntime)

        use VELOCITY_solver
        use SCALAR_solver
        use time_schemes

        implicit  none

        integer, intent(in) :: ntime
        integer             :: ns, sub_cpt


        !call generate_bubbles(1,1)

        ! Fluid Velocity
        do ns=1,nb_substep
            ! First estimation of velocity with NO scalar effect


            call update_velocity(ntime, ns)
            call solve_scalar(ns)
        end do

    end subroutine perform_multiphysics_without_bubbles


   subroutine perform_multiphysics_MHD(ntime)

        use VELOCITY_solver
        use MHD_solver
        use FRINGE_solver
        use SCALAR_solver
        use time_schemes
        use IBM_settings, only:IBM_activated
        use MHD_data, only:MHD_state
        use FRINGE_data, only:use_fringe

        use embedded_data, only: use_embedded

        use embedded_fringe_solver, only: \
          embedded_compute_fringe_force_x => compute_fringe_force_x

        use embedded_velocity_solver, only: \
          embedded_add_action_x => add_action_x, \
          embedded_update_velocity => update_velocity

        use embedded_scalar_solver, only: \
          embedded_solve_scalar => solve_scalar

        use embedded_physical_fields, only: \
            embedded_q3_z => q3_z, \
            embedded_q2_y => q2_y, \
            embedded_q1_x => q1_x, \
            embedded_q2_x => q2_x, \
            embedded_q3_x => q3_x

        use embedded_fringe_data, only: \
            embedded_f3_fringe_x => f3_fringe_x, \
            embedded_f2_fringe_x => f2_fringe_x, \
            embedded_f1_fringe_x => f1_fringe_x

        implicit  none

        integer, intent(in) :: ntime
        integer             :: ns, sub_cpt, n_iter, iter,test

        n_iter = 1

!**************************************************

          !*********** Active MHD: Lorentz forces in NS equations ********
          !*********** Use first order Euler scheme for MHD ******
          if (MHD_state.eq.2) then

            call solve_MHD(q1_x, q2_y, q3_z, fb1_MHD_x, fb2_MHD_x, fb3_MHD_x,ntime,.true.) !=> Get F(n)
            call add_action_x(fb1_MHD_x, fb2_MHD_x, fb3_MHD_x)
          endif

          do ns=1,nb_substep

              !*********** Active MHD: Lorentz forces in NS equations ********
              !*********** Use the same time advancement scheme (RK3 or AB2) for MHD as for NS ******
              if (MHD_state.eq.1) then

                call solve_MHD(q1_x, q2_y, q3_z, fb1_MHD_x, fb2_MHD_x, fb3_MHD_x,ntime,(ns.eq.nb_substep)) !=> Get F(n)
                call add_action_x(fb1_MHD_x, fb2_MHD_x, fb3_MHD_x)
              endif

              if (use_embedded) then

                  if (use_fringe) then
                      call embedded_compute_fringe_force_x(embedded_q1_x, embedded_q2_x, embedded_q3_x, embedded_f1_fringe_x, embedded_f2_fringe_x, embedded_f3_fringe_x) !=> Get F(n)
                      call embedded_add_action_x(embedded_f1_fringe_x, embedded_f2_fringe_x, embedded_f3_fringe_x)
                  endif

                  call embedded_update_velocity(ntime, ns)
                  if (SCA_state.eq.1) call embedded_solve_scalar(embedded_q1_x, embedded_q2_y, embedded_q3_z, ns)

                  call export_velocity_profile
                  if (SCA_state.eq.1) call export_sca_profile

              endif

              if (use_fringe) then
                  call compute_fringe_force_x(q1_x, q2_x, q3_x, f1_fringe_x, f2_fringe_x, f3_fringe_x,ntime,(ns.eq.nb_substep)) !=> Get F(n)
                  call add_action_x(f1_fringe_x, f2_fringe_x, f3_fringe_x)
              endif

              call update_velocity(ntime, ns)
              if(SCA_state==1) call solve_scalar(ns)

              if (use_embedded) call export_bulk_velocity

              if ((use_fringe).and.(IBM_activated)) call adjust_bulk_velocity

         enddo

!********************Procédure d'iteration **********************
!******************* nbre d'iter n_iter gt 1 ********************
!
!       do iter = 2,n_iter
!          call solve_MHD(q1_x, q2_y, q3_z, fb1_MHD_x, fb2_MHD_x, fb3_MHD_x)  !=> F_(n+1)
!          call add_action(fb1_MHD_x, fb2_MHD_x, fb3_MHD_x)
! !
!
!          do ns=1,nb_substep
!            call update_multiphysics_velocity(ntime, ns) !=> Take into account Fn and Fn+1
!            call solve_scalar(q1_x, q2_y, q3_z, ns)
!          enddo
!
!       enddo


    !*********Passive MHD: solve MHD fields without coupling with NS (Lorentz Forces = 0)********
      if (MHD_state.eq.3) then
        call solve_MHD(q1_x, q2_y, q3_z, fb1_MHD_x, fb2_MHD_x, fb3_MHD_x,ntime,.true.) !=> Get F(n)
      endif
    !**************************************************

    contains 

      subroutine export_velocity_profile()

        use embedded_physical_fields, only: \
            embedded_q1_x => q1_x

        use FRINGE_data, only: \
            mc_q1_inflow => q1_inflow, \
            mc_q2_inflow => q2_inflow, \
            mc_q3_inflow => q3_inflow

        use embedded_fringe_data, only: \
            n_interest_region, delta_activation

        use embedded_data

        use embedded_start_settings, only: index_for_output

        implicit none

        integer :: j
        integer :: i0, k0

        i0 = index_for_output

        ! Here, because the embedded channel is laminar, that implies
        ! there is no dependence upon z
        k0 = xstart_e(3)

        ! Moreover, the two channels (embedded and main channel) 
        ! SHOULD have the same number of nodes in the y-direction

        do j=xstart_e(2),xend_e(2)

          mc_q1_inflow(j,:) = embedded_q1_x(i0,j,k0)
          mc_q2_inflow(j,:) = embedded_q2_x(i0,j,k0)
          mc_q3_inflow(j,:) = embedded_q3_x(i0,j,k0)

        enddo

      end subroutine export_velocity_profile

      subroutine export_sca_profile()

        use embedded_scalar_data, only: \
            embedded_sca_x => sca_x

        use FRINGE_data, only: \
            mc_sca_inflow => sca_inflow

        use embedded_fringe_data, only: \
            n_interest_region, delta_activation

        use embedded_data

        use embedded_start_settings, only: index_for_output

        implicit none

        integer :: j
        integer :: i0, k0

        i0 = index_for_output

        ! Here, because the embedded channel is laminar, that implies
        ! there is no dependence upon z
        k0 = xstart_e(3)

        ! Moreover, the two channels (embedded and main channel) 
        ! SHOULD have the same number of nodes in the y-direction

        do j=xstart_e(2),xend_e(2)

          mc_sca_inflow(j,:) = embedded_sca_x(i0,j,k0)

        enddo

      end subroutine export_sca_profile

      subroutine export_bulk_velocity()

        use FRINGE_data, only: \
            n_interest_region, delta_activation

        use embedded_data, only: \
            u_bulk

        use decomp_2d, only: nrank 

        implicit none

        integer :: j,k, mpi_err
        integer :: i0

        real*8 :: ut1

        i0 = n_interest_region - 2*delta_activation

        ! Perform the flowrate at the inlet...
        ut1=0.d0
        do k=xstart(3),min(xend(3),n3-1)
            do j=xstart(2),min(xend(2),n2-1)
                ut1=ut1+q1_x(i0,j,k)*(Y(j+1)-Y(j))*dx3
            enddo
        enddo

        call MPI_ALLREDUCE(ut1, u_bulk, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
        u_bulk=u_bulk/(L3*L2)

        if (nrank.eq.0) write(*,*) 'new u_bulk ', u_bulk

      end subroutine export_bulk_velocity

      subroutine adjust_bulk_velocity()

        use FRINGE_data, only: \
            n_interest_region, delta_activation, q1_inflow

        use decomp_2d, only: nrank 

        implicit none

        integer :: j,k, mpi_err
        integer :: i0

        real*8 :: delta, uc_inflow, coef
        logical :: pair_n2
        real*8, dimension(n2) :: inflow_profile

        real*8 :: ut1, u_bulk, u_bulk_theo

        ! i0 = n_interest_region - 2*delta_activation
        i0 = n1m

        ! Perform the flowrate at the inlet...
        ut1=0.d0
        do k=xstart(3),min(xend(3),n3-1)
            do j=xstart(2),min(xend(2),n2-1)
                ut1=ut1+q1_x(i0,j,k)*(Y(j+1)-Y(j))*dx3
            enddo
        enddo

        call MPI_ALLREDUCE(ut1, u_bulk, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
        u_bulk=u_bulk/(L3*L2)

        ! if (nrank.eq.0) write(*,*) 'new u_bulk ', u_bulk
        ! if (nrank.eq.0) write(*,*) 'target Re_bulk ', ren*u_bulk

        ! q1_inflow(:,:) = u_bulk

        ! USING BLASIUS BOUNDARY LAYER
        delta = Y(5) ! profile on 4 nodes

        pair_n2 = (mod(n2, 2)==0)

        do j=ystart(2),n2/2
            if (Yc(j)<delta) then
                inflow_profile(j) = 2.d0 * (Yc(j)/delta) - 2.d0 * (Yc(j)/delta)**3 + 1.d0 * (Yc(j)/delta)**4
                inflow_profile(n2-j) = 2.d0 * (Yc(j)/delta) - 2.d0 * (Yc(j)/delta)**3 + 1.d0 * (Yc(j)/delta)**4
            else
                inflow_profile(j) = 1.d0
                inflow_profile(n2-j) = 1.d0
            endif
        enddo

        if (pair_n2) then
            inflow_profile(n2/2+1) = 1.d0
        endif

        ! Perform the flowrate at the inlet...
        ut1=0.d0
        do k=xstart(3),min(xend(3),n3-1)
            do j=xstart(2),min(xend(2),n2-1)
                ut1=ut1+inflow_profile(j)*(Y(j+1)-Y(j))*dx3
            enddo
        enddo

        call MPI_ALLREDUCE(ut1, u_bulk_theo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
        u_bulk_theo=u_bulk_theo/(L3*L2)

        coef = u_bulk/u_bulk_theo

        inflow_profile(:) = inflow_profile(:) * coef

        do j=xstart(2),min(xend(2),n2m)
            q1_inflow(j,:) = inflow_profile(j)
        enddo


      end subroutine adjust_bulk_velocity

    end subroutine perform_multiphysics_MHD

end module multiphysics
