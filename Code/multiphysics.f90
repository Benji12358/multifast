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
        use start_settings, only:export_outflow, import_inflow

        use embedded_data, only: use_embedded
        use embedded_fringe_data, only: \
            use_fringe_embedded => use_fringe

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

        !!!!!! This has to be commented, see below
        ! use IVELOCITY_IO_results_writer, only:VELOCITY_IO_2D_write_fields
        ! use ISCALAR_IO_results_writer, only:SCALAR_IO_2D_write_fields
        ! use IVELOCITY_IO_loader, only:VELOCITY_IO_load_2D_fields
        ! use ISCALAR_IO_loader, only:SCALAR_IO_load_2D_fields

        use FRINGE_data, only  :\
            mc_q1_inflow => q1_inflow, \
            mc_q2_inflow => q2_inflow, \
            mc_q3_inflow => q3_inflow, \
            mc_sca_inflow => sca_inflow

        use following_data, only: use_following
        use following_fringe_data, only: \
            use_fringe_following => use_fringe

        use following_fringe_solver, only: \
          following_compute_fringe_force_x => compute_fringe_force_x

        use following_velocity_solver, only: \
          following_add_action_x => add_action_x, \
          following_update_velocity => update_velocity

        use following_scalar_solver, only: \
          following_solve_scalar => solve_scalar

        use following_physical_fields, only: \
            following_q3_z => q3_z, \
            following_q2_y => q2_y, \
            following_q1_x => q1_x, \
            following_q2_x => q2_x, \
            following_q3_x => q3_x

        use following_fringe_data, only: \
            following_f3_fringe_x => f3_fringe_x, \
            following_f2_fringe_x => f2_fringe_x, \
            following_f1_fringe_x => f1_fringe_x

        implicit  none

        integer, intent(in) :: ntime
        integer             :: ns, sub_cpt, n_iter, iter,test
        logical             :: triple_decomp

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

                  if (use_fringe_embedded) then
                      call embedded_compute_fringe_force_x(embedded_q1_x, embedded_q2_x, embedded_q3_x, embedded_f1_fringe_x, embedded_f2_fringe_x, embedded_f3_fringe_x) !=> Get F(n)
                      call embedded_add_action_x(embedded_f1_fringe_x, embedded_f2_fringe_x, embedded_f3_fringe_x)
                  endif

                  call embedded_update_velocity(ntime, ns)
                  if (SCA_state.eq.1) call embedded_solve_scalar(embedded_q1_x, embedded_q2_y, embedded_q3_z, ns)
                  
                  call export_velocity_profile_from_embedded
                  if (SCA_state.eq.1) call export_sca_profile_from_embedded

              endif

              !!!!!! This was commented because it multiplied the total computational time by a factor 6-7
              !!!!!! This is very slow and the communication time between procs to write a 2D slice is not efficient
              !!!!!! Therefore, all lines of code linked to 'export_outflow' has been commented in the ENTIRE code
              ! if (import_inflow) then
              !   call VELOCITY_IO_load_2D_fields(ntime, ns, mc_q1_inflow, mc_q2_inflow, mc_q3_inflow)
              !   if(SCA_state==1) call SCALAR_IO_load_2D_fields(ntime, ns, mc_sca_inflow)
              ! endif

              if (use_fringe) then
                  call compute_fringe_force_x(q1_x, q2_x, q3_x, f1_fringe_x, f2_fringe_x, f3_fringe_x,ntime,(ns.eq.nb_substep)) !=> Get F(n)
                  call add_action_x(f1_fringe_x, f2_fringe_x, f3_fringe_x)
              endif

              call update_velocity(ntime, ns)
              if(SCA_state==1) call solve_scalar(ns)

              ! Either, the embedded channel is used and so the flowrate has to be changed in this channel
              if ((use_fringe).and.(IBM_activated).and.(use_embedded)) call export_bulk_velocity

              ! Or, the flowrate has to be changed at the inlet of the main channel
              if ((use_fringe).and.(IBM_activated).and.(.not.use_embedded)) call adjust_bulk_velocity

              !!!!!! This was commented because it multiplied the total computational time by a factor 6-7
              !!!!!! This is very slow and the communication time between procs to write a 2D slice is not efficient
              !!!!!! Therefore, all lines of code linked to 'export_outflow' has been commented in the ENTIRE code
              ! if (export_outflow) then
              !   call VELOCITY_IO_2D_write_fields(ntime, ns)
              !   if(SCA_state==1) call SCALAR_IO_2D_write_fields(ntime, ns)
              ! endif

              if (use_following) then

                  call export_velocity_profile_to_following
                  if (SCA_state.eq.1) call export_sca_profile_to_following

                  if (use_fringe_following) then
                      call following_compute_fringe_force_x(following_q1_x, following_q2_x, following_q3_x, following_f1_fringe_x, following_f2_fringe_x, following_f3_fringe_x) !=> Get F(n)
                      call following_add_action_x(following_f1_fringe_x, following_f2_fringe_x, following_f3_fringe_x)
                  endif

                  call following_update_velocity(ntime, ns)
                  if (SCA_state.eq.1) call following_solve_scalar(following_q1_x, following_q2_y, following_q3_z, ns)

              endif

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

      subroutine export_velocity_profile_from_embedded()

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

      end subroutine export_velocity_profile_from_embedded

      subroutine export_sca_profile_from_embedded()

        use FRINGE_data, only: \
            mc_sca_inflow => sca_inflow

        use embedded_fringe_data, only: \
            n_interest_region, delta_activation

        use embedded_scalar_data, only: \
            embedded_sca_x => sca_x

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

      end subroutine export_sca_profile_from_embedded

      subroutine export_velocity_profile_to_following()

        use following_FRINGE_data, only: \
            following_q1_inflow => q1_inflow, \
            following_q2_inflow => q2_inflow, \
            following_q3_inflow => q3_inflow

        use following_data

        implicit none

        integer :: j, k
        integer :: i0

        i0 = n1m-1

        ! Here, the following channel is NOT laminar, that implies
        ! there is A dependence upon z

        ! Moreover, the two channels (following and main channel) 
        ! SHOULD have the same number of nodes in the y-direction and in the z-direction

        do j=xstart_f(2),xend_f(2)
            do k=xstart_f(3),xend_f(3)

                following_q1_inflow(j,k) = q1_x(i0,j,k)
                following_q2_inflow(j,k) = q2_x(i0,j,k)
                following_q3_inflow(j,k) = q3_x(i0,j,k)

            enddo
        enddo

      end subroutine export_velocity_profile_to_following

      subroutine export_sca_profile_to_following()

        use following_FRINGE_data, only: \
            following_sca_inflow => sca_inflow

        use following_data

        implicit none

        integer :: j, k
        integer :: i0

        i0 = n1m-1

        ! Here, the following channel is NOT laminar, that implies
        ! there is A dependence upon z

        ! Moreover, the two channels (embedded and main channel) 
        ! SHOULD have the same number of nodes in the y-direction

        do j=xstart_f(2),xend_f(2)
            do k=xstart_f(3),xend_f(3)

                following_sca_inflow(j,k) = sca_x(i0,j,k)

            enddo
        enddo

      end subroutine export_sca_profile_to_following

      subroutine export_bulk_velocity()

        use embedded_data, only: \
            u_bulk

        use decomp_2d, only: nrank 

        implicit none

        integer :: j,k, mpi_err
        integer :: i0

        real*8 :: ut1, coef

        i0 = n1m-1

        ! Perform the flowrate at the inlet...
        ut1=0.d0
        u_bulk=0.d0
        do k=xstart(3),min(xend(3),n3-1)
            do j=xstart(2),min(xend(2),n2-1)
                ut1=ut1+q1_x(i0,j,k)*cell_size_Y(j)*dx3
            enddo
        enddo

        call MPI_ALLREDUCE(ut1, u_bulk, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
        u_bulk=u_bulk/(L3*2.d0)

        if (nrank.eq.0) write(*,*) 'new Re_bulk ', ren*u_bulk

      end subroutine export_bulk_velocity

      subroutine adjust_bulk_velocity()

        use FRINGE_data, only: \
            n_interest_region, delta_activation, q1_inflow, square_q1_inflow, u_bulk_theo

        use decomp_2d, only: nrank 

        implicit none

        integer :: j,k, mpi_err
        integer :: i0

        real*8 :: delta, uc_inflow, coef
        logical :: pair_n2
        real*8, dimension(n2) :: inflow_profile

        real*8 :: ut1, u_bulk

        ! i0 = n_interest_region - 2*delta_activation
        i0 = n1m-1

        ! Perform the flowrate at the inlet...
        ut1=0.d0
        u_bulk=0.d0
        do k=xstart(3),min(xend(3),n3-1)
            do j=xstart(2),min(xend(2),n2-1)
                ut1=ut1+q1_x(i0,j,k)*cell_size_Y(j)*dx3
            enddo
        enddo

        call MPI_ALLREDUCE(ut1, u_bulk, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
        u_bulk=u_bulk/(L3*2.d0)

        coef = u_bulk/u_bulk_theo

        do j=xstart(2),min(xend(2),n2m)
            q1_inflow(j,:) = square_q1_inflow(j) * coef
        enddo


      end subroutine adjust_bulk_velocity

    end subroutine perform_multiphysics_MHD

end module multiphysics
