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
            call solve_scalar(q1_x, q2_y, q3_z, ns)
        end do

    end subroutine perform_multiphysics_without_bubbles


   subroutine perform_multiphysics_MHD(ntime)

        use VELOCITY_solver
        use MHD_solver
        use SCALAR_solver
        use time_schemes
        use MHD_data, only:MHD_state

        implicit  none

        integer, intent(in) :: ntime
        integer             :: ns, sub_cpt, n_iter, iter,test

        n_iter = 1

!**************************************************

          !*********** Active MHD: Lorentz forces in NS equations ********
          !*********** Use first order Euler scheme for MHD ******
          if (MHD_state.eq.2) then

            call solve_MHD(q1_x, q2_y, q3_z, fb1_MHD_x, fb2_MHD_x, fb3_MHD_x,ntime,.true.) !=> Get F(n)
            call add_action(fb1_MHD_x, fb2_MHD_x, fb3_MHD_x)
          endif

          do ns=1,nb_substep

              !*********** Active MHD: Lorentz forces in NS equations ********
              !*********** Use the same time advancement scheme (RK3 or AB2) for MHD as for NS ******
              if (MHD_state.eq.1) then

                call solve_MHD(q1_x, q2_y, q3_z, fb1_MHD_x, fb2_MHD_x, fb3_MHD_x,ntime,(ns.eq.nb_substep)) !=> Get F(n)
                call add_action(fb1_MHD_x, fb2_MHD_x, fb3_MHD_x)
              endif

            call update_velocity(ntime, ns)

            if(SCA_state==1) call solve_scalar(q1_x, q2_y, q3_z, ns)

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

    end subroutine perform_multiphysics_MHD

end module multiphysics
