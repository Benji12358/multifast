module VELOCITY_bc_controller
    implicit none

contains

    subroutine set_numerical_boundaries()

        use boundaries
        implicit none

        ! ATTENTION
        select case (BC3)

            case (NOSLIP)
                NS_DEF_BC3=Dirichlet
                NS_PR_BC3=Dirichlet
                POISSON_VEL_BC3=antisymetric
                POISSON_PR_BC3=symetric

                NS_P13_BC3=Dirichlet
                NS_P23_BC3=Dirichlet
                NS_P33_BC3=Dirichlet

                NS_Q1_BC3=Dirichlet
                NS_Q2_BC3=Dirichlet
                NS_Q3_BC3=Dirichlet

            case (FREESLIP)

                NS_Q1_BC3=symetric      ! dQ1/dx3=0.d0
                NS_Q2_BC3=symetric      ! dQ2/dx3=0.d0
                NS_Q3_BC3=antisymetric  ! Q3=0.d0   (no penetration)

                NS_P13_BC3=antisymetric ! S*Ŝ=Ŝ
                NS_P23_BC3=antisymetric ! S*Ŝ=Ŝ
                NS_P33_BC3=symetric     ! Ŝ*Ŝ=S


            case (UNBOUNDED)

                NS_DEF_BC3=periodic
                NS_PR_BC3=periodic
                POISSON_VEL_BC3=periodic
                POISSON_PR_BC3=periodic

                NS_P13_BC3=periodic
                NS_P23_BC3=periodic
                NS_P33_BC3=periodic

                NS_Q1_BC3=periodic
                NS_Q2_BC3=periodic
                NS_Q3_BC3=periodic

            case (FRINGE)

                NS_DEF_BC3=periodic
                NS_PR_BC3=periodic
                POISSON_VEL_BC3=periodic
                POISSON_PR_BC3=periodic

                NS_P13_BC3=periodic
                NS_P23_BC3=periodic
                NS_P33_BC3=periodic

                NS_Q1_BC3=periodic
                NS_Q2_BC3=periodic
                NS_Q3_BC3=periodic

        end select

        select case (BC2)
            case (NOSLIP)
                NS_DEF_BC2=Dirichlet
                NS_PR_BC2=Dirichlet
                POISSON_VEL_BC2=antisymetric
                POISSON_PR_BC2=symetric

                NS_P12_BC2=Dirichlet
                NS_P22_BC2=Dirichlet
                NS_P23_BC2=Dirichlet

                NS_Q1_BC2=Dirichlet
                NS_Q2_BC2=Dirichlet
                NS_Q3_BC2=Dirichlet

            case (FREESLIP)

                NS_Q1_BC2=symetric      ! dQ1/dx2=0.d0
                NS_Q2_BC2=antisymetric  ! Q2=0.d0
                NS_Q3_BC2=symetric      ! dQ3/dx2=0.d0

                NS_P12_BC2=antisymetric ! S*Ŝ=Ŝ
                NS_P22_BC2=symetric     ! Ŝ*Ŝ=S
                NS_P23_BC2=antisymetric ! S*Ŝ=Ŝ

            case (UNBOUNDED)

                NS_DEF_BC2=periodic
                NS_PR_BC2=periodic
                POISSON_VEL_BC2=periodic

                NS_P12_BC2=periodic
                NS_P22_BC2=periodic
                NS_P23_BC2=periodic

                NS_Q1_BC2=periodic
                NS_Q2_BC2=periodic
                NS_Q3_BC2=periodic

        end select

        select case (BC1)
            case (NOSLIP)
                NS_DEF_BC1=Dirichlet
                NS_PR_BC1=Dirichlet
                POISSON_VEL_BC1=antisymetric
                POISSON_PR_BC1=symetric

                NS_P11_BC1=Dirichlet
                NS_P12_BC1=Dirichlet
                NS_P13_BC1=Dirichlet

                NS_Q1_BC1=Dirichlet
                NS_Q2_BC1=Dirichlet
                NS_Q3_BC1=Dirichlet

            case (FREESLIP)

                NS_Q1_BC1=antisymetric  ! Q1=0.d0
                NS_Q2_BC1=symetric      ! dQ2/dx1=0.d0
                NS_Q3_BC1=symetric      ! dQ3/dx1=0.d0

                NS_P11_BC1=symetric     ! Ŝ*Ŝ=S
                NS_P12_BC1=antisymetric ! Ŝ*S=Ŝ
                NS_P13_BC1=antisymetric ! Ŝ*S=Ŝ

            case (UNBOUNDED)

                NS_DEF_BC1=periodic
                NS_PR_BC1=periodic
                POISSON_VEL_BC1=periodic
                POISSON_PR_BC1=periodic

                NS_P11_BC1=periodic
                NS_P12_BC1=periodic
                NS_P13_BC1=periodic

                NS_Q1_BC1=periodic
                NS_Q2_BC1=periodic
                NS_Q3_BC1=periodic

            case (PSEUDO_PERIODIC)

                NS_DEF_BC1=Dirichlet
                NS_PR_BC1=Dirichlet
                POISSON_VEL_BC1=antisymetric
                POISSON_PR_BC1=symetric

                NS_P11_BC1=periodic
                NS_P12_BC1=periodic
                NS_P13_BC1=periodic

                NS_Q1_BC1=Dirichlet
                NS_Q2_BC1=Dirichlet
                NS_Q3_BC1=Dirichlet

            case (OPEN)

                NS_DEF_BC1=Dirichlet
                NS_PR_BC1=Dirichlet
                POISSON_VEL_BC1=antisymetric
                POISSON_PR_BC1=symetric

                NS_P11_BC1=Dirichlet
                NS_P12_BC1=Dirichlet
                NS_P13_BC1=Dirichlet

                NS_Q1_BC1=Dirichlet
                NS_Q2_BC1=Dirichlet
                NS_Q3_BC1=Dirichlet

            case (FRINGE)

                NS_DEF_BC1=periodic
                NS_PR_BC1=periodic
                POISSON_VEL_BC1=periodic
                POISSON_PR_BC1=periodic

                NS_P11_BC1=periodic
                NS_P12_BC1=periodic
                NS_P13_BC1=periodic

                NS_Q1_BC1=periodic
                NS_Q2_BC1=periodic
                NS_Q3_BC1=periodic

        end select

        !write(*,*)'TRANSPORT_Q2S_BC2', TRANSPORT_Q2S_BC2, Dirichlet

    end subroutine set_numerical_boundaries

    subroutine apply_BC3

        use mesh
        use physical_fields
        use boundaries
        use decomp_2d

        implicit none

        ! ATTENTION
        if (BC3==NOSLIP) then

            q3_wall30(zstart(1):zend(1), zstart(2):zend(2))=0.d0
            q2_wall30(zstart(1):zend(1), zstart(2):zend(2))=0.d0
            q1_wall30(zstart(1):zend(1), zstart(2):zend(2))=0.d0
            q3_z(zstart(1):zend(1), zstart(2):zend(2), 1) =0.d0

            q3_wall31(zstart(1):zend(1), zstart(2):zend(2))=0.d0
            q2_wall31(zstart(1):zend(1), zstart(2):zend(2))=0.d0
            q1_wall31(zstart(1):zend(1), zstart(2):zend(2))=0.d0
            q3_z(zstart(1):zend(1), zstart(2):zend(2), n3) =0.d0

        end if

        ! ATTENTION
        if (BC3==FREESLIP) then

            q3_wall30(zstart(1):zend(1), zstart(2):zend(2))=0.d0                                             ! No penetration
            q2_wall30(zstart(1):zend(1), zstart(2):zend(2))=q2_z(zstart(1):zend(1), zstart(2):zend(2), 1)     ! Neumann
            q1_wall30(zstart(1):zend(1), zstart(2):zend(2))=q1_z(zstart(1):zend(1), zstart(2):zend(2), 1)     ! Neumann

            q3_wall31(zstart(1):zend(1), zstart(2):zend(2))=0.d0                                             ! No penetration
            q2_wall31(zstart(1):zend(1), zstart(2):zend(2))=q2_z(zstart(1):zend(1), zstart(2):zend(2), n3-1)  ! Neumann
            q1_wall31(zstart(1):zend(1), zstart(2):zend(2))=q1_z(zstart(1):zend(1), zstart(2):zend(2), n3-1)  ! Neumann

        end if


        return
    end subroutine


    subroutine apply_BC3_ibm

        use IBM_data
        use mesh 
        use boundaries

        implicit none

        ! ATTENTION
        if (BC3==NOSLIP) then

            q3_wall30_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2))=0.d0
            q2_wall30_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2))=0.d0
            q1_wall30_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2))=0.d0

            if (k_start_ibm_2nd_fc.eq.1) q3_z_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2), 1) = 0.d0

            q3_wall31_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2))=0.d0
            q2_wall31_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2))=0.d0
            q1_wall31_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2))=0.d0
            if (k_end_ibm_2nd_fc.eq.n3_ibm) q3_z_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2), n3_ibm) = 0.d0

        end if

        ! ATTENTION
        if (BC3==FREESLIP) then

            q3_wall30_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2))=0.d0                                             ! No penetration
            if (k_start_ibm_2nd_fc.eq.1) q2_wall30_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2))=q2_z_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2), 1)     ! Neumann
            if (k_start_ibm_2nd_fc.eq.1) q1_wall30_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2))=q1_z_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2), 1)     ! Neumann

            q3_wall31_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2))=0.d0                                             ! No penetration
            if (k_end_ibm_2nd_fc.ge.(n3_ibm-1)) q2_wall31_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2))=q2_z_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2), n3_ibm-1)  ! Neumann
            if (k_end_ibm_2nd_fc.ge.(n3_ibm-1)) q1_wall31_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2))=q1_z_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2), n3_ibm-1)  ! Neumann

        end if


        return
    end subroutine

    subroutine apply_BC2

        use mesh
        use physical_fields
        use boundaries
        use decomp_2d

        use run_ctxt_data, only: ntime, t
        use blow_settings
        use twave_settings
        use DNS_settings

        implicit none

        integer :: k,i, s
        integer :: s_xst, s_xen, s_zst, s_zen
        real*8  :: amp,kappa,omega
        real*8  :: cf_amp,cf_kappa,cf_omega


        if (BC2==NOSLIP) then

            q3_wall20(ystart(1):yend(1), ystart(3):yend(3)) =0.d0
            q2_wall20(ystart(1):yend(1), ystart(3):yend(3)) =0.d0
            q2_y(ystart(1):yend(1), 1, ystart(3):yend(3))   =0.d0
            q1_wall20(ystart(1):yend(1), ystart(3):yend(3)) =0.d0

            q3_wall21(ystart(1):yend(1), ystart(3):yend(3)) =0.d0
            q2_wall21(ystart(1):yend(1), ystart(3):yend(3)) =0.d0
            q2_y(ystart(1):yend(1), n2, ystart(3):yend(3))  =0.d0
            q1_wall21(ystart(1):yend(1), ystart(3):yend(3)) =0.d0

            ! Blowing: At this time, blowing is only in the Y direction -----------------
            do s = 1, nb_slows

                s_xst=slows(s)%xst*(n3m-1)+1
                s_xen=slows(s)%xen*(n3m-1)+1
                s_zst=slows(s)%zst*(n1m-1)+1
                s_zen=slows(s)%zen*(n1m-1)+1

                do k=ystart(3),min(yend(3), n3m)
                    do i=ystart(1),min(yend(1), n1m)

                        if ((k>=s_xst).and.(k<=s_xen).and.(i>=s_zst).and.(i<=s_zen)) then

                            q2_wall20(i,k)  = slows(s)%blowing
                            q2_wall21(i,k)  = q2_wall20(i,k)

                            q2_y(i, 1, k)  =q2_wall20(i, k)
                            q2_y(i, n2, k) =q2_wall21(i, k)

                        end if

                    enddo
                enddo

            end do

            if(ntime.eq.blow_end) then
                call disable_blowing
            endif

        end if

        !Travelling wave Boundary condition for spanwise wall velocity
        if(twave_on==1) then

            ! To convert from inner to outer units in case travelling
            ! wave parameters are given in inner units
            if(inner_units==1) then
                cf_amp   = Re_tau/ren
                cf_omega = Re_tau**2/ren    
                cf_kappa = Re_tau
            else if(inner_units==0) then
                cf_amp   = 1.0d0
                cf_omega = 1.0d0    
                cf_kappa = 1.0d0
            end if
                
            amp   = Twave%Amp  *cf_amp
            kappa = Twave%kappa*cf_kappa
            omega = Twave%omega*cf_omega

            if(tstart .gt. t) then
                write(*,*) "****************ERROR*****************"
                write(*,*) "      tstart is greater than t        "
                write(*,*) " check the travelling wave input file "
                write(*,*) "            ABORTING CODE             "
                write(*,*) "**************************************"
                STOP
            end if

            ! Homogoneous spanwise wall oscillations W=A*sin(omega*t)
            if(kappa==0.0d0 .and. omega/=0.0d0) then    

                if(streamwise==1) then
                    do k=ystart(3),yend(3)
                        do i=ystart(1),yend(1)
                            q3_wall20(i,k) = amp*dsin(omega*(t-tstart)) !tstart is subtracted to start the oscillations from W=0
                        end do
                    end do
                            q3_wall21 = q3_wall20
                elseif(streamwise==3) then
                    do k=ystart(3),yend(3)
                        do i=ystart(1),yend(1)
                            q1_wall20(i,k) = amp*dsin(omega*(t-tstart)) !tstart is subtracted to start the oscillations from W=0
                        end do
                    end do
                            q1_wall21 = q1_wall20
                end if

            ! Homogeneous spanwise steady oscillations W=A*sin(kappa*x)
            elseif(kappa/=0.0d0 .and. omega==0.0d0) then    

                if(streamwise==1) then
                    do k=ystart(3),yend(3)
                        do i=ystart(1),yend(1)
                            q3_wall20(i,k) = amp*dsin(kappa*z(i))
                        end do
                    end do
                            q3_wall21 = q3_wall20
                elseif(streamwise==3) then
                    do k=ystart(3),yend(3)
                        do i=ystart(1),yend(1)
                            q1_wall20(i,k) = amp*dsin(kappa*x(k))
                        end do
                    end do
                            q1_wall21 = q1_wall20
                end if

            ! Travelling wave of spanwise wall velocity modulated in streamwise direction
            ! W=A*sin(kappa*x-omega*t)
            elseif(kappa/=0.0d0 .and. omega/=0.0d0) then    

                if(streamwise==1) then
                    do k=ystart(3),yend(3)
                        do i=ystart(1),yend(1)
                            q3_wall20(i,k) = amp*dsin(kappa*z(i)-omega*(t)) !tstart is subtracted to start the oscillations from W=0
                        end do
                    end do
                            q3_wall21 = q3_wall20
                elseif(streamwise==3) then
                    do k=ystart(3),yend(3)
                        do i=ystart(1),yend(1)
                            q1_wall20(i,k) = amp*dsin(kappa*x(k)-omega*(t)) !tstart is subtracted to start the oscillations from W=0
                        end do
                    end do
                            q1_wall21 = q1_wall20
                end if

            end if

        end if

    ! ATTENTION
    if (BC2==FREESLIP) then

      q3_wall20(ystart(1):yend(1), ystart(3):yend(3))=q3_y(ystart(1):yend(1), 1, ystart(3):yend(3))       ! Neumann
      q2_wall20(ystart(1):yend(1), ystart(3):yend(3))=0.d0                                                ! No penetration
      q1_wall20(ystart(1):yend(1), ystart(3):yend(3))=q1_y(ystart(1):yend(1), 1, ystart(3):yend(3))       ! Neumann

      q3_wall21(ystart(1):yend(1), ystart(3):yend(3))=q3_y(ystart(1):yend(1), n2-1, ystart(3):yend(3))    ! Neumann
      q2_wall21(ystart(1):yend(1), ystart(3):yend(3))=0.d0                                                ! No penetration
      q1_wall21(ystart(1):yend(1), ystart(3):yend(3))=q1_y(ystart(1):yend(1), n2-1, ystart(3):yend(3))    ! Neumann

     end if


        return
    end subroutine

    subroutine apply_BC2_ibm

        use IBM_data
        use mesh 
        use boundaries

        implicit none

        integer :: k,i, s
        integer :: s_xst, s_xen, s_zst, s_zen
        real*8  :: amp,kappa,omega
        real*8  :: cf_amp,cf_kappa,cf_omega


        if (BC2==NOSLIP) then

            q3_wall20_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3)) =0.d0
            q2_wall20_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3)) =0.d0
            if (j_start_ibm_2nd_fc.eq.1) q2_y_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), 1, ystart_ibm_fcm(3):yend_ibm_fcm(3))   =0.d0
            q1_wall20_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3)) =0.d0

            q3_wall21_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3)) =0.d0
            q2_wall21_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3)) =0.d0
            if (j_end_ibm_2nd_fc.eq.n2_ibm) q2_y_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), n2_ibm, ystart_ibm_fcm(3):yend_ibm_fcm(3))  =0.d0
            q1_wall21_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3)) =0.d0

        endif

        ! ATTENTION
        if (BC2==FREESLIP) then

            if (j_start_ibm_2nd_fc.eq.1) q3_wall20_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3))=q3_y_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), 1, ystart_ibm_fcm(3):yend_ibm_fcm(3))       ! Neumann
            q2_wall20_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3))=0.d0                                                ! No penetration
            if (j_start_ibm_2nd_fc.eq.1) q1_wall20_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3))=q1_y_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), 1, ystart_ibm_fcm(3):yend_ibm_fcm(3))       ! Neumann

            if (j_end_ibm_2nd_fc.ge.(n2_ibm-1)) q3_wall21_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3))=q3_y_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), n2_ibm-1, ystart_ibm_fcm(3):yend_ibm_fcm(3))    ! Neumann
            q2_wall21_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3))=0.d0                                                ! No penetration
            if (j_end_ibm_2nd_fc.ge.(n2_ibm-1)) q1_wall21_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3))=q1_y_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), n2_ibm-1, ystart_ibm_fcm(3):yend_ibm_fcm(3))    ! Neumann

        end if


        return
    end subroutine

    subroutine apply_BC1

        use mesh
        use physical_fields
        use boundaries
        use decomp_2d

        implicit none
        integer :: k,i,j, s
        real*8, dimension(xstart(2):xend(2), xstart(3):xend(3))     :: q1w_fluc, q2w_fluc, q3w_fluc
        real*8, parameter                                           :: noise1=0.1d0 !ATTENTION
        integer, parameter  :: order=0
        real*8  :: q1o(-1:1), q2o(-1:1), q3o(-1:1), df,d2f

        ! ATTENTION
        if (BC1==NOSLIP) then

            q3_wall10(xstart(2):xend(2), xstart(3):xend(3)) =0.d0
            q2_wall10(xstart(2):xend(2), xstart(3):xend(3)) =0.d0
            q1_wall10(xstart(2):xend(2), xstart(3):xend(3)) =0.d0
            q1_x(1, xstart(2):xend(2), xstart(3):xend(3))   =0.d0

            q3_wall11(xstart(2):xend(2), xstart(3):xend(3)) =0.d0
            q2_wall11(xstart(2):xend(2), xstart(3):xend(3)) =0.d0
            q1_wall11(xstart(2):xend(2), xstart(3):xend(3)) =0.d0
            q1_x(n1, xstart(2):xend(2), xstart(3):xend(3))  =0.d0

        end if

        if (BC1==FREESLIP) then

            q3_wall10(xstart(2):xend(2), xstart(3):xend(3))=q3_x(1, xstart(2):xend(2), xstart(3):xend(3))     ! Neumann
            q2_wall10(xstart(2):xend(2), xstart(3):xend(3))=q2_x(1, xstart(2):xend(2), xstart(3):xend(3))     ! Neumann
            q1_wall10(xstart(2):xend(2), xstart(3):xend(3))=0.d0                                             ! No penetration

            q3_wall11(xstart(2):xend(2), xstart(3):xend(3))=q3_x(n1-1, xstart(2):xend(2), xstart(3):xend(3))  ! Neumann
            q2_wall11(xstart(2):xend(2), xstart(3):xend(3))=q2_x(n1-1, xstart(2):xend(2), xstart(3):xend(3))  ! Neumann
            q1_wall11(xstart(2):xend(2), xstart(3):xend(3))=0.d0                                             ! No penetration

        end if


        return
    end subroutine

    subroutine apply_BC1_ibm

        use IBM_data
        use mesh 
        use boundaries

        implicit none

        ! ATTENTION
        if (BC1==NOSLIP) then

            q3_wall10_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)) =0.d0
            q2_wall10_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)) =0.d0
            q1_wall10_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)) =0.d0
            if (i_start_ibm_2nd_fc.eq.1) q1_x_ibm(1, xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3))   =0.d0

            q3_wall11_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)) =0.d0
            q2_wall11_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)) =0.d0
            q1_wall11_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)) =0.d0
            if (i_end_ibm_2nd_fc.eq.n1_ibm) q1_x_ibm(n1_ibm, xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3))  =0.d0

        end if

        if (BC1==FREESLIP) then

            if (i_start_ibm_2nd_fc.eq.1) q3_wall10_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3))=q3_x_ibm(1, xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3))     ! Neumann
            if (i_start_ibm_2nd_fc.eq.1) q2_wall10_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3))=q2_x_ibm(1, xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3))     ! Neumann
            q1_wall10_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3))=0.d0                                             ! No penetration

            if (i_end_ibm_2nd_fc.ge.(n1_ibm-1)) q3_wall11_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3))=q3_x_ibm(n1_ibm-1, xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3))  ! Neumann
            if (i_end_ibm_2nd_fc.ge.(n1_ibm-1)) q2_wall11_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3))=q2_x_ibm(n1_ibm-1, xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3))  ! Neumann
            q1_wall11_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3))=0.d0                                             ! No penetration

        end if


        return
    end subroutine

    subroutine disable_blowing

        use boundaries

        implicit none

        slows(:)%blowing=0.d0

        return
    end subroutine

end module VELOCITY_bc_controller
