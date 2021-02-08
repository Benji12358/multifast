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

    subroutine apply_BC2

        use mesh
        use physical_fields
        use boundaries
        use decomp_2d

        use run_ctxt_data, only: ntime
        use blow_settings

        implicit none

        integer :: k,i, s
        integer :: s_xst, s_xen, s_zst, s_zen


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

        if (BC1==FRINGE) then

            q3_wall10(xstart(2):xend(2), xstart(3):xend(3))=q3_x(1, xstart(2):xend(2), xstart(3):xend(3))     ! Neumann
            q2_wall10(xstart(2):xend(2), xstart(3):xend(3))=q2_x(1, xstart(2):xend(2), xstart(3):xend(3))     ! Neumann
            q1_wall10(xstart(2):xend(2), xstart(3):xend(3))=q1_x(1, xstart(2):xend(2), xstart(3):xend(3))     ! No penetration

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
