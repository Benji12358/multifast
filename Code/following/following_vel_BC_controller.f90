module following_velocity_bc_controller
    use following_data
    implicit none

contains

    subroutine apply_BC3

        use following_mesh
        use following_physical_fields
        use boundaries

        implicit none

        ! ATTENTION
        if (BC3==NOSLIP) then

            q3_wall30(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2))=0.d0
            q2_wall30(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2))=0.d0
            q1_wall30(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2))=0.d0
            q3_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), 1) =0.d0

            q3_wall31(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2))=0.d0
            q2_wall31(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2))=0.d0
            q1_wall31(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2))=0.d0
            q3_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), n3) =0.d0

        end if

        ! ATTENTION
        if (BC3==FREESLIP) then

            q3_wall30(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2))=0.d0                                             ! No penetration
            q2_wall30(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2))=q2_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), 1)     ! Neumann
            q1_wall30(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2))=q1_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), 1)     ! Neumann

            q3_wall31(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2))=0.d0                                             ! No penetration
            q2_wall31(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2))=q2_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), n3-1)  ! Neumann
            q1_wall31(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2))=q1_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), n3-1)  ! Neumann

        end if


        return
    end subroutine

    subroutine apply_BC2

        use following_mesh
        use following_physical_fields
        use boundaries

        implicit none


        if (BC2==NOSLIP) then

            q3_wall20(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3)) =0.d0
            q2_wall20(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3)) =0.d0
            q2_y(ystart_f(1):yend_f(1), 1, ystart_f(3):yend_f(3))   =0.d0
            q1_wall20(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3)) =0.d0

            q3_wall21(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3)) =0.d0
            q2_wall21(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3)) =0.d0
            q2_y(ystart_f(1):yend_f(1), n2, ystart_f(3):yend_f(3))  =0.d0
            q1_wall21(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3)) =0.d0

        end if

        ! ATTENTION
        if (BC2==FREESLIP) then

          q3_wall20(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3))=q3_y(ystart_f(1):yend_f(1), 1, ystart_f(3):yend_f(3))       ! Neumann
          q2_wall20(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3))=0.d0                                                ! No penetration
          q1_wall20(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3))=q1_y(ystart_f(1):yend_f(1), 1, ystart_f(3):yend_f(3))       ! Neumann

          q3_wall21(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3))=q3_y(ystart_f(1):yend_f(1), n2-1, ystart_f(3):yend_f(3))    ! Neumann
          q2_wall21(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3))=0.d0                                                ! No penetration
          q1_wall21(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3))=q1_y(ystart_f(1):yend_f(1), n2-1, ystart_f(3):yend_f(3))    ! Neumann

         end if


        return
    end subroutine

    subroutine apply_BC1

        use following_mesh
        use following_physical_fields
        use boundaries

        implicit none

        ! ATTENTION
        if (BC1==NOSLIP) then

            q3_wall10(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)) =0.d0
            q2_wall10(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)) =0.d0
            q1_wall10(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)) =0.d0
            q1_x(1, xstart_f(2):xend_f(2), xstart_f(3):xend_f(3))   =0.d0

            q3_wall11(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)) =0.d0
            q2_wall11(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)) =0.d0
            q1_wall11(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)) =0.d0
            q1_x(n1, xstart_f(2):xend_f(2), xstart_f(3):xend_f(3))  =0.d0

        end if

        if (BC1==FREESLIP) then

            q3_wall10(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3))=q3_x(1, xstart_f(2):xend_f(2), xstart_f(3):xend_f(3))     ! Neumann
            q2_wall10(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3))=q2_x(1, xstart_f(2):xend_f(2), xstart_f(3):xend_f(3))     ! Neumann
            q1_wall10(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3))=0.d0                                             ! No penetration

            q3_wall11(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3))=q3_x(n1-1, xstart_f(2):xend_f(2), xstart_f(3):xend_f(3))  ! Neumann
            q2_wall11(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3))=q2_x(n1-1, xstart_f(2):xend_f(2), xstart_f(3):xend_f(3))  ! Neumann
            q1_wall11(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3))=0.d0                                             ! No penetration

        end if


        return
    end subroutine

end module following_velocity_bc_controller