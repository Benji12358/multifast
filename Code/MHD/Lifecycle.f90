module MHD_init
    implicit none

contains

!    subroutine read_settings()
!        implicit none
!    end subroutine read_settings

    subroutine allocate_data()

        use decomp_2d

        use MHD_data
        implicit none

            allocate(fb1_MHD_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            allocate(fb2_MHD_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            allocate(fb3_MHD_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
!            allocate(fb1_MHD_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
!            allocate(fb2_MHD_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
!            allocate(fb3_MHD_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
!            allocate(fb3_MHD_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

            allocate(phi_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            allocate(phi_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            allocate(phi_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

            allocate(gradphi1_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            allocate(gradphi2_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            allocate(gradphi3_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

            allocate(ucrossB1_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
!            allocate(ucrossB1_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            allocate(ucrossB2_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
!            allocate(ucrossB3_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            allocate(ucrossB3_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

            allocate(B01_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
!            allocate(B01_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            allocate(B02_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
!            allocate(B03_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            allocate(B01c_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            allocate(B02c_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            allocate(B03c_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            allocate(B03_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

            allocate(A1_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            allocate(A1_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            allocate(A2_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            allocate(A3_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            allocate(A3_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

            allocate(RHS_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
            allocate(RHS_z2(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))


        fb1_MHD_x = 0.d0
        fb2_MHD_x = 0.d0
        fb3_MHD_x = 0.d0

 !       fb1_MHD_y = 0.d0
 !       fb2_MHD_y = 0.d0
 !       fb3_MHD_y = 0.d0
 !       fb3_MHD_z = 0.d0

        phi_x = 0.d0
        phi_y = 0.d0
        phi_z = 0.d0

        gradphi1_x = 0.d0
        gradphi2_y = 0.d0
        gradphi3_z = 0.d0

        ucrossB1_x = 0.d0
 !       ucrossB1_y = 0.d0
        ucrossB2_y = 0.d0
!        ucrossB3_y = 0.d0
        ucrossB3_z = 0.d0

        B01_x = 0.d0 !Magnetic_field_unit_vector_x
!       B01_y = 0.d0 !Magnetic_field_unit_vector_x
        B02_y = 0.d0 !Magnetic_field_unit_vector_y
!        B03_y = 0.d0 !Magnetic_field_unit_vector_z
        B03_z = 0.d0 !Magnetic_field_unit_vector_z
        B01c_y = 0.d0
        B02c_y = 0.d0
        B03c_y = 0.d0

        A1_x = 0.d0
        A1_y = 0.d0
        A2_y = 0.d0
        A3_y = 0.d0
        A3_z = 0.d0

        RHS_z = 0.d0
        RHS_z2 = 0.d0


    end subroutine allocate_data

    subroutine initialize()
        use MHD_data
        use MHD_field_generator
        implicit none

        call allocate_data

        call generate_total_magnetic_field(B01_x, B02_y, B03_z, B01c_y, B02c_y, B03c_y)
!
! **** A ACtiver pour la sortie des Resultats ****
!
        contains
!
            subroutine init_workspace()
!
                use mpi
!
!
                use file_copy
!
                use COMMON_workspace_view
                use MHD_workspace_view
!
                implicit none
!
!
                results_path            =trim(COMMON_results_path)
                results3D_path          =trim(results_path)//"3D/"
!
!                sensor_file             =trim(COMMON_log_path)//'Scalar/AA_v.csv'

!
            end subroutine init_workspace

    end subroutine initialize

end module MHD_init


module MHD_final
    implicit none

contains

    subroutine deallocate_data()

        use MHD_data

        implicit none

        deallocate(fb1_MHD_x)
        deallocate(fb2_MHD_x)
        deallocate(fb3_MHD_x)

        deallocate(phi_x)
        deallocate(phi_y)
        deallocate(phi_z)

        deallocate(gradphi1_x)
        deallocate(gradphi2_y)
        deallocate(gradphi3_z)

        deallocate(ucrossB1_x)
        deallocate(ucrossB2_y)
        deallocate(ucrossB3_z)

        deallocate(B01_x)
        deallocate(B02_y)
        deallocate(B01c_y)
        deallocate(B02c_y)
        deallocate(B03c_y)
        deallocate(B03_z)

        deallocate(A1_x)
        deallocate(A1_y)
        deallocate(A2_y)
        deallocate(A3_y)
        deallocate(A3_z)

        deallocate(RHS_z)
        deallocate(RHS_z2)


    end subroutine deallocate_data

    subroutine finalize()
        implicit none

        call deallocate_data

    end subroutine finalize

end module MHD_final
