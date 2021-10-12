module FRINGE_init
    implicit none

contains

    subroutine allocate_data()

        use decomp_2d
        use DNS_settings
        use FRINGE_data
        implicit none

        allocate(f1_fringe_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        allocate(f2_fringe_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        allocate(f3_fringe_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))

        allocate(f1_fringe_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
        allocate(f2_fringe_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
        allocate(f3_fringe_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

        if (streamwise==1) then
            allocate(q1_inflow(xstart(2):xend(2), xstart(3):xend(3)))
            allocate(q2_inflow(xstart(2):xend(2), xstart(3):xend(3)))
            allocate(q3_inflow(xstart(2):xend(2), xstart(3):xend(3)))
        endif

        if (streamwise==3) then
            allocate(q1_inflow(zstart(1):zend(1), zstart(2):zend(2)))
            allocate(q2_inflow(zstart(1):zend(1), zstart(2):zend(2)))
            allocate(q3_inflow(zstart(1):zend(1), zstart(2):zend(2)))
        endif

        allocate(lambda_x(xstart(1):xend(1)))
        allocate(lambda_z(zstart(3):zend(3)))
        

        f1_fringe_x = 0.d0
        f2_fringe_x = 0.d0
        f3_fringe_x = 0.d0

        f1_fringe_z = 0.d0
        f2_fringe_z = 0.d0
        f3_fringe_z = 0.d0

        q1_inflow = 0.d0
        q2_inflow = 0.d0
        q3_inflow = 0.d0

        lambda_x = 0.d0
        lambda_z = 0.d0

    end subroutine allocate_data

    subroutine initialize()
        use FRINGE_data
        use FRINGE_initializer
        implicit none

        call allocate_data

        call get_inflow()
!
        ! contains
!
!             subroutine init_workspace()
! !
!                 use mpi
! !
! !
!                 use file_copy
! !
!                 use COMMON_workspace_view
!                 use MHD_workspace_view
! !
!                 implicit none
! !
! !
!                 results_path            =trim(COMMON_results_path)
!                 results3D_path          =trim(results_path)//"3D/"
!
!
!            end subroutine init_workspace

    end subroutine initialize

end module FRINGE_init


module FRINGE_final
    implicit none

contains

    subroutine deallocate_data()

        use FRINGE_data

        implicit none

        deallocate(f1_fringe_x)
        deallocate(f2_fringe_x)
        deallocate(f3_fringe_x)

        deallocate(q1_inflow)
        deallocate(q2_inflow)
        deallocate(q3_inflow)

        deallocate(lambda_x)


    end subroutine deallocate_data

    subroutine finalize()
        implicit none

        call deallocate_data

    end subroutine finalize

end module FRINGE_final
