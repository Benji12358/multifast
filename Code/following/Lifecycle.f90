module following_velocity_init
    implicit none

contains

    subroutine load_initial_state()

        use following_start_settings

        use following_common_workspace_view
        use following_velocity_loader
        use following_velocity_bc_controller
        use following_turbulence_generator

        use run_ctxt_data

        implicit none
        logical             :: fexist(4)
        logical             :: previousRHS_are_available=.false.
        logical             :: start_from_coarse_file=.false.

        character*200       :: VELOCITY_external_fields_path
        character*20        :: filename
        character*10 tmp_str

        fexist=.true.

        select case (start_source_type)
            case (NO_SOURCE)

                call init_turbulent_field(vper)

                if(nrank==0) write(*,*) 'Running a new simulation from a generated turbulent field' , vper
                if(nrank==0) write(*,*) 'The Euler scheme will be used for the first iteration'

                call apply_BC1
                call apply_BC2

            case (HDF5_FILE)

                write(tmp_str, "(i10)")start_it
                filename="field"//trim(adjustl(tmp_str))

                VELOCITY_external_fields_path=trim(external_fields_path)//"/following/"//trim(filename)

                if(nrank==0) write(*,*) 'VELOCITY_external_fields_path ', trim(VELOCITY_external_fields_path)
                if(nrank==0) write(*,*) 'start_from_coarse_file ', start_from_coarse_file

                call load_fields(VELOCITY_external_fields_path, start_from_coarse_file, fexist)

        end select

        call transpose_y_to_x(q1_y, q1_x, decomp_following)
        call transpose_y_to_z(q1_y, q1_z, decomp_following)

        call transpose_y_to_x(q2_y, q2_x, decomp_following)
        call transpose_y_to_z(q2_y, q2_z, decomp_following)

        call transpose_y_to_x(q3_y, q3_x, decomp_following)
        call transpose_y_to_z(q3_y, q3_z, decomp_following)


    end subroutine load_initial_state

    subroutine allocate_data()

        use decomp_2d
        use following_physical_fields
        use following_data

        implicit none

        ! Inner values (in x,y,z decomposition configuration)-----------------
        allocate(q3_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        q3_x=0.d0
        allocate(q2_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        q2_x=0.d0
        allocate(q1_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        q1_x=0.d0
        allocate(divu_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        divu_x=0.d0
        allocate(dp_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        dp_x=0.d0

        allocate(q3_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
        q3_y=0.d0
        allocate(q2_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
        q2_y=0.d0
        allocate(q1_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
        q1_y=0.d0
        allocate(divu_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
        divu_y=0.d0
        allocate(dp_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
        dp_y=0.d0

        allocate(q3_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
        q3_z=0.d0
        allocate(q2_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
        q2_z=0.d0
        allocate(q1_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
        q1_z=0.d0

        allocate(divu_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
        divu_z=0.d0
        allocate(dp_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
        dp_z=0.d0
        allocate(source_term(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
        source_term=0.d0

        allocate(dphidx1_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        dphidx1_x=0.d0
        allocate(dphidx2_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        dphidx2_x=0.d0
        allocate(dphidx3_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        dphidx3_x=0.d0

        allocate(dphidx1_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
        dphidx1_y=0.d0
        allocate(dphidx2_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
        dphidx2_y=0.d0
        allocate(dphidx3_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
        dphidx3_y=0.d0

        allocate(dphidx1_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
        dphidx1_z=0.d0
        allocate(dphidx2_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
        dphidx2_z=0.d0
        allocate(dphidx3_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
        dphidx3_z=0.d0

        ! Wall values allocation ---------------------------------------------
        !ATTENTION
        allocate(q3_wall10(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        q3_wall10=0.d0
        allocate(q2_wall10(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        q2_wall10=0.d0
        allocate(q1_wall10(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        q1_wall10=0.d0

        allocate(q3_wall11(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        q3_wall11=0.d0
        allocate(q2_wall11(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        q2_wall11=0.d0
        allocate(q1_wall11(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        q1_wall11=0.d0

        allocate(q3_wall20(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3)))
        q3_wall20=0.d0
        allocate(q2_wall20(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3)))
        q2_wall20=0.d0
        allocate(q1_wall20(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3)))
        q1_wall20=0.d0

        allocate(q3_wall21(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3)))
        q3_wall21=0.d0
        allocate(q2_wall21(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3)))
        q2_wall21=0.d0
        allocate(q1_wall21(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3)))
        q1_wall21=0.d0

        !ATTENTION
        allocate(q3_wall30(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2)))
        q3_wall30=0.d0
        allocate(q2_wall30(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2)))
        q2_wall30=0.d0
        allocate(q1_wall30(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2)))
        q1_wall30=0.d0

        allocate(q3_wall31(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2)))
        q3_wall31=0.d0
        allocate(q2_wall31(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2)))
        q2_wall31=0.d0
        allocate(q1_wall31(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2)))
        q1_wall31=0.d0


    end subroutine allocate_data

    subroutine initialize()

        use following_velocity_bc_controller
        use following_velocity_solver, only: SOLVER_init => init

        implicit none

        call allocate_data

        call SOLVER_init

        call load_initial_state

        contains

    end subroutine initialize

end module following_velocity_init


module following_velocity_final
    implicit none

contains

    subroutine deallocate_data()

        use following_physical_fields
        use following_velocity_solver

        implicit none

        !!! PHYSICAL FIELDS !!!!
        ! Inner values (in x,y,z decomposition configuration)-----------------
        deallocate(q3_x)
        deallocate(q2_x)
        deallocate(q1_x)
        deallocate(divu_x)
        deallocate(dp_x)

        deallocate(q3_y)
        deallocate(q2_y)
        deallocate(q1_y)
        deallocate(divu_y)
        deallocate(dp_y)

        deallocate(q3_z)
        deallocate(q2_z)
        deallocate(q1_z)

        deallocate(divu_z)
        deallocate(dp_z)

        deallocate(dphidx1_x)
        deallocate(dphidx2_y)
        deallocate(dphidx2_x)
        deallocate(dphidx2_z)
        deallocate(dphidx3_z)
        deallocate(dphidx3_x)

        ! Wall values deallocation ---------------------------------------------
        deallocate(q3_wall10)
        deallocate(q2_wall10)
        deallocate(q1_wall10)

        deallocate(q3_wall11)
        deallocate(q2_wall11)
        deallocate(q1_wall11)

        deallocate(q3_wall20)
        deallocate(q2_wall20)
        deallocate(q1_wall20)

        deallocate(q3_wall21)
        deallocate(q2_wall21)
        deallocate(q1_wall21)

        deallocate(q3_wall30)
        deallocate(q2_wall30)
        deallocate(q1_wall30)

        deallocate(q3_wall31)
        deallocate(q2_wall31)
        deallocate(q1_wall31)

    end subroutine deallocate_data

    subroutine finalize()

        use following_velocity_solver, SOLVER_finalize=>finalize
        implicit none

        call SOLVER_finalize
        call deallocate_data

    end subroutine finalize

end module following_velocity_final

module following_fringe_init
    implicit none

contains

    subroutine allocate_data()

        use decomp_2d
        use following_data
        use following_fringe_data
        implicit none

        allocate(f1_fringe_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        allocate(f2_fringe_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        allocate(f3_fringe_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))

        allocate(sca_inflow(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        allocate(q1_inflow(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        allocate(q2_inflow(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        allocate(q3_inflow(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))

        allocate(lambda_x(xstart_f(1):xend_f(1)))        

        f1_fringe_x = 0.d0
        f2_fringe_x = 0.d0
        f3_fringe_x = 0.d0

        sca_inflow = 0.d0
        q1_inflow = 0.d0
        q2_inflow = 0.d0
        q3_inflow = 0.d0

        lambda_x = 0.d0

    end subroutine allocate_data

    subroutine initialize()
        use following_fringe_data
        use following_fringe_initializer
        implicit none

        call allocate_data

        call get_inflow()

    end subroutine initialize

end module following_fringe_init


module following_fringe_final
    implicit none

contains

    subroutine deallocate_data()

        use following_fringe_data

        implicit none

        deallocate(f1_fringe_x)
        deallocate(f2_fringe_x)
        deallocate(f3_fringe_x)

        deallocate(lambda_x)

    end subroutine deallocate_data

    subroutine finalize()
        implicit none

        call deallocate_data

    end subroutine finalize

end module following_fringe_final

module following_scalar_init
    implicit none

contains

    subroutine load_initial_state()

        use run_ctxt_data

        use following_start_settings

        use following_common_workspace_view
        use following_scalar_loader
        use following_scalar_field_generator

        implicit none
        logical             :: fexist(4),start_from_coarse_file

        character*200       :: SCALAR_external_fields_path
        character*20        :: filename
        character*10 tmp_str

        fexist(4)=.false.
        start_from_coarse_file=.false.

        if ((start_source_type==NO_SOURCE).or.(reset_scalar_field)) then

            call generate_fields(sca_x(:,:,:), sca_y(:,:,:), sca_z(:,:,:), -delta_T, delta_T)

        else if (start_source_type==HDF5_FILE) then

            write(tmp_str, "(i10)")start_it
            filename="field"//trim(adjustl(tmp_str))

            SCALAR_external_fields_path=trim(external_fields_path)//"/following/"//trim(filename)

            if(nrank==0) write(*,*) 'SCALAR_external_fields_path ', trim(SCALAR_external_fields_path)

            call load_fields(SCALAR_external_fields_path, start_from_coarse_file, fexist)
            if (nrank==0) write(*,*) 'bef init_wall_sca'
            call init_wall_sca()

        endif

    end subroutine load_initial_state

    subroutine allocate_data()

        use decomp_2d
        use following_data
        use following_scalar_data
        implicit none


        allocate(sca_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        sca_x=0.d0

        allocate(sca_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
        sca_y=0.d0

        allocate(sca_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
        sca_z=0.d0

        ! Wall values allocation ---------------------------------------------
        allocate(sca_wall10(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        sca_wall10=0.d0
        allocate(sca_wall11(xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
        sca_wall11=0.d0

        allocate(sca_wall20(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3)))
        sca_wall20=0.d0
        allocate(sca_wall21(ystart_f(1):yend_f(1), ystart_f(3):yend_f(3)))
        sca_wall21=0.d0

        allocate(sca_wall30(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2)))
        sca_wall30=0.d0
        allocate(sca_wall31(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2)))
        sca_wall31=0.d0

    end subroutine allocate_data

    subroutine initialize()

        use following_scalar_solver, only: SOLVER_init => init

        implicit none

        call allocate_data

        call SOLVER_init

        call load_initial_state

    end subroutine initialize

end module following_scalar_init


module following_scalar_final
    implicit none

contains

    subroutine deallocate_data()

        use following_scalar_data
        use following_scalar_solver, only: previousRHS1, previousRHS2

        implicit none

        deallocate(sca_x)
        deallocate(sca_y)
        deallocate(sca_z)

        deallocate(sca_wall10)
        deallocate(sca_wall11)

        deallocate(sca_wall20)
        deallocate(sca_wall21)

        deallocate(sca_wall30)
        deallocate(sca_wall31)

        deallocate(previousRHS1)
        deallocate(previousRHS2)

    end subroutine deallocate_data

    subroutine finalize()
        implicit none

        call deallocate_data

    end subroutine finalize

end module following_scalar_final
