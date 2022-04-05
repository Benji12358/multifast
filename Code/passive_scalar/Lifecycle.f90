module SCALAR_init
    implicit none

contains

    subroutine read_settings()
        implicit none
    end subroutine read_settings

    subroutine load_initial_state()

        use run_ctxt_data

        use start_settings

        use SCALAR_workspace_view
        use SCALAR_loader
        use SCALAR_field_generator

        implicit none
        logical             :: fexist(4)

        character*200       :: SCALAR_external_fields_path
        character*20        :: filename
        character*10 tmp_str

        fexist(4)=.false.

        if ((run_ctxt==CONTINUE_FROM_PREVIOUS_RUN).or.(run_ctxt==RECOVERY_A_RUN)) then
            write(*,*) 'here'
            fexist=.true.

            call load_fields(recovery_fields_dir, .false., fexist)
            if (nrank==0) write(*,*)'fexist', fexist
            if (.not. fexist(4)) then
                call generate_fields(sca_x(:,:,:), sca_y(:,:,:), sca_z(:,:,:), -delta_T, delta_T)
            endif

        else

            if ((start_source_type==NO_SOURCE).or.(reset_scalar_field)) then

            ! select case (start_source_type)
            !     case (NO_SOURCE)

                if ((init_type==CLASSIC_INIT).or.(init_type==KAWAMURA_INIT).or.(init_type==CONSTANT_HEAT_FLUX)) then

                    call generate_fields(sca_x(:,:,:), sca_y(:,:,:), sca_z(:,:,:), -delta_T, delta_T)

                elseif (init_type==INIT_FROM_FILE) then

                    call generate_from_file(sca_x(:,:,:), sca_y(:,:,:), sca_z(:,:,:))

                endif

            else if (start_source_type==HDF5_FILE) then

                ! case (HDF5_FILE)

                write(tmp_str, "(i10)")start_it
                filename="field"//trim(adjustl(tmp_str))

                SCALAR_external_fields_path=trim(external_fields_path)//"/3D/"//trim(filename)

                if(nrank==0) write(*,*) 'SCALAR_external_fields_path ', trim(SCALAR_external_fields_path)


                call load_fields(SCALAR_external_fields_path, start_from_coarse_file, fexist)

                if (.not. fexist(4)) then
                    call generate_fields(sca_x(:,:,:), sca_y(:,:,:), sca_z(:,:,:), -delta_T, delta_T)
                endif


            ! end select
            endif

        end if

    end subroutine load_initial_state

    subroutine allocate_data()

        use decomp_2d

        use SCALAR_data
        implicit none


        allocate(sca_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        sca_x=0.d0

        allocate(sca_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
        sca_y=0.d0

        allocate(sca_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
        sca_z=0.d0

        ! Wall values allocation ---------------------------------------------
        allocate(sca_wall10(xstart(2):xend(2), xstart(3):xend(3)))
        sca_wall10=0.d0
        allocate(sca_wall11(xstart(2):xend(2), xstart(3):xend(3)))
        sca_wall11=0.d0

        allocate(sca_wall20(ystart(1):yend(1), ystart(3):yend(3)))
        sca_wall20=0.d0
        allocate(sca_wall21(ystart(1):yend(1), ystart(3):yend(3)))
        sca_wall21=0.d0

        allocate(sca_wall30(zstart(1):zend(1), zstart(2):zend(2)))
        sca_wall30=0.d0
        allocate(sca_wall31(zstart(1):zend(1), zstart(2):zend(2)))
        sca_wall31=0.d0

    end subroutine allocate_data

    subroutine initialize()

        use SCALAR_bc_controller
        use SCALAR_solver, only: SOLVER_init => init
        use SCALAR_inout_flow, only: OPEN_init=>inoutflow_init

        implicit none

        call allocate_data

        call set_numerical_boundaries
        call init_workspace
        call SOLVER_init
        call OPEN_init

!        call set_numerical_boundaries
!        call init_workspace
!        call SOLVER_init


        call load_initial_state

        contains

            subroutine init_workspace()

                use mpi


                use file_copy

                use COMMON_workspace_view
                use SCALAR_workspace_view

                implicit none


                results_path            =trim(COMMON_results_path)
                results3D_path          =trim(results_path)//"3D/"

                sensor_file             =trim(COMMON_log_path)//'Scalar/AA_v.csv'
                recovery_RHS_dir       =trim(COMMON_recovery_path)//"Scalar/previousRHS"
                recovery_fields_dir  =trim(COMMON_recovery_path)//"Scalar/Fields"
                recovery_sensor_file    =trim(COMMON_recovery_path)//'Scalar/AA_v.csv'

            end subroutine init_workspace

    end subroutine initialize

end module SCALAR_init

module SCALAR_recovery
    implicit none

contains

    subroutine save_state()


        use MPI
        use decomp_2d

        use SCALAR_workspace_view
        use SCALAR_dao
        use SCALAR_solver, SOVLER_save_state=>save_state

        implicit none

        call write_fields(recovery_fields_dir)
        call SOVLER_save_state

    end subroutine save_state

end module SCALAR_recovery


module SCALAR_final
    implicit none

contains

    subroutine deallocate_data()

        use SCALAR_data
        use SCALAR_solver, only: previousRHS1, previousRHS2, conv4o, diffo

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

        if (allocated(diffo)) deallocate(diffo)
        if (allocated(conv4o)) deallocate(conv4o)

    end subroutine deallocate_data

    subroutine finalize()
        implicit none

        call deallocate_data

    end subroutine finalize

end module SCALAR_final
