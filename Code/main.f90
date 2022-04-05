module CORE_data
    implicit none
    integer     :: prow, pcol
    integer     ::nb_iteration, save_3D_frequency, save_recovery_frequency
    integer     :: nb_iteration_remaining
    real*8      :: t0, t1, t2
    real*8      :: max_runtime
end module CORE_data


module CORE_init

    use CORE_data

    implicit none

contains

    subroutine read_settings()

        use decomp_2d

        use mesh
        use numerical_methods_settings
        use boundaries
        use start_settings, only: n1c, n2c, n3c, start_from_coarse_file
        use fringe_data

        use COMMON_workspace_view

        use ISCALAR_IO_settings
        use IMHD_IO_settings
        use IFRINGE_IO_settings
        use IVELOCITY_IO_settings

        use CORE_IO_settings

        implicit none

        integer n,k
        integer start_from_coarse_file_int, IBM_activated_int, generic_poisson_flag
        integer :: s

        character*20, parameter :: NO_FILE="NONE", DEFAULT_PATH="DEFAULT_PATH"


        call read_start_settings
        call read_discretization_settings
        call read_domain_settings
        call read_global_settings(nb_iteration, save_3D_frequency, save_recovery_frequency)
        call read_IBM_settings
        call read_animation_settings

        call SCALAR_IO_read_settings
        call MHD_IO_read_settings
        call FRINGE_IO_read_settings
        call VELOCITY_IO_read_blowing_settings
        call VELOCITY_IO_read_twave_settings
        call read_bubble_settings

        call read_vortices_settings
        call read_counterrotating_vortices_settings

        !! nrank not yet defined, equals 0 for every process. Will be in decomp_2d_init
        if (start_from_coarse_file) then
            if( ((n1/=n1c).and.(n1/=2*n1c-1)).or. &
                ((n2/=n2c).and.(n2/=2*n2c-1)).or. &
                ((n3/=n3c).and.(n3/=2*n3c-1))) then

                if (nrank==0) write(*,*)'ERROR: the interpolation from a coarse field needs for each direction:'
                if (nrank==0) write(*,*)'Nfine = 2*Ncoarse - 1    or    Nfine = Ncoarse'
                call exit
            endif
        endif

        if (use_generic_poisson.and.(mesh_type==ORLANDI_MESH)) then
            if (nrank==0) write(*,*)'ERROR: GENERIC POISSON AND ORLANDI MESH ARE NOT COMPATIBLE!'
            if (nrank==0) write(*,*)'Choose Lamballais mesh or choose the Poisson physical solver'
            call exit
        endif

        if ((.not. use_generic_poisson).and.((BC1/=UNBOUNDED).or.(BC2/=NOSLIP).or.(BC3/=UNBOUNDED)).and.((BC1/=FRINGE))) then
            if (nrank==0) write(*,*)'ERROR: The physial solver is only suited for (BC1-BC2-BC3)=0-2-0 or with the use of the Fringe function'
            if (nrank==0) write(*,*)'Choose the Lamballais purely spectral solver or change the set of boundary conditions'
            call exit
        endif





    end subroutine read_settings

    subroutine load_initial_state()

        use DNS_settings
        use run_ctxt_data
        use start_settings

        use time_writer

        use COMMON_workspace_view

        implicit none
        integer             :: recovery_run_ctxt_id=22
        character*10 tmp_str
        character*200       :: timefile_path
        character*20        :: timefile_name

        ! Getting run context
        open(recovery_run_ctxt_id, file=trim(COMMON_recovery_run_ctxt_file))
        read(recovery_run_ctxt_id,*)run_ctxt

        if (run_ctxt==RECOVERY_A_RUN) then
            read(recovery_run_ctxt_id,*)nb_iteration
        else
            read(recovery_run_ctxt_id,*)
        end if

        if ((run_ctxt==RECOVERY_A_RUN).or.(run_ctxt==CONTINUE_FROM_PREVIOUS_RUN)) then
            read(recovery_run_ctxt_id,*)first_it
            read(recovery_run_ctxt_id,*)t
        endif

        close(recovery_run_ctxt_id)


        ! If the simulation start from previous fields, time information are read in the external simulation time file
        if ((run_ctxt==NEW_SIMULATION).and.(start_source_type==HDF5_FILE)) then

            write(tmp_str, "(i10)")start_it
            timefile_name="time"//trim(adjustl(tmp_str))
            timefile_path=trim(external_fields_path)//"/"//trim(timefile_name)


            call read_timefile(timefile_path, first_it, t)

        end if

        tinit=t

    end subroutine load_initial_state

    subroutine allocate_data()

        use decomp_2d
        use DNS_settings

        use IBM_data
        use IBM_settings
        use mesh

        use bubble_parallel_data
        use subdomains_view

        implicit none

        if (IBM_activated) then
            allocate(vel_term1(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
            vel_term1=0.d0
            allocate(vel_term2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
            vel_term2=0.d0
            allocate(vel_term3(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
            vel_term3=0.d0
        endif

        allocate(IBM_mask1_x(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
        IBM_mask1_x=0.d0
        allocate(flow_rate_IBM(n1))
        flow_rate_IBM=0.d0
        allocate(kin1_IBM(n1))
        kin1_IBM=0.d0
        allocate(kin2_IBM(n1))
        kin2_IBM=0.d0
        allocate(kin3_IBM(n1))
        kin3_IBM=0.d0
        allocate(q1max_IBM(n1))
        q1max_IBM=0.d0
        allocate(q2max_IBM(n1))
        q2max_IBM=0.d0
        allocate(q3max_IBM(n1))
        q3max_IBM=0.d0
        allocate(prmax_IBM(n1))
        prmax_IBM=0.d0
        allocate(q1min_IBM(n1))
        q1min_IBM=0.d0
        allocate(q2min_IBM(n1))
        q2min_IBM=0.d0
        allocate(q3min_IBM(n1))
        q3min_IBM=0.d0
        allocate(prmin_IBM(n1))
        prmin_IBM=0.d0

        ! IBM_mask1_x is used even when IBM is not activated and is always allocated
        if (IBM_activated) then
            allocate(IBM_mask2_x(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
            IBM_mask2_x=0.d0
            allocate(IBM_mask3_x(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
            IBM_mask3_x=0.d0

            if (interpol_type==ANTISYMMETRIC_INTERPOL) then

                allocate(IBM_modulation_x(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
                IBM_modulation_x = 0.d0
                allocate(IBM_modulation_y(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
                IBM_modulation_y = 0.d0
                allocate(IBM_modulation_z(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
                IBM_modulation_z = 0.d0
            endif

        end if


        ! Bubble_parallel_data
        ! A process can post all of its bubble (in case of exporting bubble for example)
!        allocate(send_list(0:nproc-1, MAX_BUBBLES, BUBBLE_SIZE))
!        allocate(send_nb(0:nproc-1))
!        send_nb=0
!        allocate(receive_nb(0:nproc-1))
!        receive_nb=0

        ! MPI global decomp data
        allocate(xpos_start_glob(0:nproc-1, 3))
        xpos_start_glob=0.d0
        allocate(ypos_start_glob(0:nproc-1, 3))
        ypos_start_glob=0.d0
        allocate(zpos_start_glob(0:nproc-1, 3))
        zpos_start_glob=0.d0

        allocate(xstart_glob(0:nproc-1, 3))
        xstart_glob=0
        allocate(ystart_glob(0:nproc-1, 3))
        ystart_glob=0
        allocate(zstart_glob(0:nproc-1, 3))
        zstart_glob=0

        allocate(xpos_end_glob(0:nproc-1, 3))
        xpos_end_glob=0.d0
        allocate(ypos_end_glob(0:nproc-1, 3))
        ypos_end_glob=0.d0
        allocate(zpos_end_glob(0:nproc-1, 3))
        zpos_end_glob=0.d0

        allocate(xend_glob(0:nproc-1, 3))
        xend_glob=0
        allocate(yend_glob(0:nproc-1, 3))
        yend_glob=0
        allocate(zend_glob(0:nproc-1, 3))
        zend_glob=0


    end subroutine allocate_data

    subroutine initialize()

        use mpi
        use decomp_2d

        use run_ctxt_data
        use mesh
        use numerical_methods_settings
        use time_schemes
        use inflow_settings
        use mesh_generator
        use subdomains_view
        use schemes_loader
        use flow_buffer_data
        use flow_buffer_handler, only: OPEN_BUFFER_init=>init

        use IVELOCITY_OPEN
        use ISCALAR_OPEN

        use CORE_IO_settings

        implicit none
        character(50)   :: progamm_param
        integer         :: mpi_err
        integer         :: testmode_flag

        ! Getting the command line parameters----------------------------------------------
        call getarg(2, progamm_param)
        read(progamm_param,*)max_runtime

        call getarg(3, progamm_param)
        read(progamm_param,*)prow

        call getarg(4, progamm_param)
        read(progamm_param,*)pcol

        call getarg(5, progamm_param)
        read(progamm_param,*)testmode_flag

        test_mode=(testmode_flag==1)

        call MPI_BARRIER(MPI_COMM_WORLD , mpi_err)

        call init_workspace

        call read_settings
        call decomp_2d_init(n1, n2, n3, prow, pcol)

        call resume_settings

        call allocate_data

!        call load_initial_state

        call generate_mesh(n1, n2, n3, L1, L2, L3, stretch_Y)
        call init_subdomains_view

        call configure_temporal_schemes(time_advancement, ga, ro, nb_substep)
        call configure_schemes(schemes_configuration)

        if (inout_newversion) then

            if(outflow_buff>0)  &
                call OPEN_BUFFER_init(WRITE_MODE, outflow_buff, COMMON_recovery_outflow_dir, COMMON_outflow_path)
            if(inflow_buff>0)  &
                call OPEN_BUFFER_init(READ_MODE, inflow_buff, COMMON_recovery_outflow_dir, COMMON_inflow_path)
        else


            call VELOCITY_OPEN_init(nb_substep, outflow_buff, inflow_buff)
            call SCALAR_OPEN_init(nb_substep, outflow_buff, inflow_buff)

        endif




        call load_initial_state

        contains

            subroutine init_workspace()

                use COMMON_workspace_view

                implicit none

                ! Create current simulation arborescence -------------------------
                call getarg(1, progamm_param)
                read(progamm_param,*)COMMON_simulation_name

                ! Setting workspace paths ----------------------------------------
                COMMON_simulation_path=trim(COMMON_simulation_name)//"/"

                COMMON_log_path                =trim(COMMON_simulation_path)//"Log/"

                COMMON_results_path            =trim(COMMON_simulation_path)//"Results/"
                COMMON_results3D_path          =trim(COMMON_results_path)//"3D/"
                COMMON_outflow_path            =trim(COMMON_results_path)//"Outflow/"

                COMMON_recovery_path           =trim(COMMON_simulation_path)//".recovery/"
                COMMON_recovery_outflow_dir    =trim(COMMON_recovery_path)//"Outflow/"
                COMMON_recovery_run_ctxt_file  =trim(COMMON_recovery_path)//"run_ctxt.d"


                COMMON_lapY_file           =trim(COMMON_simulation_path)//"Input/lapY_matrix.dat"


                COMMON_settings_path           =trim(COMMON_simulation_path)//"Input/Settings/"

                COMMON_snapshot_path=trim(COMMON_results_path)//"Snapshots/"
                COMMON_anim2D_path          =trim(COMMON_results_path)//"Animation/Anim2D/"
                COMMON_anim2D_bubble_path   =trim(COMMON_results_path)//"Animation/Bubbles/"



            end subroutine init_workspace

    end subroutine initialize

end module CORE_init

module CORE_recovery
    implicit none

contains

    subroutine load_state()
        implicit none

    end subroutine load_state

    subroutine save_state()

        use time_writer

        use COMMON_workspace_view
        use numerical_methods_settings, only: inout_newversion
        use flow_buffer_handler, only: OPEN_BUFFER_save_state=>save_state

        use IVELOCITY_LIFE
        use ISCALAR_LIFE

        use CORE_data
        use SCALAR_data, only : SCA_state
        use inflow_settings

        implicit none
        integer             :: recovery_run_ctxt_id=22
        character*10 tmp_str

                ! Replace the simulation-recovery file (velocity + RHS)

        if ((inout_newversion)  .and. (outflow_buff>0)) call OPEN_BUFFER_save_state
!         if (inout_newversion) call OPEN_BUFFER_save_state
        call VELOCITY_save_state

        write(tmp_str, "(i10)")ntime
        call write_timefile( trim(COMMON_recovery_path)//'time' )


        if (SCA_state/=0) call SCALAR_save_state


        ! Setting recovery information  -----------------------------------
        if (nrank==0) then
            open(recovery_run_ctxt_id, file=trim(COMMON_recovery_run_ctxt_file))
            write(recovery_run_ctxt_id,*)run_status
            write(recovery_run_ctxt_id,*)nb_iteration_remaining
            write(recovery_run_ctxt_id,*)ntime
            write(recovery_run_ctxt_id,*)t
            close(recovery_run_ctxt_id)

        end if

    end subroutine save_state

end module CORE_recovery


module CORE_final
    implicit none

contains

    subroutine CORE_deallocate_data()

        use IBM_data
        use IBM_settings
        use subdomains_view

        implicit none

        deallocate(IBM_mask1_x)
        deallocate(flow_rate_IBM)
        deallocate(kin1_IBM)
        deallocate(kin2_IBM)
        deallocate(kin3_IBM)
        deallocate(q1max_IBM)
        deallocate(q2max_IBM)
        deallocate(q3max_IBM)
        deallocate(prmax_IBM)
        deallocate(q1min_IBM)
        deallocate(q2min_IBM)
        deallocate(q3min_IBM)
        deallocate(prmin_IBM)

        ! IBM_mask1_x is used even when IBM is not activated and is always allocated
        if (IBM_activated) then
            deallocate(IBM_mask2_x)
            deallocate(IBM_mask3_x)

            deallocate(vel_term1)
            deallocate(vel_term2)
            deallocate(vel_term3)
        end if


        ! MPI global decomp data
        deallocate(xpos_start_glob)
        deallocate(ypos_start_glob)
        deallocate(zpos_start_glob)

        deallocate(xstart_glob)
        deallocate(ystart_glob)
        deallocate(zstart_glob)

        deallocate(xpos_end_glob)
        deallocate(ypos_end_glob)
        deallocate(zpos_end_glob)

        deallocate(xend_glob)
        deallocate(yend_glob)
        deallocate(zend_glob)

    end subroutine CORE_deallocate_data

    subroutine finalize()

        use mpi
        use decomp_2d

        use physical_fields
        use DNS_settings

        use VELOCITY_final, VELOCITY_finalize=>finalize
        use SCALAR_final, SCALAR_finalize=>finalize
        use MHD_final, MHD_finalize=>finalize
        use FRINGE_final, FRINGE_finalize=>finalize

        use CORE_data
        use MHD_data, only: MHD_state
        use SCALAR_data, only : SCA_state

        implicit none

        integer :: mpi_err
        real*8  :: time_per_it
        real*8  :: t12, t02
        real*8  :: t02_glob, t12_glob

        t2 = MPI_WTIME()

        t02=t2-t0
        t12=t2-t1
        call MPI_ALLREDUCE (t02, t02_glob, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (t12, t12_glob, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)

        time_per_it=t12_glob/((last_it-first_it+1)-nb_iteration_remaining)

        if(nrank==0) write(*,*)
        if(nrank==0) write(*,*)'Nb de procs:', nproc
        if(nrank==0) write(*,*)'Temps execution:', t02_glob
        if(nrank==0) write(*,*)'Temps boucle:', t12_glob
        if(nrank==0) write(*,*)'Temps par it:', time_per_it

        ! Close all modules
        call VELOCITY_finalize
        if (SCA_state/=0) call SCALAR_finalize
        if (MHD_state/=0) call MHD_finalize

        call CORE_deallocate_data

    end subroutine finalize

end module CORE_final



program main

    use mpi

    use run_ctxt_data
    use time_schemes
    use time_writer
    use mesh
    use anim2D
    use DNS_settings
    use inflow_settings
    use numerical_methods_settings, only: inout_newversion
    use flow_buffer_handler, only: OPEN_BUFFER_write_cursorfile=>write_cursorfile

    use IVELOCITY_IO
    use IVELOCITY_LIFE
    use ISCALAR_IO
    use ISCALAR_LIFE

    use IMHD_IO
    use IMHD_LIFE

    use IFRINGE_IO
    use IFRINGE_LIFE

    use IBM
    use CORE_IO_settings
    use CORE_data
    use CORE_init, CORE_initialize=>initialize
    use CORE_recovery
    use CORE_final
    use CORE_diagnostic
    use CORE_log_writers
    use anim2D_writer


    use multiphysics
    use IVELOCITY_OPEN
    use ISCALAR_OPEN

    use Bubble_generator

    use MHD_data, only: MHD_state, MHD_export_3D
    use FRINGE_data, only: use_fringe

    use flow_buffer_handler_test

    implicit none

    integer ::i
    integer :: mpi_err

    character*10 tmp_str
    character*80 time_filename

    logical :: error

! ***********************************************************************************
! ******************************** Initialization ***********************************
! ***********************************************************************************

    call MPI_INIT(mpi_err)

    t0 = MPI_WTIME()

    call CORE_initialize

    if (IBM_activated) then
        call IBM_setup
    else
        IBM_mask1_x=0.d0
        ibm_volume=0.d0
    endif
    
    call VELOCITY_initialize
    if (SCA_state/=0) call SCALAR_initialize
    if (MHD_state/=0) call MHD_initialize
    if (use_fringe) call FRINGE_initialize


    first_it=first_it+1
    last_it=first_it+nb_iteration-1
    run_status=IN_PROGRESS

    nb_iteration_remaining=nb_iteration


    if ((inflow_buff>0).and.(.not. inout_newversion)) then
    !Synchronizes 3D field with the last substep: we have to load/shift
    !the boundaries conditions associated to the 3D field
        do i = 1, nb_substep
            call VELOCITY_OPEN_get_inflow(first_it,i)
            if (SCA_state==1) call SCALAR_OPEN_get_inflow(first_it,i)
        end do

    endif





! ***********************************************************************************
! ***************************   CORE OF THE DNS *************************************
! ***********************************************************************************

    t1 = MPI_WTIME()

!    call recovery_test
!    call inflow_outflow_communication_test

!    call exit(0)

    do ntime = first_it, last_it

        nb_iteration_remaining=nb_iteration_remaining-1

        call write_sensorfile(n1/2, n2/2, n3/2)

        !!! Solving the physics : Navier-Stokes, MHD, scalar transport... !!!!
        call perform_multiphysics_MHD(ntime)

        t=t+dt

        ! Check the validity of the physical fields ---------------------------------------------
        call diagnostic_fields(error)
        if(error) call finalize

        ! Write velocity properties (kinetic energy, max velocity etc) in Log directory-----------
        call VELOCITY_IO_write_properties

        ! Save pressure gradient in TMP directory
        if (mod(ntime, save_gradP_frequency).eq.0) call VELOCITY_IO_write_gradP
        if (mod(ntime, save_gradP_frequency).eq.0) call VELOCITY_IO_write_flowrate

        ! Save velocity in Results directory------------------------------------------------------
        if(mod(ntime, save_3D_frequency).eq.0) then

            write(tmp_str, "(i10)")ntime
            time_filename="time"//trim(adjustl(tmp_str))
            call write_timefile( trim(COMMON_results_path)//time_filename )

            if ((inout_newversion) .and. (outflow_buff>0)) call OPEN_BUFFER_write_cursorfile

            call VELOCITY_IO_write_fields(ntime)
            if (SCA_state==1) call SCALAR_IO_write_fields(ntime)
            if ((MHD_state/=0).and.(MHD_export_3D.eq.1)) call MHD_IO_write_fields(ntime)
            ! if (use_fringe) call FRINGE_IO_write_fields(ntime)
            ! call CELL_CENTERED_IO_WRITE_FIELDS(ntime,t)
        endif

        ! Animation2D -----------------------------------------------------------------------------
        if (anim2D_1/=0) then
            if(mod(ntime, param_anim2D_1%step_size).eq.0) call anim2D_addframe1(COMMON_anim2D_path)
        endif
        if (anim2D_3/=0) then
            if(mod(ntime, param_anim2D_3%step_size).eq.0) call anim2D_addframe3(COMMON_anim2D_path)
        endif

!        call anim2D_bubble(ntime)

        ! Save simulation state for recovery support ----------------------------------------------
        if(mod(ntime, save_recovery_frequency).eq.0) then
            call save_state
        endif

        if (test_mode) then
            open(15, file="current_it")
            write(15,*) ntime
            close(15)
        end if


        t2 = MPI_WTIME()

        if ((t2-t0)>max_runtime) then

            if(nrank==0) write(*,*) 'Run time exceeded... Save for recovery'
            call save_state
            call finalize

        end if

    enddo


    !   *******************************************************************************
    !   ***************************** End of simulation *******************************
    !   *******************************************************************************


    ntime=last_it
    run_status=FINISHED


    call save_state
    call finalize

    call MPI_FINALIZE(mpi_err)


end
