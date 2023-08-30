module embedded_init

    use embedded_data

    implicit none

contains

    subroutine read_settings()

        use decomp_2d

        use embedded_settings

        implicit none

        call read_start_settings

        call read_embedded_domain_settings

        call read_fringe_settings

        call read_sca_settings

    end subroutine read_settings

    subroutine update_allocation()

        use embedded_data

        implicit none

        xstart_e    = decomp_embedded%xst
        xend_e      = decomp_embedded%xen
        xsize_e     = decomp_embedded%xsz

        ystart_e    = decomp_embedded%yst
        yend_e      = decomp_embedded%yen
        ysize_e     = decomp_embedded%ysz

        zstart_e    = decomp_embedded%zst
        zend_e      = decomp_embedded%zen
        zsize_e     = decomp_embedded%zsz

    end subroutine update_allocation

end module

module embedded_initializer

    use embedded_init

    implicit none

    contains

    subroutine initialize()

        use mpi
        use decomp_2d

        use embedded_mesh
        use embedded_mesh_generator
        use embedded_data
        use embedded_fringe_data, only: use_fringe

        use embedded_scalar_data, only: SCA_state

        use embedded_start_settings, only: wanted_delta

        use embedded_velocity_init, only:    &
        embedded_velocity_initialize              => initialize

        use embedded_fringe_init, only:    &
        embedded_fringe_initialize              => initialize

        use embedded_scalar_init, only:    &
        embedded_scalar_initialize              => initialize

        implicit none
        integer         :: mpi_err

        call MPI_BARRIER(MPI_COMM_WORLD , mpi_err)

        call init_workspace

        call read_settings

        if (use_embedded) then

            if (nrank.eq.0) write(*,*) n1, n2, n3 
            call decomp_info_init(n1,n2,n3,decomp_embedded)

            call generate_mesh(n1, n2, n3, L1, L2, L3, stretch_Y)

            call update_allocation

            call embedded_velocity_initialize
            if (SCA_state/=0) call embedded_scalar_initialize
            if (use_fringe) call embedded_fringe_initialize

            if (nrank.eq.0) then

                write(*,*) 'The wanted BL thickness is', wanted_delta
                ! write(*,*) 'that corresponds to x =', Zc(index_for_output)

            endif

        endif

        contains

            subroutine init_workspace()

                use COMMON_workspace_view, only: COMMON_results_path
                use embedded_common_workspace_view

                implicit none

                ! Setting workspace paths ----------------------------------------
                embedded_common_results3D_path   = trim(COMMON_results_path)//"Embedded/"

                embedded_common_snapshot_path    = trim(COMMON_results_path)//"Snapshots/"


            end subroutine init_workspace

    end subroutine

end module embedded_initializer