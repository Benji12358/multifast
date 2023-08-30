module following_init

    use following_data

    implicit none

contains

    subroutine read_settings()

        use decomp_2d

        use following_settings

        implicit none

        call read_start_settings

        call read_following_domain_settings

        call read_fringe_settings

        call read_sca_settings

    end subroutine read_settings

    subroutine update_allocation()

        use following_data

        implicit none

        xstart_f    = decomp_following%xst
        xend_f      = decomp_following%xen
        xsize_f     = decomp_following%xsz

        ystart_f    = decomp_following%yst
        yend_f      = decomp_following%yen
        ysize_f     = decomp_following%ysz

        zstart_f    = decomp_following%zst
        zend_f      = decomp_following%zen
        zsize_f     = decomp_following%zsz

    end subroutine update_allocation

end module

module following_initializer

    use following_init

    implicit none

    contains

    subroutine initialize()

        use mpi
        use decomp_2d

        use following_mesh
        use following_mesh_generator
        use following_data
        use following_fringe_data, only: use_fringe

        use following_scalar_data, only: SCA_state

        use following_start_settings, only: index_for_output

        use following_velocity_init, only:    &
        following_velocity_initialize              => initialize

        use following_fringe_init, only:    &
        following_fringe_initialize              => initialize

        use following_scalar_init, only:    &
        following_scalar_initialize              => initialize

        implicit none
        integer         :: mpi_err

        call MPI_BARRIER(MPI_COMM_WORLD , mpi_err)

        call init_workspace

        call read_settings

        if (use_following) then

            if (nrank.eq.0) write(*,*) n1, n2, n3 
            call decomp_info_init(n1,n2,n3,decomp_following)

            call generate_mesh(n1, n2, n3, L1, L2, L3, stretch_Y)

            call update_allocation

            call following_velocity_initialize
            if (SCA_state/=0) call following_scalar_initialize
            if (use_fringe) call following_fringe_initialize

        endif

        contains

            subroutine init_workspace()

                use COMMON_workspace_view, only: COMMON_results_path
                use following_common_workspace_view

                implicit none

                ! Setting workspace paths ----------------------------------------
                following_common_results3D_path   = trim(COMMON_results_path)//"following/"

                following_common_snapshot_path    = trim(COMMON_results_path)//"Snapshots/"


            end subroutine init_workspace

    end subroutine

end module following_initializer