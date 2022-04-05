module CORE_diagnostic

    use COMMON_workspace_view
    use run_ctxt_data
    use physical_fields

    use IVELOCITY_workspace_view

    implicit none

contains

    subroutine diagnostic_fields(error)

        use mesh
        use VELOCITY_quality_criteriums
        use VELOCITY_properties
        use IVELOCITY_IO_log_writer
        use snapshot_writer

        implicit none

        logical, intent(out)        :: error

        real*8      :: divg1, divg2
        real*8    ::div_max, div_diff, div_mean, cflm, vmax(3)
        integer   :: error_flag, l
        integer, parameter    :: NO_ERROR=0, MAX_DIV_REACHED=1, MAX_VELOCITY_IO_REACHED=2
        real*8, parameter    :: MAX_VEL=10000.d0, MAX_DIV=2.d-3

        ! Log
        integer :: master=0
        integer :: kinetic_file_id=201

        ! MPI
        real*8  :: glob_value, enstrophy_glob

        integer :: i,j,k

            ! The calculation stop if the velocities are diverging for numer
            ! stab conditions (courant number restrictions)

        error=.false.

        call perform_stability(q3_z, q2_y, q1_x, cflm, div_max, div_diff, div_mean)
        call perform_global_divergence(divg1, divg2, q3_z, q2_y, q1_x)
        if (nrank==master) call VELOCITY_IO_write_divergence(VELOCITY_divergence_history_file, cflm, div_max, div_mean, div_diff, divg1, divg2, ntime)


        error_flag=NO_ERROR

        call perform_maxvel(vmax)


        if(MAXVAL(vmax)>MAX_VEL)    error_flag=MAX_VELOCITY_IO_REACHED
        ! if(div_max>MAX_DIV)         error_flag=MAX_DIV_REACHED

        select case (error_flag)

            case (MAX_DIV_REACHED)
                if(nrank==0) write(*,*)'calculation diverged for  dmax=', div_max

                call create_snapshot(COMMON_snapshot_path, "DIVERGENCE", q1_y, "W", 2)
                call create_snapshot(COMMON_snapshot_path, "DIVERGENCE", q2_y, "V", 2)
                call create_snapshot(COMMON_snapshot_path, "DIVERGENCE", q3_y, "U", 2)

            case (MAX_VELOCITY_IO_REACHED)
                if(nrank==0) write(*,*)'*** calculation diverged in time for the vel field ***'
                if(nrank==0) write(*,*)'  vmax(1,2,3 ',(vmax(l),l=1,3)

                call create_snapshot(COMMON_snapshot_path, "DIVERGENCE", q1_y, "W", 2)
                call create_snapshot(COMMON_snapshot_path, "DIVERGENCE", q2_y, "V", 2)
                call create_snapshot(COMMON_snapshot_path, "DIVERGENCE", q3_y, "U", 2)

            case default

        end select

        if (error_flag/=NO_ERROR) error=.true.

    end subroutine diagnostic_fields

end module CORE_diagnostic
