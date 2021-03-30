module FRINGE_settings
    implicit none

contains

    subroutine read_settings()

        use COMMON_workspace_view
        use FRINGE_data
        use DNS_settings
        use mesh

        implicit none
        integer             :: fringe_state, n_interest_region

        open(15,file=trim(COMMON_settings_path)//'fringe.d')

        read(15,*) fringe_state  ! 0: no fringe ; 1: fringe activated
        read(15,*) fringe_length
        read(15,*) delta_rise
        read(15,*) delta_fall
        read(15,*) delta_activation
        read(15,*) max_strength_damping
        read(15,*)
        read(15,*) inflow_type ! 0 : poiseuille inflow, 1 : square inflow, 2 : inflow from file
        read(15,*) 
        read(15,*) 
        read(15,*)

        close(15)

        use_fringe=(fringe_state==1)

        if (use_fringe) then

            if (streamwise==1) then

                L1 = L1*(1+fringe_length)
                n_interest_region = int(n1/(1+fringe_length))
                n_fringe_region = n1 - n_interest_region
                n_fringe_end = n1
            endif

            if (streamwise==3) then

                L3 = L3*(1+fringe_length)
                n_interest_region = int(n3/(1+fringe_length))
                n_fringe_region = n3 - n_interest_region
                n_fringe_end = n3
            endif
    
            ! update the number of cells
            n1m=n1-1                !number of spanwise cells
            n2m=n2-1                !number of normal cells
            n3m=n3-1                !number of streamwise cells

            ! At this stage, L1/L3/n1 and n3 have not been updated
            n_delta_rise = int(delta_rise*n_fringe_region)
            n_delta_fall = int(delta_fall*n_fringe_region)
            n_fringe_start = n_interest_region + delta_activation
            n_delta_fall = 1

        endif


    end subroutine read_settings

end module FRINGE_settings

module FRINGE_dao

    use mpi
    use decomp_2d
    use decomp_2d_io
    use FRINGE_data

    implicit none

contains

    subroutine fringe_infos()
        implicit none
        if (nrank==0) then
            write(*,*) 'FRINGE INFOS___________________________'
            write(*,*)'fringe_length:',fringe_length
            write(*,*)'n_interest_region:',n_fringe_start-delta_activation, 'n_fringe_region:',n_fringe_region
            write(*,*)'n_fringe_start:', n_fringe_start, 'n_fringe_end:', n_fringe_end
            write(*,*)'n_delta_rise', n_delta_rise
            write(*,*)'n_delta_fall', n_delta_fall
            write(*,*) '_______________________________________________________'
        endif
    end subroutine fringe_infos

    subroutine write_fields(fields_dir)
!
        use physical_fields
        use mesh
        use DNS_settings
        use HDF5_IO
!
        implicit none
        character(*)    :: fields_dir
!
!
        character(200)    :: file_path
!
!
        if (streamwise==1) then
            file_path=trim(fields_dir)//"/f3_fringe"
            call hdf_write_3Dfield(file_path, f3_fringe_x(:,:,:), "f3_fringe", nx_global, ny_global, nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
    !
            file_path=trim(fields_dir)//"/f2_fringe"
            call hdf_write_3Dfield(file_path, f2_fringe_x(:,:,:), "f2_fringe", nx_global, ny_global, nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
    !
            file_path=trim(fields_dir)//"/f1_fringe"
            call hdf_write_3Dfield(file_path, f1_fringe_x(:,:,:), "f1_fringe", nx_global, ny_global, nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
        endif
!
        if (streamwise==3) then
            file_path=trim(fields_dir)//"/f3_fringe"
            call hdf_write_3Dfield(file_path, f3_fringe_z(:,:,:), "f3_fringe", nx_global, ny_global, nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
    !
            file_path=trim(fields_dir)//"/f2_fringe"
            call hdf_write_3Dfield(file_path, f2_fringe_z(:,:,:), "f2_fringe", nx_global, ny_global, nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
    !
            file_path=trim(fields_dir)//"/f1_fringe"
            call hdf_write_3Dfield(file_path, f1_fringe_z(:,:,:), "f1_fringe", nx_global, ny_global, nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
        endif

    end subroutine write_fields

end module FRINGE_dao



module FRINGE_results_writer

    use mpi
    use decomp_2d
    use decomp_2d_io
    use FRINGE_data

    implicit none

contains

    ! Export all physical field (velocity and scalar) in HDF5 format in individual files
    ! More generic, better than 2decomp format (not generic)
    subroutine write_fields(it)
!
        use run_ctxt_data
        use COMMON_workspace_view
        use FRINGE_dao, only: DAO_write_fields=>write_fields
!
        implicit none
!
        integer, intent(in) :: it
!
        character*10 tmp_str
        character*80 result_dir_path
!
        write(tmp_str, "(i10)")it
        result_dir_path=trim(COMMON_results3D_path)//'field'//trim(adjustl(tmp_str))
!
        call DAO_write_fields(result_dir_path)
!
    end subroutine write_fields

end module FRINGE_results_writer
