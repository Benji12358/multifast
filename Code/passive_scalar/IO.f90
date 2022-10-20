module SCALAR_settings
    implicit none

contains

    subroutine read_settings()

        use COMMON_workspace_view
        use DNS_settings
        use boundaries

        use SCALAR_data

        implicit none
        integer :: reset_int


        ! Scalar *********************************************************************************
        open(15,file=trim(COMMON_settings_path)//'scalar.d')

        read(15,*) SCA_state
        read(15,*) reset_int
        read(15,*) init_type ! CLASSIC_INIT=0, INIT_FROM_FILE=1
        read(15,*) SCA_BC1, SCA_BC2, SCA_BC3
        read(15,*) prandtl
        read(15,*) delta_T
        read(15,*) heat_flux
        heat_flux = heat_flux*h_height
        reset_scalar_field = (reset_int==1)

        close(15)



    end subroutine read_settings


end module SCALAR_settings

module SCALAR_dao

    use mpi
    use decomp_2d
    use decomp_2d_io
    use SCALAR_data

    implicit none

contains


    ! Read all physical field (velocity and scalar) in HDF5 format in individual files
    subroutine read_fields(fields_dir, sca_x_XYZ, decomp_XYZ, files_exist)

        use HDF5_IO
        use physical_fields
        implicit none
        character(*)   :: fields_dir
        character(200)    :: file_path
        logical, optional   :: files_exist(4)
        logical             :: fexist

        type(DECOMP_INFO)   :: decomp_XYZ
        real(mytype), dimension(decomp_XYZ%xst(1):decomp_XYZ%xen(1), decomp_XYZ%xst(2):decomp_XYZ%xen(2),     &
        decomp_XYZ%xst(3):decomp_XYZ%xen(3))   :: sca_x_XYZ

        file_path=trim(fields_dir)//"/sca1"

        inquire( file=trim(file_path)//".h5", exist=fexist)

        if (fexist) then
            call hdf_read_3Dfield(file_path, sca_x_XYZ, "sca1", decomp_XYZ%xsz(1),decomp_XYZ%ysz(2),decomp_XYZ%zsz(3), decomp_XYZ%xst(1),decomp_XYZ%xen(1),  decomp_XYZ%xst(2),decomp_XYZ%xen(2), decomp_XYZ%xst(3),decomp_XYZ%xen(3))

            ! file_path=trim(fields_dir)//"/ScaWall"

            ! call hdf_read_2Dfield(file_path, sca_wall10(:,:), "Wall10/sca", decomp_XYZ%yen(2), decomp_XYZ%zen(3), decomp_XYZ%xst(2),decomp_XYZ%xen(2), decomp_XYZ%xst(3),decomp_XYZ%xen(3))
            ! call hdf_read_2Dfield(file_path, sca_wall11(:,:), "Wall11/sca", decomp_XYZ%yen(2), decomp_XYZ%zen(3), decomp_XYZ%xst(2),decomp_XYZ%xen(2), decomp_XYZ%xst(3),decomp_XYZ%xen(3))

            ! call hdf_read_2Dfield(file_path, sca_wall20(:,:),"Wall20/sca", decomp_XYZ%xen(1), decomp_XYZ%zen(3), decomp_XYZ%yst(1),decomp_XYZ%yen(1), decomp_XYZ%yst(3),decomp_XYZ%yen(3))
            ! call hdf_read_2Dfield(file_path, sca_wall21(:,:),"Wall21/sca", decomp_XYZ%xen(1), decomp_XYZ%zen(3), decomp_XYZ%yst(1),decomp_XYZ%yen(1), decomp_XYZ%yst(3),decomp_XYZ%yen(3))

            ! call hdf_read_2Dfield(file_path, sca_wall30(:,:),"Wall30/sca", decomp_XYZ%xen(1), decomp_XYZ%yen(2), decomp_XYZ%zst(1),decomp_XYZ%zen(1), decomp_XYZ%zst(2),decomp_XYZ%zen(2))
            ! call hdf_read_2Dfield(file_path, sca_wall31(:,:),"Wall31/sca", decomp_XYZ%xen(1), decomp_XYZ%yen(2), decomp_XYZ%zst(1),decomp_XYZ%zen(1), decomp_XYZ%zst(2),decomp_XYZ%zen(2))

        endif

        if (present(files_exist)) then
            files_exist(4)=fexist
        endif


    end subroutine read_fields

    subroutine write_fields(fields_dir)

        use physical_fields
        use mesh
        use DNS_settings
        use HDF5_IO

        implicit none
        character(*)    :: fields_dir


        character(200)    :: file_path


        file_path=trim(fields_dir)//"/sca1"
        call hdf_write_3Dfield(file_path, sca_x(:,:,:), "sca1", nx_global, ny_global, nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))

        ! ! Save boundary conditions
        ! file_path=trim(fields_dir)//"/ScaWall"
        ! if(nrank==0)  call hdf_create_emptyfile(file_path)
        ! if(nrank==0)  call hdf_addgroup(file_path, "Wall10")
        ! call hdf_add_2Dfield(file_path, sca_wall10(:,:), "Wall10/sca", ny_global, nz_global, xstart(2),xend(2), xstart(3),xend(3))
        ! if(nrank==0)  call hdf_addgroup(file_path, "Wall11")
        ! call hdf_add_2Dfield(file_path, sca_wall11(:,:), "Wall11/sca", ny_global, nz_global, xstart(2),xend(2), xstart(3),xend(3))

        ! if(nrank==0)  call hdf_addgroup(file_path, "Wall20")
        ! call hdf_add_2Dfield(file_path, sca_wall20(:,:), "Wall20/sca", nx_global, nz_global, ystart(1),yend(1), ystart(3),yend(3))
        ! if(nrank==0)  call hdf_addgroup(file_path, "Wall21")
        ! call hdf_add_2Dfield(file_path, sca_wall21(:,:), "Wall21/sca", nx_global, nz_global, ystart(1),yend(1), ystart(3),yend(3))

        ! if(nrank==0)  call hdf_addgroup(file_path, "Wall30")
        ! call hdf_add_2Dfield(file_path, sca_wall30(:,:), "Wall30/sca", nx_global, ny_global, zstart(1),zend(1), zstart(2),zend(2))
        ! if(nrank==0)  call hdf_addgroup(file_path, "Wall31")
        ! call hdf_add_2Dfield(file_path, sca_wall31(:,:), "Wall31/sca", nx_global, ny_global, zstart(1),zend(1), zstart(2),zend(2))

    end subroutine write_fields

end module SCALAR_dao

module SCALAR_loader

    use mpi
    use decomp_2d
    use decomp_2d_io
    use SCALAR_data

    implicit none

contains


    subroutine load_fields(file, fill_from_coarse, fexist)

        use SCALAR_data
        use mesh
        use DNS_settings
        use start_settings, only:half_length, previous_fringe_start

        use SCALAR_dao, only: DAO_read_fields=>read_fields

        implicit none

        character(*), intent(in)        :: file
        logical, intent(in)             :: fill_from_coarse
        logical, optional               :: fexist(4)

        type(DECOMP_INFO)   :: decomp_XYZ
        integer             :: n1c, n2c, n3c        ! Coarse mesh resolution
        real*8, dimension(:,:,:), allocatable   :: sca_x_XYZ
        integer             :: i


        if (fill_from_coarse) then
            ! Coarse field dimensions :
!            n1c=(n1+1)/2
!            n2c=(n2+1)/2
!            n3c=(n3+1)/2
!
!            call decomp_info_init(n1c, n2c, n3c, decomp_XYZ)
!
!            allocate(sca_x_XYZ(decomp_XYZ%xst(1):decomp_XYZ%xen(1), decomp_XYZ%xst(2):decomp_XYZ%xen(2), decomp_XYZ%xst(3):decomp_XYZ%xen(3)))

            ! IMPLEMENT !!!!!

        else


            if (half_length.eq.1) then

                call decomp_info_init(previous_fringe_start, n2, n3, decomp_XYZ)
                call DAO_read_fields(file, sca_x(:previous_fringe_start,:,:), decomp_XYZ, fexist)

                do i=previous_fringe_start,n1
                    sca_x(i,:,:) = sca_x(previous_fringe_start,:,:)
                enddo

            else 

                call decomp_info_init(n1, n2, n3, decomp_XYZ)
                call DAO_read_fields(file, sca_x(:,:,:), decomp_XYZ, fexist)

            endif
            ! Spread to all transpositions

            call transpose_x_to_y(sca_x(:,:,:), sca_y(:,:,:))
            call transpose_y_to_z(sca_y(:,:,:), sca_z(:,:,:))

        end if

    end subroutine load_fields

    subroutine init_wall_sca()
        use SCALAR_data

        implicit none

        if(nrank==0) write(*,*)'SCALAR_generate_fields at wall'

        sca_wall10=0.d0
        sca_wall11=0.d0

        if (init_type==KAWAMURA_INIT) then
            sca_wall20=0.d0
            sca_wall21=0.d0
        else if (init_type==CONSTANT_HEAT_FLUX) then
            sca_wall20=-delta_T
            sca_wall21=-delta_T
        else
            sca_wall20=-delta_T
            sca_wall21=delta_T
        endif

        sca_wall30=0.d0
        sca_wall31=0.d0

    end subroutine init_wall_sca

end module SCALAR_loader

module SCALAR_results_writer

    use mpi
    use decomp_2d
    use decomp_2d_io
    use SCALAR_data

    implicit none

contains

    ! Export all physical field (velocity and scalar) in HDF5 format in individual files
    ! More generic, better than 2decomp format (not generic)
    subroutine write_fields(it)

        use run_ctxt_data
        use ISCALAR_workspace_view
        use SCALAR_dao, only: DAO_write_fields=>write_fields

        implicit none

        integer, intent(in) :: it

        character*10 tmp_str
        character*80 result_dir_path

        write(tmp_str, "(i10)")it
        result_dir_path=trim(SCALAR_results3D_path)//'field'//trim(adjustl(tmp_str))
        call DAO_write_fields(result_dir_path)

    end subroutine write_fields

end module SCALAR_results_writer
