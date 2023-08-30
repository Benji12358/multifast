module following_settings
    use COMMON_workspace_view, only: COMMON_settings_path 
    implicit none

contains

    subroutine read_start_settings()

        use following_start_settings
        use following_data

        implicit none

        integer             :: reading_format
        integer             :: use_following_int 

        ! Start **********************************************************************************
        open(15,file=trim(COMMON_settings_path)//'following/start.d')

        read(15,*) use_following_int
        read(15,*) start_source_type
        read(15,*) index_for_output
        read(15,*)
        read(15,'(a)') external_fields_path  ! path of the read field  file (cha.rea)
        read(15,*) start_it
        read(15,*) vper
        read(15,*) reading_format
        close(15)

        use_following = (use_following_int==1)


    end subroutine read_start_settings

    subroutine read_following_domain_settings()

        use mathematical_constants

        use following_mesh

        implicit none
        integer     :: pi_scale_flag

        ! Computational domain *******************************************************************
        open(15,file=trim(COMMON_settings_path)//'following/computational_domain.d')

        read(15,*) n1,n2,n3
        read(15,*) L3, L2, L1
        read(15,*) stretch_Y, mesh_type
        read(15,*) pi_scale_flag

        close(15)

        n1m=n1-1                !number of spanwise cells
        n2m=n2-1                !number of normal cells
        n3m=n3-1                !number of streamwise cells

        if (pi_scale_flag.eq.1) then
            L1=L1*pi        !spanwise size of the whole calcul box
            L3=L3*pi        !streamwise size of the whole calcul box
        else
            L1=L1           !spanwise size of the whole calcul box
            L3=L3           !streamwise size of the whole calcul box
        endif

    end subroutine read_following_domain_settings

    subroutine read_fringe_settings()

        use COMMON_workspace_view
        use following_fringe_data
        use following_mesh
        use following_data, only: u_bulk

        implicit none
        integer             :: fringe_state

        open(15,file=trim(COMMON_settings_path)//'following/fringe.d')

        read(15,*) fringe_state  ! 0: no fringe ; 1: fringe activated
        read(15,*) fringe_length
        read(15,*) delta_rise
        read(15,*) delta_fall
        read(15,*) delta_activation
        read(15,*) max_strength_damping
        read(15,*)
        read(15,*) inflow_type ! 0 : poiseuille inflow, 1 : square inflow, 2 : inflow from file
        read(15,*) u_bulk

        close(15)

        use_fringe=(fringe_state==1)

        if (use_fringe) then

            L1 = L1*(1+fringe_length)
            n_interest_region = int(n1/(1+fringe_length))
            n_fringe_region = n1 - n_interest_region - 1
            ! n_fringe_end = n1
    
            ! update the number of cells
            n1m=n1-1                !number of spanwise cells
            n2m=n2-1                !number of normal cells
            n3m=n3-1                !number of streamwise cells

            ! At this stage, L1/L3/n1 and n3 have not been updated
            n_delta_rise = int(delta_rise*n_fringe_region)
            n_delta_fall = int(delta_fall*n_fringe_region)
            
            ! with fringe at the start
            n_fringe_start = 1
            n_fringe_end = n_fringe_region !+ delta_activation
            ! n_interest_region_start = n_fringe_region + 1

        endif

        if (fringe_state==0) u_bulk=1 ! reset u_bulk value if fringe not used

    end subroutine read_fringe_settings

    subroutine read_sca_settings()

        use following_common_workspace_view
        use DNS_settings
        use boundaries

        use following_scalar_data

        implicit none
        integer :: reset_int

        ! Scalar *********************************************************************************
        open(15,file=trim(COMMON_settings_path)//'scalar.d')

        read(15,*) SCA_state
        read(15,*) reset_int
        read(15,*) init_type
        read(15,*) 
        read(15,*) prandtl
        read(15,*) delta_T
        read(15,*) heat_flux

        heat_flux = heat_flux*h_height
        reset_scalar_field = (reset_int==1)

        close(15)



    end subroutine read_sca_settings

end module following_settings

module following_snapshot_writer

    use mpi
    use decomp_2d
    use following_data

    implicit none

contains

    subroutine create_snapshot(snaps_dir, snap_dir, field, field_name, pencil)

        use following_mesh
        use HDF5_IO

        implicit none
        real*8, dimension(:,:,:)    :: field
        character(*)                :: snaps_dir, snap_dir, field_name
        integer                     :: pencil

        integer                     :: mpi_err


        character(200)    :: file_path, snap_path
        logical             :: snap_dir_exist


        snap_path=trim(snaps_dir)//"/"//trim(snap_dir)

        file_path=trim(snap_path)//"/"//trim(field_name)

        if(nrank==0)  call hdf_create_emptyfile(file_path)
        !if (pencil==1) call hdf_add_3Dfield(file_path, field(:,:,:), trim(field_name), n1, n2, n3, 1,xsize(1),1,xsize(2),1,xsize(3))
        if (pencil==1) call hdf_add_3Dfield(file_path, field(:,:,:), trim(field_name), n1, n2, n3, xstart_f(1),xend_f(1),xstart_f(2),xend_f(2),xstart_f(3),xend_f(3))
        !if (pencil==2) call hdf_add_3Dfield(file_path, field(:,:,:), trim(field_name), n1, n2, n3, 1,ysize(1),1,ysize(2),1,ysize(3))
        if (pencil==2) call hdf_add_3Dfield(file_path, field(:,:,:), trim(field_name), n1, n2, n3, ystart_f(1),yend_f(1),ystart_f(2),yend_f(2),ystart_f(3),yend_f(3))
        !if (pencil==3) call hdf_add_3Dfield(file_path, field(:,:,:), trim(field_name), n1, n2, n3, 1,zsize(1),1,zsize(2),1,zsize(3))
        if (pencil==3) call hdf_add_3Dfield(file_path, field(:,:,:), trim(field_name), n1, n2, n3, zstart_f(1),zend_f(1),zstart_f(2),zend_f(2),zstart_f(3),zend_f(3))

        call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

    end subroutine create_snapshot

end module following_snapshot_writer


! module FRINGE_results_writer

!     use mpi
!     use decomp_2d
!     use decomp_2d_io
!     use FRINGE_data

!     implicit none

! contains

!     ! Export all physical field (velocity and scalar) in HDF5 format in individual files
!     ! More generic, better than 2decomp format (not generic)
!     subroutine write_fields(it)
! !
!         use run_ctxt_data
!         use COMMON_workspace_view
!         use FRINGE_dao, only: DAO_write_fields=>write_fields
! !
!         implicit none
! !
!         integer, intent(in) :: it
! !
!         character*10 tmp_str
!         character*80 result_dir_path
! !
!         write(tmp_str, "(i10)")it
!         result_dir_path=trim(COMMON_results3D_path)//'field'//trim(adjustl(tmp_str))
! !
!         call DAO_write_fields(result_dir_path)
! !
!     end subroutine write_fields

! end module FRINGE_results_writer

module following_velocity_dao

    use mpi
    use decomp_2d
    use decomp_2d_io
    use following_data

    implicit none

contains

    ! Export all velocity fields in HDF5 format in individual files
    ! More generic, better than 2decomp format (not generic)
    subroutine write_fields(fields_dir)

        use following_physical_fields
        use following_mesh
        use DNS_settings
        use HDF5_IO
        use following_data

        implicit none
        character(*)    :: fields_dir


        character(200)    :: file_path

        file_path=trim(fields_dir)//"/U"
        call hdf_write_3Dfield(file_path, q3_z(:,:,:), "U", n1, n2, n3, zstart_f(1),zend_f(1),zstart_f(2),zend_f(2),zstart_f(3),zend_f(3))

        file_path=trim(fields_dir)//"/V"
        call hdf_write_3Dfield(file_path, q2_y(:,:,:), "V", n1, n2, n3, ystart_f(1),yend_f(1),ystart_f(2),yend_f(2),ystart_f(3),yend_f(3))

        file_path=trim(fields_dir)//"/W"
        call hdf_write_3Dfield(file_path, q1_x(:,:,:), "W", n1, n2, n3, xstart_f(1),xend_f(1),xstart_f(2),xend_f(2),xstart_f(3),xend_f(3))

        file_path=trim(fields_dir)//"/P"
        call hdf_write_3Dfield(file_path, dp_x(:,:,:), "P", n1, n2, n3, xstart_f(1),xend_f(1),xstart_f(2),xend_f(2),xstart_f(3),xend_f(3))

    end subroutine write_fields


    ! Read all velocity fields in HDF5 format in individual files
    subroutine read_fields(fields_dir, u_x_XYZ, v_x_XYZ, w_x_XYZ, dp_x_XYZ, decomp_XYZ, files_exist)

        use HDF5_IO
        use following_physical_fields
        implicit none
        character(*)   :: fields_dir
        character(200)    :: file_path
        logical, optional   :: files_exist(4)
        logical             :: fexist

        type(DECOMP_INFO)   :: decomp_XYZ
        real(mytype), dimension(decomp_XYZ%xst(1):decomp_XYZ%xen(1), decomp_XYZ%xst(2):decomp_XYZ%xen(2),     &
        decomp_XYZ%xst(3):decomp_XYZ%xen(3))   :: u_x_XYZ, v_x_XYZ, w_x_XYZ, dp_x_XYZ

        ! WARNING : NOT VALIDATED FOR COARSE FIELDS
        file_path=trim(fields_dir)//"/U"
        call hdf_read_3Dfield(file_path, u_x_XYZ, "U", decomp_XYZ%xsz(1),decomp_XYZ%ysz(2),decomp_XYZ%zsz(3), decomp_XYZ%xst(1),decomp_XYZ%xen(1),  decomp_XYZ%xst(2),decomp_XYZ%xen(2), decomp_XYZ%xst(3),decomp_XYZ%xen(3))

        file_path=trim(fields_dir)//"/V"
        call hdf_read_3Dfield(file_path, v_x_XYZ, "V", decomp_XYZ%xsz(1),decomp_XYZ%ysz(2),decomp_XYZ%zsz(3), decomp_XYZ%xst(1),decomp_XYZ%xen(1),  decomp_XYZ%xst(2),decomp_XYZ%xen(2), decomp_XYZ%xst(3),decomp_XYZ%xen(3))

        file_path=trim(fields_dir)//"/W"
        call hdf_read_3Dfield(file_path, w_x_XYZ, "W", decomp_XYZ%xsz(1),decomp_XYZ%ysz(2),decomp_XYZ%zsz(3), decomp_XYZ%xst(1),decomp_XYZ%xen(1),  decomp_XYZ%xst(2),decomp_XYZ%xen(2), decomp_XYZ%xst(3),decomp_XYZ%xen(3))

        file_path=trim(fields_dir)//"/P"
        call hdf_read_3Dfield(file_path, dp_x_XYZ, "P", decomp_XYZ%xsz(1),decomp_XYZ%ysz(2),decomp_XYZ%zsz(3), decomp_XYZ%xst(1),decomp_XYZ%xen(1),  decomp_XYZ%xst(2),decomp_XYZ%xen(2), decomp_XYZ%xst(3),decomp_XYZ%xen(3))

    end subroutine read_fields

end module following_velocity_dao


module following_velocity_results_writer

    use mpi
    use decomp_2d
    use decomp_2d_io
    use following_data

    implicit none
    integer, parameter  :: SEQUENTIAL_FILE=0, DECOMP2D_FILE=1
    integer, parameter  :: BINARY_FILE=0, ASCII_FILE=1
    integer, parameter  :: SIMPLE_FILE=0, DIRECTORY=1



contains

    ! Export all velocity fields in HDF5 format in individual files
    ! More generic, better than 2decomp format (not generic)
    subroutine write_fields(it)
        use time_writer
        use following_common_workspace_view
        use following_velocity_dao, only: DAO_write_fields=> write_fields

        implicit none

        integer, intent(in)         :: it

        character*10 tmp_str
        character*80 result_dir_path

        write(tmp_str, "(i10)")it
        result_dir_path=trim(following_common_results3D_path)//'field'//trim(adjustl(tmp_str))

        call write_timefile(trim(result_dir_path)//"/advancement.d")
        call DAO_write_fields(result_dir_path)

    end subroutine write_fields

end module following_velocity_results_writer

module following_velocity_loader
    use mpi
    use decomp_2d
    use decomp_2d_io
    use following_data

    implicit none

contains


    subroutine load_fields(file, fill_from_coarse, fexist)

        use following_physical_fields
        use following_mesh
        use DNS_settings
        use following_start_settings

        use following_velocity_dao, only: DAO_read_fields => read_fields

        implicit none

        character(*), intent(in)        :: file
        logical, intent(in)             :: fill_from_coarse
        logical, optional               :: fexist(4)

        type(DECOMP_INFO)   :: decomp_XYZ
!        integer, intent(in)             :: n1c, n2c, n3c        ! Coarse mesh resolution
        real*8, dimension(:,:,:), allocatable   :: u_x_XYZ, v_x_XYZ, w_x_XYZ, dp_x_XYZ
        integer             :: i

        call decomp_info_init(n1, n2, n3, decomp_XYZ)
        call DAO_read_fields(file, q3_x, q2_x, q1_x, dp_x, decomp_XYZ, fexist)

        ! Spread to all transpositions
        call transpose_x_to_y(q3_x, q3_y, decomp_following)
        call transpose_y_to_z(q3_y, q3_z, decomp_following)

        call transpose_x_to_y(q2_x, q2_y, decomp_following)
        call transpose_y_to_z(q2_y, q2_z, decomp_following)

        call transpose_x_to_y(q1_x, q1_y, decomp_following)
        call transpose_y_to_z(q1_y, q1_z, decomp_following)

        call transpose_x_to_y(dp_x, dp_y, decomp_following)
        call transpose_y_to_z(dp_y, dp_z, decomp_following)

    end subroutine load_fields

end module following_velocity_loader


module following_scalar_dao

    use mpi
    use decomp_2d
    use decomp_2d_io
    use following_scalar_data

    implicit none

contains

    ! Read all physical field (velocity and scalar) in HDF5 format in individual files
    subroutine read_fields(fields_dir, sca_x_XYZ, decomp_XYZ, files_exist)

        use HDF5_IO
        implicit none
        character(*)   :: fields_dir
        character(200)    :: file_path
        logical, optional   :: files_exist(4)
        logical             :: fexist

        type(DECOMP_INFO)   :: decomp_XYZ
        real(mytype), dimension(decomp_XYZ%xst(1):decomp_XYZ%xen(1), decomp_XYZ%xst(2):decomp_XYZ%xen(2),     &
        decomp_XYZ%xst(3):decomp_XYZ%xen(3))   :: sca_x_XYZ

        file_path=trim(fields_dir)//"/sca1"
        call hdf_read_3Dfield(file_path, sca_x_XYZ, "sca1", decomp_XYZ%xsz(1),decomp_XYZ%ysz(2),decomp_XYZ%zsz(3), decomp_XYZ%xst(1),decomp_XYZ%xen(1),  decomp_XYZ%xst(2),decomp_XYZ%xen(2), decomp_XYZ%xst(3),decomp_XYZ%xen(3))

    end subroutine read_fields

    subroutine write_fields(fields_dir)

        use following_physical_fields
        use following_mesh
        use DNS_settings
        use HDF5_IO
        use following_data

        implicit none
        character(*)    :: fields_dir

        character(200)    :: file_path

        file_path=trim(fields_dir)//"/sca1"
        call hdf_write_3Dfield(file_path, sca_x(:,:,:), "sca1", n1, n2, n3, xstart_f(1),xend_f(1),xstart_f(2),xend_f(2),xstart_f(3),xend_f(3))

    end subroutine write_fields

end module following_scalar_dao

module following_scalar_loader

    use mpi
    use decomp_2d
    use decomp_2d_io
    use following_scalar_data

    implicit none

contains


    subroutine load_fields(file, fill_from_coarse, fexist)

        use following_scalar_data
        use following_mesh
        use following_scalar_dao, only: DAO_read_fields=>read_fields

        implicit none

        character(*), intent(in)        :: file
        logical, intent(in)             :: fill_from_coarse
        logical, optional               :: fexist(4)

        type(DECOMP_INFO)   :: decomp_XYZ
        integer             :: n1c, n2c, n3c        ! Coarse mesh resolution
        real*8, dimension(:,:,:), allocatable   :: sca_x_XYZ
        integer             :: i


        call decomp_info_init(n1, n2, n3, decomp_XYZ)
        call DAO_read_fields(file, sca_x(:,:,:), decomp_XYZ, fexist)

        ! Spread to all transpositions

        call transpose_x_to_y(sca_x(:,:,:), sca_y(:,:,:), decomp_XYZ)
        call transpose_y_to_z(sca_y(:,:,:), sca_z(:,:,:), decomp_XYZ)

    end subroutine load_fields

    subroutine init_wall_sca()
        use following_scalar_data

        implicit none

        if(nrank==0) write(*,*)'### Following channel SCALAR_generate_fields at wall'

        sca_wall10=0.d0
        sca_wall11=0.d0

        sca_wall20=-delta_T
        sca_wall21=delta_T

        sca_wall30=0.d0
        sca_wall31=0.d0

    end subroutine init_wall_sca

end module following_scalar_loader

module following_scalar_results_writer

    use mpi
    use decomp_2d
    use decomp_2d_io
    use following_scalar_data

    implicit none

contains

    ! Export all physical field (velocity and scalar) in HDF5 format in individual files
    ! More generic, better than 2decomp format (not generic)
    subroutine write_fields(it)

        use run_ctxt_data
        use following_common_workspace_view
        use following_scalar_dao, only: DAO_write_fields=>write_fields

        implicit none

        integer, intent(in) :: it

        character*10 tmp_str
        character*80 result_dir_path

        write(tmp_str, "(i10)")it
        result_dir_path=trim(following_common_results3D_path)//'field'//trim(adjustl(tmp_str))

        call DAO_write_fields(result_dir_path)

    end subroutine write_fields

end module following_scalar_results_writer
