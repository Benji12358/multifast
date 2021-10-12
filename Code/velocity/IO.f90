module VELOCITY_settings
    implicit none

contains

    subroutine read_blowing_settings()

        use COMMON_workspace_view
        use boundaries
        use blow_settings

        implicit none
        integer         :: s

        open(15,file=trim(COMMON_settings_path)//'blowing.d')

        read(15,*) nb_slows
        read(15,*) blow_end

        do s = 1, nb_slows
            read(15,*) slows(s)%xst, slows(s)%xen, slows(s)%zst, slows(s)%zen, slows(s)%blowing
        end do

        close(15)


    end subroutine read_blowing_settings

    subroutine read_twave_settings()

        use COMMON_workspace_view
        use boundaries
        use twave_settings

        implicit none

        open(15,file=trim(COMMON_settings_path)//'travelling_wave.d')

        read(15,*) twave_on
        !read(15,*) tstart
        read(15,*) inner_units
        read(15,*) Re_tau
        read(15,*) twave%Amp
        read(15,*) twave%kappa
        read(15,*) twave%omega

        close(15)

    end subroutine read_twave_settings


end module VELOCITY_settings

module VELOCITY_dao

    use mpi
    use decomp_2d
    use decomp_2d_io

    implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!      Below lies our implementation that write all variables       !!!!!!!!
!!!!!!!!               at the cell centered in a single hdf5               !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hdf5output(nstep,time)

    use physical_fields
    use mesh
    use DNS_settings
    use HDF5_IO
    use VELOCITY_operations
    use VELOCITY_workspace_view
    use SCALAR_data
    implicit none

    type(datafile) :: file
    integer, intent(in) :: nstep
    real*8 , intent(in) :: time

    character*500 suffix
    character*500 filename, path, fnam
    integer :: error,ppos
    integer(hsize_t) :: dimsf(3)

    integer :: n1s,n1e,n2s,n2e,n3s,n3e,ierrr,myid
    real*8, allocatable, dimension(:,:,:)  :: q1c_x, q2c_y, q3c_z, scaC_y

    call mpi_comm_rank(MPI_COMM_WORLD,myid,ierrr)

    allocate(q1c_x(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
    allocate(q2c_y(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
    allocate(q3c_z(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))

    allocate(scaC_y(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))

    ! VALGRIND
    q1c_x=0.d0
    q2c_y=0.d0
    q3c_z=0.d0

    call perform_velocity_at_center(q3_z, q2_y, q1_x, q3c_z, q2c_y, q1c_x)

    dimsf(1:3)=(/nx_global-1,ny_global-1,nz_global-1/)

    path = trim(results_path)//'cell_centered_fields/'
    call system('mkdir -p '//trim(path))
                                     
    write(suffix,'(I0.8)') nstep                                           
    filename=trim(path)//'fields_'//trim(suffix)//'.h5'

    ! Create and Open HDF5 file
    call create_file(filename,file)

    ! Write simulation parameters: Re and Mach
    call h5gcreate_f(file%id,'setup',file%current_group, ierrr)
    call write_hdf(file,ren,'Re')
    call h5gclose_f(file%current_group, ierrr)

    ! Write Grid
    call h5gcreate_f(file%id,'grid',file%current_group,ierrr)
    call write_hdf(file,nx_global-1,'nx')
    call write_hdf(file,ny_global-1,'ny')
    call write_hdf(file,nz_global-1,'nz')
    call write_hdf(file,time,'time')
    call write_hdf(file,nstep,'step')
    call write_hdf(file,Xc,'Xc')
    call write_hdf(file,Yc,'Yc')
    call write_hdf(file,Zc,'Zc')
    call h5gclose_f(file%current_group,ierrr)

    ! Write density,pressure and temperature fields 
    call h5gcreate_f(file%id,'fields',file%current_group, ierrr)

    if (SCA_state==1) then
        n1s = ystart(1); n1e=min(yend(1),nx_global-1)
        n2s = ystart(2); n2e=min(yend(2),ny_global-1)
        n3s = ystart(3); n3e=min(yend(3),nz_global-1)

        call perform_passive_scalar_at_center(sca_y, scaC_y, delta_T)

        call write_hdf(file,scaC_y(n1s:n1e,n2s:n2e,n3s:n3e),dimsf,'temperature',n1s,n1e,n2s,n2e,n3s,n3e)
    endif

    n1s = xstart(1); n1e=min(xend(1),nx_global-1)
    n2s = xstart(2); n2e=min(xend(2),ny_global-1)
    n3s = xstart(3); n3e=min(xend(3),nz_global-1)

    call write_hdf(file,dp_x(n1s:n1e,n2s:n2e,n3s:n3e),dimsf,'pressure',n1s,n1e,n2s,n2e,n3s,n3e)
    call h5gclose_f(file%current_group, ierrr)

    ! Write velocity fields
    call h5gcreate_f(file%id, "fields/velocity",file%current_group, ierrr)

    n1s = xstart(1); n1e=min(xend(1),nx_global-1)
    n2s = xstart(2); n2e=min(xend(2),ny_global-1)
    n3s = xstart(3); n3e=min(xend(3),nz_global-1)
    call write_hdf(file,q1c_x(n1s:n1e,n2s:n2e,n3s:n3e),dimsf,'W',n1s,n1e,n2s,n2e,n3s,n3e)

    n1s = ystart(1); n1e=min(yend(1),nx_global-1)
    n2s = ystart(2); n2e=min(yend(2),ny_global-1)
    n3s = ystart(3); n3e=min(yend(3),nz_global-1)
    call write_hdf(file,q2c_y(n1s:n1e,n2s:n2e,n3s:n3e),dimsf,'V',n1s,n1e,n2s,n2e,n3s,n3e)

    n1s = zstart(1); n1e=min(zend(1),nx_global-1)
    n2s = zstart(2); n2e=min(zend(2),ny_global-1)
    n3s = zstart(3); n3e=min(zend(3),nz_global-1)
    call write_hdf(file,q3c_z(n1s:n1e,n2s:n2e,n3s:n3e),dimsf,'U',n1s,n1e,n2s,n2e,n3s,n3e)

    call h5gclose_f(file%current_group, ierrr)

    ! Close HDF5 file
    call close_file(file)

    ! Write XDMF file for visualization in Paraview or Visit
    if (myid==0) then
      dimsf(1:3)=(/nz_global-1,ny_global-1,nx_global-1/)
      !write(suffix,'(I0.8)') step
      !fnam='fields_'//suffix
      filename='fields_'//suffix
      call write_xdmf(trim(path),trim(filename),time,dimsf,SCA_state)
      print '(3A)', '# Written HDF5/XDMF files to disk: ',trim(filename),'.{h5,xmf}'
    endif

    deallocate(q1c_x)
    deallocate(q2c_y)
    deallocate(q3c_z)

    deallocate(scaC_y)

  end subroutine hdf5output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Export all velocity fields in HDF5 format in individual files
    ! More generic, better than 2decomp format (not generic)
    subroutine write_fields(fields_dir)

        use physical_fields
        use mesh
        use DNS_settings
        use HDF5_IO

        implicit none
        character(*)    :: fields_dir


        character(200)    :: file_path

        file_path=trim(fields_dir)//"/U"
        call hdf_write_3Dfield(file_path, q3_z(:,:,:), "U", nx_global, ny_global, nz_global, zstart(1),zend(1),zstart(2),zend(2),zstart(3),zend(3))

        file_path=trim(fields_dir)//"/V"
        call hdf_write_3Dfield(file_path, q2_y(:,:,:), "V", nx_global, ny_global, nz_global, ystart(1),yend(1),ystart(2),yend(2),ystart(3),yend(3))

        file_path=trim(fields_dir)//"/W"
        call hdf_write_3Dfield(file_path, q1_x(:,:,:), "W", nx_global, ny_global, nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))

        file_path=trim(fields_dir)//"/P"
        call hdf_write_3Dfield(file_path, dp_x(:,:,:), "P", nx_global, ny_global, nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))

!        file_path=trim(fields_dir)//"/void"
!        call hdf_write_3Dfield(file_path, void_x(:,:,:), "void", nx_global, ny_global, nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))

        ! Save boundary conditions
        file_path=trim(fields_dir)//"/Wall"
        if(nrank==0)  call hdf_create_emptyfile(file_path)
        if(nrank==0)  call hdf_addgroup(file_path, "Wall10")
        call hdf_add_2Dfield(file_path, q1_wall10(:,:), "Wall10/q1", ny_global, nz_global, xstart(2),xend(2), xstart(3),xend(3))
        call hdf_add_2Dfield(file_path, q2_wall10(:,:), "Wall10/q2", ny_global, nz_global, xstart(2),xend(2), xstart(3),xend(3))
        call hdf_add_2Dfield(file_path, q3_wall10(:,:), "Wall10/q3", ny_global, nz_global, xstart(2),xend(2), xstart(3),xend(3))
        if(nrank==0)  call hdf_addgroup(file_path, "Wall11")
        call hdf_add_2Dfield(file_path, q1_wall11(:,:), "Wall11/q1", ny_global, nz_global, xstart(2),xend(2), xstart(3),xend(3))
        call hdf_add_2Dfield(file_path, q2_wall11(:,:), "Wall11/q2", ny_global, nz_global, xstart(2),xend(2), xstart(3),xend(3))
        call hdf_add_2Dfield(file_path, q3_wall11(:,:), "Wall11/q3", ny_global, nz_global, xstart(2),xend(2), xstart(3),xend(3))

        if(nrank==0)  call hdf_addgroup(file_path, "Wall20")
        call hdf_add_2Dfield(file_path, q1_wall20(:,:), "Wall20/q1", nx_global, nz_global, ystart(1),yend(1), ystart(3),yend(3))
        call hdf_add_2Dfield(file_path, q2_wall20(:,:), "Wall20/q2", nx_global, nz_global, ystart(1),yend(1), ystart(3),yend(3))
        call hdf_add_2Dfield(file_path, q3_wall20(:,:), "Wall20/q3", nx_global, nz_global, ystart(1),yend(1), ystart(3),yend(3))
        if(nrank==0)  call hdf_addgroup(file_path, "Wall21")
        call hdf_add_2Dfield(file_path, q1_wall21(:,:), "Wall21/q1", nx_global, nz_global, ystart(1),yend(1), ystart(3),yend(3))
        call hdf_add_2Dfield(file_path, q2_wall21(:,:), "Wall21/q2", nx_global, nz_global, ystart(1),yend(1), ystart(3),yend(3))
        call hdf_add_2Dfield(file_path, q3_wall21(:,:), "Wall21/q3", nx_global, nz_global, ystart(1),yend(1), ystart(3),yend(3))

        if(nrank==0)  call hdf_addgroup(file_path, "Wall30")
        call hdf_add_2Dfield(file_path, q1_wall30(:,:), "Wall30/q1", nx_global, ny_global, zstart(1),zend(1), zstart(2),zend(2))
        call hdf_add_2Dfield(file_path, q2_wall30(:,:), "Wall30/q2", nx_global, ny_global, zstart(1),zend(1), zstart(2),zend(2))
        call hdf_add_2Dfield(file_path, q3_wall30(:,:), "Wall30/q3", nx_global, ny_global, zstart(1),zend(1), zstart(2),zend(2))
        if(nrank==0)  call hdf_addgroup(file_path, "Wall31")
        call hdf_add_2Dfield(file_path, q1_wall31(:,:), "Wall31/q1", nx_global, ny_global, zstart(1),zend(1), zstart(2),zend(2))
        call hdf_add_2Dfield(file_path, q2_wall31(:,:), "Wall31/q2", nx_global, ny_global, zstart(1),zend(1), zstart(2),zend(2))
        call hdf_add_2Dfield(file_path, q3_wall31(:,:), "Wall31/q3", nx_global, ny_global, zstart(1),zend(1), zstart(2),zend(2))

    end subroutine write_fields


    ! Read all velocity fields in HDF5 format in individual files
    subroutine read_fields(fields_dir, u_x_XYZ, v_x_XYZ, w_x_XYZ, dp_x_XYZ, decomp_XYZ, files_exist)

        use HDF5_IO
        use physical_fields
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


        file_path=trim(fields_dir)//"/Wall"
        call hdf_read_2Dfield(file_path, q1_wall10(:,:), "Wall10/q1", decomp_XYZ%yen(2), decomp_XYZ%zen(3), decomp_XYZ%xst(2),decomp_XYZ%xen(2), decomp_XYZ%xst(3),decomp_XYZ%xen(3))
        call hdf_read_2Dfield(file_path, q2_wall10(:,:), "Wall10/q2", decomp_XYZ%yen(2), decomp_XYZ%zen(3), decomp_XYZ%xst(2),decomp_XYZ%xen(2), decomp_XYZ%xst(3),decomp_XYZ%xen(3))
        call hdf_read_2Dfield(file_path, q3_wall10(:,:), "Wall10/q3", decomp_XYZ%yen(2), decomp_XYZ%zen(3), decomp_XYZ%xst(2),decomp_XYZ%xen(2), decomp_XYZ%xst(3),decomp_XYZ%xen(3))

        call hdf_read_2Dfield(file_path, q1_wall11(:,:), "Wall11/q1", decomp_XYZ%yen(2), decomp_XYZ%zen(3), decomp_XYZ%xst(2),decomp_XYZ%xen(2), decomp_XYZ%xst(3),decomp_XYZ%xen(3))
        call hdf_read_2Dfield(file_path, q2_wall11(:,:), "Wall11/q2", decomp_XYZ%yen(2), decomp_XYZ%zen(3), decomp_XYZ%xst(2),decomp_XYZ%xen(2), decomp_XYZ%xst(3),decomp_XYZ%xen(3))
        call hdf_read_2Dfield(file_path, q3_wall11(:,:), "Wall11/q3", decomp_XYZ%yen(2), decomp_XYZ%zen(3), decomp_XYZ%xst(2),decomp_XYZ%xen(2), decomp_XYZ%xst(3),decomp_XYZ%xen(3))

        call hdf_read_2Dfield(file_path, q1_wall20(:,:), "Wall20/q1", decomp_XYZ%xen(1), decomp_XYZ%zen(3), decomp_XYZ%yst(1),decomp_XYZ%yen(1), decomp_XYZ%yst(3),decomp_XYZ%yen(3))
        call hdf_read_2Dfield(file_path, q2_wall20(:,:), "Wall20/q2", decomp_XYZ%xen(1), decomp_XYZ%zen(3), decomp_XYZ%yst(1),decomp_XYZ%yen(1), decomp_XYZ%yst(3),decomp_XYZ%yen(3))
        call hdf_read_2Dfield(file_path, q3_wall20(:,:), "Wall20/q3", decomp_XYZ%xen(1), decomp_XYZ%zen(3), decomp_XYZ%yst(1),decomp_XYZ%yen(1), decomp_XYZ%yst(3),decomp_XYZ%yen(3))

        call hdf_read_2Dfield(file_path, q1_wall21(:,:), "Wall21/q1", decomp_XYZ%xen(1), decomp_XYZ%zen(3), decomp_XYZ%yst(1),decomp_XYZ%yen(1), decomp_XYZ%yst(3),decomp_XYZ%yen(3))
        call hdf_read_2Dfield(file_path, q2_wall21(:,:), "Wall21/q2", decomp_XYZ%xen(1), decomp_XYZ%zen(3), decomp_XYZ%yst(1),decomp_XYZ%yen(1), decomp_XYZ%yst(3),decomp_XYZ%yen(3))
        call hdf_read_2Dfield(file_path, q3_wall21(:,:), "Wall21/q3", decomp_XYZ%xen(1), decomp_XYZ%zen(3), decomp_XYZ%yst(1),decomp_XYZ%yen(1), decomp_XYZ%yst(3),decomp_XYZ%yen(3))

        call hdf_read_2Dfield(file_path, q1_wall30(:,:), "Wall30/q1", decomp_XYZ%xen(1), decomp_XYZ%yen(2), decomp_XYZ%zst(1),decomp_XYZ%zen(1), decomp_XYZ%zst(2),decomp_XYZ%zen(2))
        call hdf_read_2Dfield(file_path, q2_wall30(:,:), "Wall30/q2", decomp_XYZ%xen(1), decomp_XYZ%yen(2), decomp_XYZ%zst(1),decomp_XYZ%zen(1), decomp_XYZ%zst(2),decomp_XYZ%zen(2))
        call hdf_read_2Dfield(file_path, q3_wall30(:,:), "Wall30/q3", decomp_XYZ%xen(1), decomp_XYZ%yen(2), decomp_XYZ%zst(1),decomp_XYZ%zen(1), decomp_XYZ%zst(2),decomp_XYZ%zen(2))

        call hdf_read_2Dfield(file_path, q1_wall31(:,:), "Wall31/q1", decomp_XYZ%xen(1), decomp_XYZ%yen(2), decomp_XYZ%zst(1),decomp_XYZ%zen(1), decomp_XYZ%zst(2),decomp_XYZ%zen(2))
        call hdf_read_2Dfield(file_path, q2_wall31(:,:), "Wall31/q2", decomp_XYZ%xen(1), decomp_XYZ%yen(2), decomp_XYZ%zst(1),decomp_XYZ%zen(1), decomp_XYZ%zst(2),decomp_XYZ%zen(2))
        call hdf_read_2Dfield(file_path, q3_wall31(:,:), "Wall31/q3", decomp_XYZ%xen(1), decomp_XYZ%yen(2), decomp_XYZ%zst(1),decomp_XYZ%zen(1), decomp_XYZ%zst(2),decomp_XYZ%zen(2))

    end subroutine read_fields

end module VELOCITY_dao


module VELOCITY_loader
    use mpi
    use decomp_2d
    use decomp_2d_io

    implicit none

contains



    ! ATTENTION BC3 remplacé par NS_DEF_BC3 (pareil pour les autres directions)

    ! Interpolation of data from a coarse grid ((N+1)/2 points in each direction) to the current one
    ! used for recovering a simulation from a coarse grid
    subroutine build_fine_field(q1_x_XYZ, q2_x_XYZ, q3_x_XYZ, press_3D_x_XYZ, q1_z, q2_z, q3_z, press_3D_z, decomp_XYZ)
        use schemes_interface
        use boundaries

        use mesh
!        use HDF5_IO

        implicit none

        type(DECOMP_INFO)   :: decomp_XYZ, decomp_YZ, decomp_Z

        real*8, dimension(decomp_XYZ%xst(1):decomp_XYZ%xen(1), decomp_XYZ%xst(2):decomp_XYZ%xen(2), decomp_XYZ%xst(3):decomp_XYZ%xen(3))   :: q3_x_XYZ, q2_x_XYZ, q1_x_XYZ, press_3D_x_XYZ
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))   :: q3_z, q2_z, q1_z, press_3D_z
        integer                     :: n1c, n2c, n3c

        real*8, dimension(:,:,:), allocatable   :: q1_x_YZ, q2_x_YZ, q3_x_YZ, press_3D_x_YZ
        real*8, dimension(:,:,:), allocatable   :: q1_y_YZ, q2_y_YZ, q3_y_YZ, press_3D_y_YZ
        real*8, dimension(:,:,:), allocatable   :: q1_y_Z, q2_y_Z, q3_y_Z, press_3D_y_Z
        real*8, dimension(:,:,:), allocatable   :: q1_z_Z, q2_z_Z, q3_z_Z, press_3D_z_Z
        !real*8, dimension(:,:,:), allocatable   :: q1_y

        real*8, dimension(4, decomp_XYZ%xsz(1))    :: tmp_x_c
        real*8, dimension(4, 2*decomp_XYZ%xsz(1)-1)    :: tmp_x_f
        real*8, dimension(4, decomp_XYZ%ysz(2))    :: tmp_y_c
        real*8, dimension(4, 2*decomp_XYZ%ysz(2)-1)    :: tmp_y_f
        real*8, dimension(4, decomp_XYZ%zsz(3))    :: tmp_z_c
        real*8, dimension(4, 2*decomp_XYZ%zsz(3)-1)    :: tmp_z_f

        integer :: i, j, k

!        integer :: n1s, n1e, n2s, n2e, n3s, n3e

        n1c=decomp_XYZ%xsz(1)
        n2c=decomp_XYZ%ysz(2)
        n3c=decomp_XYZ%zsz(3)


        call decomp_info_init(xsize(1), n2c, n3c, decomp_YZ)
        call decomp_info_init(xsize(1), ysize(2), n3c, decomp_Z)

        ! X INTERPOLATION_________________________________________________________________________________________________
        allocate(q1_x_YZ(decomp_YZ%xst(1):decomp_YZ%xen(1), decomp_YZ%xst(2):decomp_YZ%xen(2), decomp_YZ%xst(3):decomp_YZ%xen(3)))
        allocate(q2_x_YZ(decomp_YZ%xst(1):decomp_YZ%xen(1), decomp_YZ%xst(2):decomp_YZ%xen(2), decomp_YZ%xst(3):decomp_YZ%xen(3)))
        allocate(q3_x_YZ(decomp_YZ%xst(1):decomp_YZ%xen(1), decomp_YZ%xst(2):decomp_YZ%xen(2), decomp_YZ%xst(3):decomp_YZ%xen(3)))
        allocate(press_3D_x_YZ(decomp_YZ%xst(1):decomp_YZ%xen(1), decomp_YZ%xst(2):decomp_YZ%xen(2), decomp_YZ%xst(3):decomp_YZ%xen(3)))

        if (n1c/=n1) then
            do k=decomp_XYZ%xst(3), decomp_XYZ%xen(3)! min(n3m, decomp_XYZ%xen(3))
                do j=decomp_XYZ%xst(2), decomp_XYZ%xen(2)! min(n2, decomp_XYZ%xen(2))
                    call D0s(q1_x_XYZ(:,j,k), tmp_x_c(1, :), n1c, 0.d0, .false., NS_DEF_BC3)
                    call D0s(q2_x_XYZ(:,j,k), tmp_x_c(2, :), n1c, 0.d0, .true., NS_DEF_BC3)
                    call D0s(q3_x_XYZ(:,j,k), tmp_x_c(3, :), n1c, 0.d0, .true., NS_DEF_BC3)
                    call D0s(press_3D_x_XYZ(:,j,k), tmp_x_c(4, :), n1c, 0.d0, .true., NS_DEF_BC3)

                    do i = 1, n1c
                        if(i/=n1c) q1_x_YZ(2*i, j, k)=tmp_x_c(1, i)
                        q1_x_YZ(2*i-1, j, k)=q1_x_XYZ(i, j, k)

                        if(i/=n1c) tmp_x_f(2, 2*i)=q2_x_XYZ(i, j, k)
                        tmp_x_f(2,2*i-1)=tmp_x_c(2, i)

                        if(i/=n1c) tmp_x_f(3, 2*i)=q3_x_XYZ(i, j, k)
                        tmp_x_f(3, 2*i-1)=tmp_x_c(3, i)

                        if(i/=n1c) tmp_x_f(4, 2*i)=press_3D_x_XYZ(i, j, k)
                        tmp_x_f(4, 2*i-1)=tmp_x_c(4, i)
                    end do

                    call D0s(tmp_x_f(2, :), q2_x_YZ(:,j,k), n1, 0.d0, .false., NS_DEF_BC3)
                    call D0s(tmp_x_f(3, :), q3_x_YZ(:,j,k), n1, 0.d0, .false., NS_DEF_BC3)
                    call D0s(tmp_x_f(4, :), press_3D_x_YZ(:,j,k), n1, 0.d0, .false., NS_DEF_BC3)

                enddo
            enddo
        else
            q1_x_YZ=q1_x_XYZ
            q2_x_YZ=q2_x_XYZ
            q3_x_YZ=q3_x_XYZ
            press_3D_x_YZ=press_3D_x_XYZ
        endif

        allocate(q1_y_YZ(decomp_YZ%yst(1):decomp_YZ%yen(1), decomp_YZ%yst(2):decomp_YZ%yen(2), decomp_YZ%yst(3):decomp_YZ%yen(3)))
        allocate(q2_y_YZ(decomp_YZ%yst(1):decomp_YZ%yen(1), decomp_YZ%yst(2):decomp_YZ%yen(2), decomp_YZ%yst(3):decomp_YZ%yen(3)))
        allocate(q3_y_YZ(decomp_YZ%yst(1):decomp_YZ%yen(1), decomp_YZ%yst(2):decomp_YZ%yen(2), decomp_YZ%yst(3):decomp_YZ%yen(3)))
        allocate(press_3D_y_YZ(decomp_YZ%yst(1):decomp_YZ%yen(1), decomp_YZ%yst(2):decomp_YZ%yen(2), decomp_YZ%yst(3):decomp_YZ%yen(3)))

        call transpose_x_to_y(q1_x_YZ, q1_y_YZ, decomp_YZ)
        call transpose_x_to_y(q2_x_YZ, q2_y_YZ, decomp_YZ)
        call transpose_x_to_y(q3_x_YZ, q3_y_YZ, decomp_YZ)
        call transpose_x_to_y(press_3D_x_YZ, press_3D_y_YZ, decomp_YZ)

        deallocate(q1_x_YZ, q2_x_YZ, q3_x_YZ, press_3D_x_YZ)

        ! Y INTERPOLATION
        allocate(q1_y_Z(decomp_Z%yst(1):decomp_Z%yen(1), decomp_Z%yst(2):decomp_Z%yen(2), decomp_Z%yst(3):decomp_Z%yen(3)))
        allocate(q2_y_Z(decomp_Z%yst(1):decomp_Z%yen(1), decomp_Z%yst(2):decomp_Z%yen(2), decomp_Z%yst(3):decomp_Z%yen(3)))
        allocate(q3_y_Z(decomp_Z%yst(1):decomp_Z%yen(1), decomp_Z%yst(2):decomp_Z%yen(2), decomp_Z%yst(3):decomp_Z%yen(3)))
        allocate(press_3D_y_Z(decomp_Z%yst(1):decomp_Z%yen(1), decomp_Z%yst(2):decomp_Z%yen(2), decomp_Z%yst(3):decomp_Z%yen(3)))

        if (n2c/=n2) then
            do k=decomp_YZ%yst(3), decomp_YZ%yen(3) !min(n3m, decomp_YZ%yen(3))
                do i=ystart(1), yend(1) !min(n1m, yend(1))

                    call D0s(q3_y_YZ(i,:,k), tmp_y_c(3, :), n2c, 0.d0, .true., NS_DEF_BC2)
                    call D0s(q1_y_YZ(i,:,k), tmp_y_c(1, :), n2c, 0.d0, .true., NS_DEF_BC2)
                    call D0s(press_3D_y_YZ(i,:,k), tmp_y_c(4, :), n2c, 0.d0, .true., NS_DEF_BC2)

                    tmp_y_c(:, 1)=0.d0
                    tmp_y_c(:, n2c)=0.d0
                    do j = 1, n2c-1
                        tmp_y_f(3, 2*j)=q3_y_YZ(i, j, k)
                        tmp_y_f(3, 2*j-1)=tmp_y_c(3, j)

                        tmp_y_f(1, 2*j)=q1_y_YZ(i, j, k)
                        tmp_y_f(1, 2*j-1)=tmp_y_c(1, j)

                        tmp_y_f(4, 2*j)=press_3D_y_YZ(i, j, k)
                        tmp_y_f(4, 2*j-1)=tmp_y_c(4, j)
                    end do

                    tmp_y_f(3, 2*n2c-1)=tmp_y_c(3, n2c)
                    tmp_y_f(1, 2*n2c-1)=tmp_y_c(1, n2c)
                    tmp_y_f(4, 2*n2c-1)=tmp_y_c(4, n2c)

                    call D0s(tmp_y_f(3, :), q3_y_Z(i,:,k), n2, 0.d0, .false., NS_DEF_BC2)
                    call D0s(tmp_y_f(1, :), q1_y_Z(i,:,k), n2, 0.d0, .false., NS_DEF_BC2)
                    call D0s(tmp_y_f(4, :), press_3D_y_Z(i,:,k), n2, 0.d0, .false., NS_DEF_BC2)

                end do
            end do

            do k=decomp_YZ%yst(3), decomp_YZ%yen(3) !min(n3m, decomp_YZ%yen(3))
                do i=ystart(1), yend(1) !min(n1m, yend(1))

                    call D0s(q2_y_YZ(i,:,k), tmp_y_c(2, :), n2c, 0.d0, .false., NS_DEF_BC2)

                    do j = 1, n2c-1
                        q2_y_Z(i, 2*j-1, k)=q2_y_YZ(i, j, k)
                        q2_y_Z(i, 2*j, k)=tmp_y_c(2, j)
                    end do

                    q2_y_Z(i, 2*n2c-1, k)=q2_y_YZ(i, n2c, k)

                end do
            end do
        else
            q1_y_Z=q1_y_YZ
            q2_y_Z=q2_y_YZ
            q3_y_Z=q3_y_YZ
            press_3D_y_Z=press_3D_y_YZ
        endif

        allocate(q1_z_Z(decomp_Z%zst(1):decomp_Z%zen(1), decomp_Z%zst(2):decomp_Z%zen(2), decomp_Z%zst(3):decomp_Z%zen(3)))
        allocate(q2_z_Z(decomp_Z%zst(1):decomp_Z%zen(1), decomp_Z%zst(2):decomp_Z%zen(2), decomp_Z%zst(3):decomp_Z%zen(3)))
        allocate(q3_z_Z(decomp_Z%zst(1):decomp_Z%zen(1), decomp_Z%zst(2):decomp_Z%zen(2), decomp_Z%zst(3):decomp_Z%zen(3)))
        allocate(press_3D_z_Z(decomp_Z%zst(1):decomp_Z%zen(1), decomp_Z%zst(2):decomp_Z%zen(2), decomp_Z%zst(3):decomp_Z%zen(3)))

        call transpose_y_to_z(q1_y_Z, q1_z_Z, decomp_Z)
        call transpose_y_to_z(q2_y_Z, q2_z_Z, decomp_Z)
        call transpose_y_to_z(q3_y_Z, q3_z_Z, decomp_Z)
        call transpose_y_to_z(press_3D_y_Z, press_3D_z_Z, decomp_Z)

        deallocate(q1_y_Z, q2_y_Z, q3_y_Z, press_3D_y_Z)
        deallocate(q1_y_YZ, q2_y_YZ, q3_y_YZ, press_3D_y_YZ)

        ! Z INTERPOLATION ______________________________________________________

        if(n3c/=n3) then
            do j = decomp_Z%zst(2), decomp_Z%zen(2) ! min(n2, decomp_Z%zen(2))
                do i = decomp_Z%zst(1), decomp_Z%zen(1) ! min(n1-1, decomp_Z%zen(1))
                    call D0s(q1_z_Z(i,j,:), tmp_z_c(1, :), n3c, 0.d0, .true., NS_DEF_BC1)
                    call D0s(q2_z_Z(i,j,:), tmp_z_c(2, :), n3c, 0.d0, .true., NS_DEF_BC1)
                    call D0s(q3_z_Z(i,j,:), tmp_z_c(3, :), n3c, 0.d0, .false., NS_DEF_BC1)
                    call D0s(press_3D_z_Z(i,j,:), tmp_z_c(4, :), n3c, 0.d0, .true., NS_DEF_BC1)

                    do k = 1, n3c

                        if(k/=n3c) tmp_z_f(1, 2*k)=q1_z_Z(i, j, k)
                        tmp_z_f(1, 2*k-1)=tmp_z_c(1, k)

                        if(k/=n3c) tmp_z_f(2, 2*k)=q2_z_Z(i, j, k)
                        tmp_z_f(2, 2*k-1)=tmp_z_c(2, k)

                        if(k/=n3c) q3_z(i, j, 2*k)=tmp_z_c(3, k)
                        q3_z(i, j, 2*k-1)=q3_z_Z(i, j, k)

                        if(k/=n3c) tmp_z_f(4, 2*k)=press_3D_z_Z(i, j, k)
                        tmp_z_f(4, 2*k-1)=tmp_z_c(4, k)

                    end do

                    call D0s(tmp_z_f(1, :), q1_z(i,j,:), n3, 0.d0, .false., NS_DEF_BC1)
                    call D0s(tmp_z_f(2, :), q2_z(i,j,:), n3, 0.d0, .false., NS_DEF_BC1)
                    call D0s(tmp_z_f(4, :), press_3D_z(i,j,:), n3, 0.d0, .false., NS_DEF_BC1)

                end do
            end do
        else
            q1_z=q1_z_Z
            q2_z=q2_z_Z
            q3_z=q3_z_Z
            press_3D_z=press_3D_z_Z
        endif

!        if(nrank==0)  call hdf_create_file_with_3Dmesh('Fine_fields', Z, Y,X, "x1", "x2", "x3", n1, n2,n3)
!
!        if(nrank==0)  call hdf_addgroup('Fine_fields', "q")
!
!        n1s=zstart(1)
!        n1e=zend(1)
!        n2s=zstart(2)
!        n2e=zend(2)
!        n3s=zstart(3)
!        n3e=zend(3)
!
!        call hdf_add_3Dfield('Fine_fields', q1_z(n1s:n1e,n2s:n2e,n3s:n3e), "q/q1", n1, n2, n3, n1s, n1e, n2s, n2e, n3s, n3e)
!        call hdf_add_3Dfield('Fine_fields', q2_z(n1s:n1e,n2s:n2e,n3s:n3e), "q/q2", n1, n2, n3, n1s, n1e, n2s, n2e, n3s, n3e)
!        call hdf_add_3Dfield('Fine_fields', q3_z(n1s:n1e,n2s:n2e,n3s:n3e), "q/q3", n1, n2, n3, n1s, n1e, n2s, n2e, n3s, n3e)
!        call hdf_add_3Dfield('Fine_fields', press_3D_z(n1s:n1e,n2s:n2e,n3s:n3e), "pr", n1, n2, n3, n1s, n1e, n2s, n2e, n3s, n3e)

        deallocate(q1_z_Z, q2_z_Z, q3_z_Z, press_3D_z_Z)


    end subroutine build_fine_field


    subroutine load_fields(file, fill_from_coarse, fexist)

        use physical_fields
        use mesh
        use DNS_settings
        use start_settings


        use VELOCITY_dao, only: DAO_read_fields => read_fields

        implicit none

        character(*), intent(in)        :: file
        logical, intent(in)             :: fill_from_coarse
        logical, optional               :: fexist(4)

        type(DECOMP_INFO)   :: decomp_XYZ
!        integer, intent(in)             :: n1c, n2c, n3c        ! Coarse mesh resolution
        real*8, dimension(:,:,:), allocatable   :: u_x_XYZ, v_x_XYZ, w_x_XYZ, dp_x_XYZ
        integer             :: i

        if (fill_from_coarse) then
            ! Coarse field dimensions :
!            n1c=(n1+1)/2
!            n2c=(n2+1)/2
!            n3c=(n3+1)/2

            call decomp_info_init(n1c, n2c, n3c, decomp_XYZ)

            allocate(w_x_XYZ(decomp_XYZ%xst(1):decomp_XYZ%xen(1), decomp_XYZ%xst(2):decomp_XYZ%xen(2), decomp_XYZ%xst(3):decomp_XYZ%xen(3)))
            allocate(v_x_XYZ(decomp_XYZ%xst(1):decomp_XYZ%xen(1), decomp_XYZ%xst(2):decomp_XYZ%xen(2), decomp_XYZ%xst(3):decomp_XYZ%xen(3)))
            allocate(u_x_XYZ(decomp_XYZ%xst(1):decomp_XYZ%xen(1), decomp_XYZ%xst(2):decomp_XYZ%xen(2), decomp_XYZ%xst(3):decomp_XYZ%xen(3)))
            allocate(dp_x_XYZ(decomp_XYZ%xst(1):decomp_XYZ%xen(1), decomp_XYZ%xst(2):decomp_XYZ%xen(2), decomp_XYZ%xst(3):decomp_XYZ%xen(3)))

            u_x_XYZ=2.d0
            v_x_XYZ=0.d0
            w_x_XYZ=0.d0
            dp_x_XYZ=0.d0

            ! Read coarse fields
            call DAO_read_fields(file, u_x_XYZ, v_x_XYZ, w_x_XYZ, dp_x_XYZ, decomp_XYZ)
            ! Build the fine field from coarse ones
            call build_fine_field(w_x_XYZ, v_x_XYZ, u_x_XYZ, dp_x_XYZ, q1_z, q2_z, q3_z, dp_z, decomp_XYZ)

            ! Spread to all transpositions
            call transpose_z_to_y(q3_z, q3_y)
            call transpose_y_to_x(q3_y, q3_x)

            call transpose_z_to_y(q2_z, q2_y)
            call transpose_y_to_x(q2_y, q2_x)

            call transpose_z_to_y(q1_z, q1_y)
            call transpose_y_to_x(q1_y, q1_x)

            call transpose_z_to_y(dp_z, dp_y)
            call transpose_y_to_x(dp_y, dp_x)

        else

            if (half_length.eq.1) then

                call decomp_info_init(previous_fringe_start, n2, n3, decomp_XYZ)
                call DAO_read_fields(file, q3_x(:previous_fringe_start,:,:), q2_x(:previous_fringe_start,:,:), q1_x(:previous_fringe_start,:,:), dp_x(:previous_fringe_start,:,:), decomp_XYZ, fexist)

                do i=previous_fringe_start,n1
                    q1_x(i,:,:) = q1_x(previous_fringe_start,:,:)
                    q2_x(i,:,:) = q2_x(previous_fringe_start,:,:)
                    q3_x(i,:,:) = q3_x(previous_fringe_start,:,:)
                    dp_x(i,:,:) = dp_x(previous_fringe_start,:,:)
                enddo

            else 

                call decomp_info_init(n1, n2, n3, decomp_XYZ)
                call DAO_read_fields(file, q3_x, q2_x, q1_x, dp_x, decomp_XYZ, fexist)

            endif

            ! Spread to all transpositions
            call transpose_x_to_y(q3_x, q3_y)
            call transpose_y_to_z(q3_y, q3_z)

            call transpose_x_to_y(q2_x, q2_y)
            call transpose_y_to_z(q2_y, q2_z)

            call transpose_x_to_y(q1_x, q1_y)
            call transpose_y_to_z(q1_y, q1_z)

            call transpose_x_to_y(dp_x, dp_y)
            call transpose_y_to_z(dp_y, dp_z)

        end if




    end subroutine load_fields

end module VELOCITY_loader


module VELOCITY_results_writer

    use mpi
    use decomp_2d
    use decomp_2d_io

    implicit none
    integer, parameter  :: SEQUENTIAL_FILE=0, DECOMP2D_FILE=1
    integer, parameter  :: BINARY_FILE=0, ASCII_FILE=1
    integer, parameter  :: SIMPLE_FILE=0, DIRECTORY=1



contains

    ! Export all velocity fields in HDF5 format in individual files
    ! More generic, better than 2decomp format (not generic)
    subroutine write_fields(it)
        use time_writer
        use VELOCITY_workspace_view
        use VELOCITY_dao, only: DAO_write_fields=> write_fields

        implicit none

        integer, intent(in)         :: it

        character*10 tmp_str
        character*80 result_dir_path

        write(tmp_str, "(i10)")it
        result_dir_path=trim(results3D_path)//'field'//trim(adjustl(tmp_str))

        call write_timefile(trim(result_dir_path)//"/advancement.d")
        call DAO_write_fields(result_dir_path)

    end subroutine write_fields

end module VELOCITY_results_writer

module VELOCITY_log_writer

    use decomp_2d

    implicit none

contains

    subroutine export_divergence_file(file_path, cflmax, divmax, div_mean, divdiff, divg1, divg2, ntime)
        implicit none
        real*8  :: cflmax, divmax, divdiff, divg1, divg2, div_mean
        integer :: ntime, divergence_file_id
        character(*)    :: file_path

        divergence_file_id=400
!7698    format(i7,',',x,8e18.9)

        open(divergence_file_id, file=trim(file_path), position="append")

        write(divergence_file_id,*) ntime, cflmax, divmax, div_mean, divdiff, divg1, divg2
        close(divergence_file_id)

    end subroutine export_divergence_file

    subroutine export_kinetic(file_path, flow_rate, spanwise_flow_rate, kinetic_energy, enstrophy, ntime)

        implicit none
        real*8          :: kinetic_energy(3), enstrophy, flow_rate, spanwise_flow_rate
        integer         :: ntime, kinetic_file_id
        character(*)    :: file_path

        kinetic_file_id=400
!7699    format(i7,',',x,7e16.9)

        open(kinetic_file_id, file=trim(file_path), position="append")

        write(kinetic_file_id,*) ntime, sum(kinetic_energy), flow_rate, spanwise_flow_rate, kinetic_energy, enstrophy
        close(kinetic_file_id)

    end subroutine export_kinetic

    subroutine export_kinetic_IBM(kinetic_dir, flow_rate, k1,k2,k3, n, ntime)
        implicit none
        integer                 :: n, ntime
        real*8, dimension(n)    :: flow_rate, k1,k2,k3
        character(*)            :: kinetic_dir

        integer                 :: i, kinetic_file_id
        character(20)           :: tmp_str
        character(200)          :: kinetic_file

        kinetic_file_id=1201

        write(tmp_str, "(i10)")ntime
        kinetic_file=trim(kinetic_dir)//'flow_rate'//trim(adjustl(tmp_str))//".csv"

        open(kinetic_file_id, file=kinetic_file)

        do i = 1, n
            write(kinetic_file_id,*)i, flow_rate(i), k1(i), k2(i), k3(i)
        end do

        close(kinetic_file_id)

    end subroutine export_kinetic_IBM

    subroutine export_velmax_IBM(velmax_dir, q1_max, q2_max, q3_max, pr_max, q1_min, q2_min, q3_min, pr_min, n, ntime)
        use IBM_data
        implicit none
        integer                 :: n, ntime
        real*8, dimension(n)    :: q1_max, q2_max, q3_max, pr_max, mask1_max, mask2_max, mask3_max
        real*8, dimension(n)    :: q1_min, q2_min, q3_min, pr_min
        character(*)            :: velmax_dir

        integer                 :: i, velmax_file_id
        character(20)           :: tmp_str
        character(200)          :: velmax_file

        velmax_file_id=1201
7700    format(i7,',',x,11e17.8)

        write(tmp_str, "(i10)")ntime
        velmax_file=trim(velmax_dir)//'velmax'//trim(adjustl(tmp_str))//".csv"

        open(velmax_file_id, file=velmax_file)

        do i = 1, n
            mask1_max(i)=maxval(IBM_mask1_x(i,:,:))
            mask2_max(i)=maxval(IBM_mask2_x(i,:,:))
            mask3_max(i)=maxval(IBM_mask3_x(i,:,:))
        end do

        do i = 1, n
            write(velmax_file_id,7700)i, q1_max(i), q2_max(i), q3_max(i), pr_max(i), q1_min(i), q2_min(i), q3_min(i), pr_min(i), &
            mask1_max(i), mask2_max(i), mask3_max(i)
        end do

        close(velmax_file_id)

    end subroutine export_velmax_IBM


    subroutine export_properties()

        use COMMON_workspace_view
        use physical_fields
        use run_ctxt_data

        use IBM_settings
        use IBM_data

        use VELOCITY_workspace_view
        use VELOCITY_operations
        use VELOCITY_properties
        implicit none
        character(200)  :: flowrate_file
        integer         :: flowrate_file_id, i

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))      :: q1c_x
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))      :: q2c_y
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))      :: q3c_z

        real*8      :: kinetic_energy(3), enstrophy, flow_rate, spanwise_flow_rate

        call perform_kinetic(q3_z, q2_y, q1_x, flow_rate,spanwise_flow_rate, kinetic_energy, enstrophy, streamwise)
        if(nrank==0) call export_kinetic(kinetic_history_file, flow_rate,spanwise_flow_rate, kinetic_energy, enstrophy, ntime)

        if (IBM_activated) then
            if(mod(ntime, nb_flow_rate).eq.0) then

                call perform_velocity_at_center(q3_z, q2_y, q1_x, q3c_z, q2c_y, q1c_x)
                call perform_kinetic_IBM(q3c_z, q2c_y, q1c_x, flow_rate_IBM, kin1_IBM, kin2_IBM, kin3_IBM)
                call perform_maxvel_IBM(q3_x, q2_x, q1_x, dp_x, q1max_IBM, q2max_IBM, q3max_IBM, prmax_IBM, q1min_IBM, q2min_IBM, q3min_IBM, prmin_IBM)


                if(nrank==0) call export_kinetic_IBM(trim(COMMON_log_path)//'/IBM/Flow_rate/', flow_rate_IBM, kin1_IBM, kin2_IBM, kin3_IBM, n1-1, ntime)
                if(nrank==0) call export_velmax_IBM(trim(COMMON_log_path)//'/IBM/Vel_max/', q1max_IBM, q2max_IBM, q3max_IBM, prmax_IBM, q1min_IBM, q2min_IBM, q3min_IBM, prmin_IBM, n1-1, ntime)

            endif
        endif

    end subroutine export_properties

    subroutine export_gradP

        use run_ctxt_data
        use physical_fields, only: dP_streamwise, dP_spanwise
        use VELOCITY_solver, only: Fext1, Fext3
        use DNS_settings, only:streamwise
        use mesh
        use mpi

        implicit none

        integer :: i, j, k, mpierr
        real*8  :: Fstream, Fspan, stot

        Fstream=0.d0
        Fspan=0.d0
        stot=0.d0

        if (streamwise==1) then
            do i=1, n1m
                do j =xstart(2), min(n2m, xend(2))
                    do k= xstart(3), min (n3m, xend(3))
                        Fstream=Fstream+Fext1(i,j,k)*dx1*dx3*cell_size_Y(j)
                        Fspan=Fspan+Fext3(i,j,k)*dx1*dx3*cell_size_Y(j)
                    enddo
                enddo
            enddo

        elseif (streamwise==3) then
            do i=1, n1m
                do j =xstart(2), min(n2m, xend(2))
                    do k= xstart(3), min (n3m, xend(3))
                        Fstream=Fstream+Fext3(i,j,k)*dx1*dx3*cell_size_Y(j)
                        Fspan=Fspan+Fext1(i,j,k)*dx1*dx3*cell_size_Y(j)
                    enddo
                enddo
            enddo
        endif

        call MPI_ALLREDUCE (Fstream, stot, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpierr)
        Fstream=stot/(L1*L2*L3)

        stot=0.d0

        call MPI_ALLREDUCE (Fspan, stot, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpierr)
        Fspan=stot/(L1*L2*L3)

        if(nrank==0) then
            open(3698, file="gradP.dat", position="append")
            write(3698,*) ntime, t, -dP_streamwise +Fstream, -dP_streamwise, -dP_spanwise+ Fspan, -dP_spanwise
            close(3698)
        endif

    end subroutine export_gradP

end module VELOCITY_log_writer
