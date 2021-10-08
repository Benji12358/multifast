

module HDF5_IO
    USE HDF5
    use MPI
    implicit none
  
    public :: create_file,read_file,close_file,write_hdf,read_hdf,write_xdmf,datafile

    type :: datafile
      integer(kind=hid_t) :: id
      integer(HID_T) :: current_group
    end type datafile
    integer :: error
    private :: error
      interface write_hdf
        module procedure write_string,write_int,write_real,write_1dreal,write_3dreal
      end interface write_hdf

      interface read_hdf
        module procedure read_int,read_real,read_1dreal,read_3dreal
      end interface read_hdf

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!      Below lies our implementation that write all variables       !!!!!!!!
!!!!!!!!               at the cell centered in a single hdf5               !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalizeHDF5()
    ! Close HDF5 environment
    call h5close_f(error)
  end subroutine finalizeHDF5

  ! Subroutine to create HDF5 file with collective IO
  subroutine create_file(filename,file)
    implicit none
    type(datafile), intent(out) :: file
    character(*), intent(in) :: filename
    integer(kind=hid_t) :: p_id, f_id

    ! MPI IO/Lustre file striping
    integer :: info         ! MPI IO Info
    integer :: lcount       ! Lustre count size
    integer :: lsize        ! Lustre stripe size
    character(len=1024) :: clcount, clsize ! Strings of LFS
    
    ! Initialize HDF5 environment
    call h5open_f(error)

    ! Set a stripe count of 4 and a stripe size of 4MB
    lcount = 4
    lsize  = 4 * 1024 * 1024
    write(clcount, '(I4)') lcount
    write(clsize, '(I8)') lsize

    call mpi_info_create(info, error)
    call mpi_info_set(info, "striping_factor", trim(clcount), error)
    call mpi_info_set(info, "striping_unit", trim(clsize), error)

    ! Set up the access properties
    call h5pcreate_f(H5P_FILE_ACCESS_F, p_id, error)
    call h5pset_fapl_mpio_f(p_id, MPI_COMM_WORLD, info, error)

    ! Create and Open HDF5 file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, f_id, error, &
                     access_prp = p_id)
    if (error .ne. 0) then
      write(0,*) 'Unable to open: ', trim(filename), ': ', error
      call mpi_abort(MPI_COMM_WORLD, 1, error)
    endif

    call h5pclose_f(p_id, error)

    file%current_group=0
    file%id = f_id

  end subroutine create_file

  ! Subroutine to read HDF5 file with collective IO
  subroutine read_file(filename,file)
    implicit none
    type(datafile), intent(out) :: file
    character(*), intent(in) :: filename
    integer(kind=hid_t) :: p_id, f_id

    ! MPI IO/Lustre file striping
    integer :: info         ! MPI IO Info
    integer :: lcount       ! Lustre count size
    integer :: lsize        ! Lustre stripe size
    character(len=1024) :: clcount, clsize ! Strings of LFS
    
    ! Initialize HDF5 environment
    call h5open_f(error)

    ! Set a stripe count of 4 and a stripe size of 4MB
    lcount = 4
    lsize  = 4 * 1024 * 1024
    write(clcount, '(I4)') lcount
    write(clsize, '(I8)') lsize

    call mpi_info_create(info, error)
    call mpi_info_set(info, "striping_factor", trim(clcount), error)
    call mpi_info_set(info, "striping_unit", trim(clsize), error)

    ! Set up the access properties
    call h5pcreate_f(H5P_FILE_ACCESS_F, p_id, error)
    call h5pset_fapl_mpio_f(p_id, MPI_COMM_WORLD, info, error)

    ! Open HDF5 file
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, f_id, error, &
                   access_prp = p_id)
    if (error .ne. 0) then
      write(0,*) 'Unable to open: ', trim(filename), ': ', error
      call mpi_abort(MPI_COMM_WORLD, 1, error)
    endif

    call h5pclose_f(p_id, error)

    file%current_group=0
    file%id = f_id

  end subroutine read_file

  ! Subroutine to close HDF5 file
  subroutine close_file(file)
    type(datafile), intent(in) :: file
    ! Close HDF5 file
    call h5fclose_f(file%id, error)
  end subroutine close_file

  ! Subroutine to write string attributes
  subroutine write_string(file,data,dname)
    implicit none
    character(*), intent(in) :: data
    character(*), intent(in) :: dname
    type(datafile), intent(in) :: file
    integer :: error
    integer(HID_T) :: filespace,memspace,dset_id,attr_id,group_id,dspace_id,dtype_id
    integer(hsize_t) :: dims(1)

    dims(1) = 1
    group_id = file%current_group

    ! Create scalar dataspace
    call h5screate_f(H5S_SCALAR_F, dspace_id, error)

    ! Create the datatype
    call h5tcopy_f(H5T_NATIVE_CHARACTER, dtype_id, error)
    call h5tset_size_f(dtype_id, int(len(trim(data)), kind=size_t), error)

    ! Create the attribute
    call h5acreate_f(group_id, dname,dtype_id , dspace_id, attr_id, error)
    call h5awrite_f(attr_id, dtype_id, data, dims, error)
    call h5aclose_f(attr_id, error)

    call h5sclose_f(dspace_id, error)
    call h5tclose_f(dtype_id, error)

  end subroutine write_string

  ! Subroutine to write integer attributes
  subroutine write_int(file,data,dname)
    implicit none
    integer(KIND=4), intent(in) :: data
    character(*), intent(in) :: dname
    type(datafile), intent(in) :: file
    integer :: error
    integer(HID_T) :: filespace,memspace,dset_id,attr_id,group_id,dspace_id
    integer(hsize_t) :: dims(1)

    dims(1) = 1
    group_id = file%current_group

    ! Create scalar dataspace
    call h5screate_f(H5S_SCALAR_F, dspace_id, error)

    ! Create the attribute
    call h5acreate_f(group_id, dname, H5T_NATIVE_INTEGER, dspace_id, attr_id, error)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, data, dims, error)
    call h5aclose_f(attr_id, error)

    call h5sclose_f(dspace_id, error)

  end subroutine write_int

  ! Subroutine to read integer attributes
  subroutine read_int(file,data,dname)
    implicit none
    integer, intent(out) :: data
    character(*), intent(in) :: dname
    type(datafile), intent(in) :: file
    integer :: error
    integer(HID_T) :: filespace,memspace,dset_id,attr_id,group_id,dspace_id
    integer(hsize_t) :: ndims(1)

    ndims(1) = 1
    group_id = file%current_group   

    call h5aopen_f(group_id, dname, attr_id, error)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, data,ndims, error)
    call h5aclose_f(attr_id, error)

  end subroutine read_int

  ! Subroutine to write real attributes
  subroutine write_real(file,data,dname)
    implicit none
    real(KIND=8), intent(in) :: data
    character(*), intent(in) :: dname
    type(datafile), intent(in) :: file
    integer :: error
    integer(HID_T) :: filespace,memspace,dset_id,attr_id,group_id,dspace_id
    integer(hsize_t) :: dims(1)

    dims(1) = 1
    group_id = file%current_group

    ! Create scalar dataspace
    call h5screate_f(H5S_SCALAR_F, dspace_id, error)

    ! Create the attribute
    call h5acreate_f(group_id, dname, H5T_NATIVE_DOUBLE, dspace_id, attr_id, error)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, data, dims, error)
    call h5aclose_f(attr_id, error)

    call h5sclose_f(dspace_id, error)

  end subroutine write_real

  ! Subroutine to read real attributes
  subroutine read_real(file,data,dname)
    implicit none
    double precision, intent(out) :: data
    character(*), intent(in) :: dname
    type(datafile), intent(in) :: file
    integer :: error
    integer(HID_T) :: filespace,memspace,dset_id,attr_id,group_id,dspace_id
    integer(hsize_t) :: dims(1)

    dims(1) = 1
    group_id = file%current_group   

    call h5aopen_f(group_id, dname, attr_id, error)
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, data,dims, error)
    call h5aclose_f(attr_id, error)

  end subroutine read_real

  ! Subroutine to write real 1d arrays
  subroutine write_1dreal(file,data,dname)
    implicit none
    integer, parameter :: rank=1
    double precision, intent(in) :: data(:)
    character(*), intent(in) :: dname
    type(datafile), intent(in) :: file

    integer(HSIZE_T) :: dims(rank)
    integer :: error

    integer(HID_T) :: filespace,d_id,group_id

    dims = shape(data)
    group_id = file%current_group

    call h5screate_simple_f(rank, dims, filespace, error)
    call h5dcreate_f(group_id, dname, H5T_NATIVE_DOUBLE, filespace, d_id, error)
    call h5dwrite_f(d_id, H5T_NATIVE_DOUBLE, data, dims, error)
    call h5dclose_f(d_id, error)
    call h5sclose_f(filespace,error)

  end subroutine write_1dreal

  ! Subroutine to read real 1d arrays
  subroutine read_1dreal(file,data,dname)
    implicit none
    integer, parameter :: rank=1
    double precision, intent(out) :: data(:)
    character(*), intent(in) :: dname
    type(datafile), intent(in) :: file

    integer(HSIZE_T) :: ndims(rank)
    integer :: error

    integer(HID_T) :: filespace,d_id,group_id

    ndims = shape(data)
    group_id = file%current_group

    call h5dopen_f(group_id, dname, d_id, error)
    call h5dread_f(d_id, H5T_NATIVE_DOUBLE, data, ndims, error)
    call h5dclose_f(d_id, error)

  end subroutine read_1dreal
  
    ! Subroutine to write real 3d arrays 
  subroutine write_3dreal(file,data,dimsf,dname,n1s,n1e,n2s,n2e,n3s,n3e)
    implicit none
    integer, parameter :: ndims=3
    real*8, intent(in) :: data(:,:,:)
    character(*), intent(in) :: dname
    type(datafile), intent(in) :: file
    integer(HSIZE_T), dimension(ndims), intent(in) :: dimsf
    integer, intent(in) :: n1s, n1e, n2s, n2e, n3s, n3e
    integer(kind=hid_t) :: x_id, d_id, group_id
    integer(kind=hid_t) :: memspace, filespace
    ! Local hyper slab info
    integer(kind=hsize_t) :: count(ndims), offset(ndims)

    group_id = file%current_group

    ! Create the data space for the  dataset.
    call h5screate_simple_f(ndims, dimsf, filespace, error)
        
    call h5dcreate_f(group_id, dname, H5T_IEEE_F64LE, filespace, d_id, &
                     error)

    ! Each process defines dataset in memory and writes it to the hyperslab
    ! in the file.
    count(1)  = n1e-n1s+1
    count(2)  = n2e-n2s+1
    count(3)  = n3e-n3s+1
    offset(1) = n1s-1
    offset(2) = n2s-1
    offset(3) = n3s-1

    call h5screate_simple_f(ndims, count, memspace, error)    

    ! Select hyperslab in the file.
    call h5dget_space_f(d_id, filespace, error)

    call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error)

    ! Create a data transfer property
    call h5pcreate_f(H5P_DATASET_XFER_F, x_id, error)
        
    call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, error)
        
    ! Write the data
    call h5dwrite_f(d_id, H5T_IEEE_F64LE, data, dimsf, error,    &
                    file_space_id=filespace, mem_space_id=memspace, &
                    xfer_prp=x_id)                  

    ! Close everything and exit
    call h5dclose_f(d_id, error)

        
    call h5sclose_f(filespace, error)
 
    
    call h5sclose_f(memspace, error)

        
    call h5pclose_f(x_id, error) 

  end subroutine write_3dreal
 

  ! Subroutine to read real 3d arrays 
  subroutine read_3dreal(file,data,dimsf,dname,n1s,n1e,n2s,n2e,n3s,n3e)
    implicit none
    integer, parameter :: ndims=3
    real*8, intent(out) :: data(:,:,:)
    character(*), intent(in) :: dname
    type(datafile), intent(in) :: file
    integer(HSIZE_T), dimension(ndims), intent(in) :: dimsf
    integer, intent(in) :: n1s, n1e, n2s, n2e, n3s, n3e
    integer(kind=hid_t) :: x_id, d_id, group_id
    integer(kind=hid_t) :: memspace, filespace
    ! Local hyper slab info
    integer(kind=hsize_t) :: count(ndims), offset(ndims)

    group_id = file%current_group

    ! Create the data space for the  dataset.
    call h5screate_simple_f(ndims, dimsf, filespace, error)      
        
    call h5dopen_f(group_id, dname, d_id, error)
    
    ! Each process defines dataset in memory and writes it to the hyperslab
    ! in the file.
    count(1)  = n1e-n1s+1
    count(2)  = n2e-n2s+1
    count(3)  = n3e-n3s+1
    offset(1) = n1s-1
    offset(2) = n2s-1
    offset(3) = n3s-1

    call h5screate_simple_f(ndims, count, memspace, error)
    
    ! Select hyperslab in the file.
    call h5dget_space_f(d_id, filespace, error)
    
    call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error) 
    
    ! Create a data transfer property
    call h5pcreate_f(H5P_DATASET_XFER_F, x_id, error)
    
    call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, error)       
    
    ! Write the data
    call h5dread_f(d_id, H5T_IEEE_F64LE, data, dimsf, error,    &
                   file_space_id=filespace, mem_space_id=memspace, &
                   xfer_prp=x_id)                 

    ! Close everything and exit
    call h5dclose_f(d_id, error)     
    
    call h5sclose_f(filespace, error)   
    
    call h5sclose_f(memspace, error)
   
    call h5pclose_f(x_id, error) 


  end subroutine read_3dreal

  ! Subroutine to write XDMF file for visualization in Paraview or Visit
  subroutine write_xdmf(path,filename,time,dimsf,SCA_state)
    integer, parameter :: lun=42
    character(*), intent(in) :: filename, path
    double precision, intent(in) :: time
    integer(kind=hsize_t), intent(in) :: dimsf(3)
    integer :: SCA_state

    open(lun,file=trim(path)//trim(filename)//'.xmf',action="write")
    write(lun,'(A)') '<?xml version="1.0" ?>'
    write(lun,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(lun,'(A)') '<Xdmf Version="2.0">'
    write(lun,'(A)') '<Domain>'
    write(lun,'(A)') '<Grid Name="mesh" GridType="Uniform">'
    write(lun,'(A,3I6,A)') '<Topology TopologyType="3DRectMesh" Dimensions="',dimsf(1:3),'"/>'

    write(lun,'(A)') '<Geometry GeometryType="VXVYVZ">'

    write(lun,'(A,I6,A)') ' <DataItem Dimensions="',dimsf(3),'" Name="Xc" NumberType="Float" Precision="8" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/grid/Zc'
    write(lun,'(A)') ' </DataItem>'

    write(lun,'(A,I6,A)') ' <DataItem Dimensions="',dimsf(2),'" Name="Yc" NumberType="Float" Precision="8" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/grid/Yc'
    write(lun,'(A)') ' </DataItem>'

    write(lun,'(A,I6,A)') ' <DataItem Dimensions="',dimsf(1),'" Name="Zc" NumberType="Float" Precision="8" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/grid/Xc'
    write(lun,'(A)') ' </DataItem>'

    write(lun,'(A)') '</Geometry>'

    write(lun,'(A,1E11.4,A)') '<Time Value="',time,'" />'

    write(lun,'(A)') '<Attribute Name="pressure" AttributeType="Scalar" Center="Node">'
    write(lun,'(A,3I6,A)') '  <DataItem Dimensions="',dimsf(1:3),'" NumberType="Float" Precision="8" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/fields/pressure'
    write(lun,'(A)') ' </DataItem>'
    write(lun,'(A)') '</Attribute>'

    if (SCA_state==1) then
        write(lun,'(A)') '<Attribute Name="temperature" AttributeType="Scalar" Center="Node">'
        write(lun,'(A,3I6,A)') '  <DataItem Dimensions="',dimsf(1:3),'" NumberType="Float" Precision="8" Format="HDF">'
        write(lun,'(3A)') '  ',filename//'.h5',':/fields/temperature'
        write(lun,'(A)') ' </DataItem>'
        write(lun,'(A)') '</Attribute>'
    endif

    write(lun,'(A)') '<Attribute Name="U" AttributeType="Scalar" Center="Node">'

    write(lun,'(A,3I6,A)') '  <DataItem Dimensions="',dimsf(1:3),'" NumberType="Float" Precision="8" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/fields/velocity/U'
    write(lun,'(A)') '  </DataItem>'
    write(lun,'(A)') '</Attribute>'

    write(lun,'(A)') '<Attribute Name="V" AttributeType="Scalar" Center="Node">'

    write(lun,'(A,3I6,A)') '  <DataItem Dimensions="',dimsf(1:3),'" NumberType="Float" Precision="8" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/fields/velocity/V'
    write(lun,'(A)') '  </DataItem>'
    write(lun,'(A)') '</Attribute>'

    write(lun,'(A)') '<Attribute Name="W" AttributeType="Scalar" Center="Node">'

    write(lun,'(A,3I6,A)') '  <DataItem Dimensions="',dimsf(1:3),'" NumberType="Float" Precision="8" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/fields/velocity/W'
    write(lun,'(A)') '  </DataItem>'
    write(lun,'(A)') '</Attribute>'


    write(lun,'(A)') '<Attribute Name="velocity" AttributeType="Vector" Center="Node">'
    write(lun,'(A,4I6,A)') '<DataItem ItemType="Function" Dimensions="',dimsf(1:3),3,'" Function="JOIN($0 , $1, $2)">'

    write(lun,'(A,3I6,A)') '  <DataItem Dimensions="',dimsf(1:3),'" NumberType="Float" Precision="8" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/fields/velocity/U'
    write(lun,'(A)') '  </DataItem>'

    write(lun,'(A,3I6,A)') '  <DataItem Dimensions="',dimsf(1:3),'" NumberType="Float" Precision="8" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/fields/velocity/V'
    write(lun,'(A)') '  </DataItem>'

    write(lun,'(A,3I6,A)') '  <DataItem Dimensions="',dimsf(1:3),'" NumberType="Float" Precision="8" Format="HDF">'
    write(lun,'(3A)') '  ',filename//'.h5',':/fields/velocity/W'
    write(lun,'(A)') '  </DataItem>'

    write(lun,'(A)') '</DataItem>'
    write(lun,'(A)') '</Attribute>'


    write(lun,'(A)') '</Grid>'
    write(lun,'(A)') '</Domain>'
    write(lun,'(A)') '</Xdmf>'
    close(lun)

  end subroutine write_xdmf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine hdf_create_emptyfile(filename)

        implicit none
        ! Dataset names
        character(*), intent(in)        :: filename

        integer(HID_T)                  :: file_id                          ! File identifier

        integer                         :: error

        call h5open_f(error)

        ! Create the file collectively.
        call h5fcreate_f(trim(filename)//".h5", H5F_ACC_TRUNC_F, file_id, error)
        call h5fclose_f(file_id, error)

        call h5close_f(error)

    end subroutine hdf_create_emptyfile


    subroutine hdf_create_file_with_2Dmesh(filename, x,y, x_name, y_name, n1,n2)

        implicit none
        ! Dataset names
        character(*), intent(in)        :: filename, x_name, y_name  ! File name
        integer, intent(in)             :: n1, n2

        integer(HID_T)                  :: file_id                  ! File identifier
        integer(HID_T)                  :: dsetx_id, dsety_id       ! Dataset identifier
        integer(HID_T)                  :: spacex_id, spacey_id     ! Dataspace identifier x
        integer(HSIZE_T), dimension(1)  :: dimsx, dimsy             ! Dataset dimensions.
        real*8, dimension(n1)           :: x
        real*8, dimension(n2)           :: y

        integer                         :: error

        dimsx = (/n1/)
        dimsy = (/n2/)

        call h5open_f(error)

        ! Create the file collectively.
        call h5fcreate_f(trim(filename)//".h5", H5F_ACC_TRUNC_F, file_id, error)

        ! Create the data space for the  dataset.
        call h5screate_simple_f(1, dimsx, spacex_id, error)
        call h5screate_simple_f(1, dimsy, spacey_id, error)

        ! Create the dataset with default properties.
        call h5dcreate_f(file_id, x_name, H5T_NATIVE_DOUBLE, spacex_id, dsetx_id, error)
        call h5dcreate_f(file_id, y_name, H5T_NATIVE_DOUBLE, spacey_id, dsety_id, error)


        ! Write the dataset independently.
        call h5dwrite_f(dsetx_id, H5T_NATIVE_DOUBLE, x, dimsx, error)
        call h5dwrite_f(dsety_id, H5T_NATIVE_DOUBLE, y, dimsy, error)


        ! Close dataspaces.
        call h5sclose_f(spacex_id, error)
        call h5sclose_f(spacey_id, error)


        ! Close the dataset and property list.
        call h5dclose_f(dsetx_id, error)
        call h5dclose_f(dsety_id, error)


        ! Close the file.
        call h5fclose_f(file_id, error)

        call h5close_f(error)

    end subroutine hdf_create_file_with_2Dmesh

    subroutine hdf_addgroup(filename, groupname)

        implicit none
        ! Dataset names
        character(*), intent(in)        :: filename, groupname

        integer(HID_T)                  :: file_id, group_id                          ! File identifier

        integer                         :: error

        call h5open_f(error)

        call h5fopen_f (trim(filename)//".h5", H5F_ACC_RDWR_F, file_id, error)
        call h5gcreate_f(file_id, trim(groupname), group_id, error)
        call h5gclose_f(group_id, error)
        call h5fclose_f(file_id, error)

        call h5close_f(error)

    end subroutine hdf_addgroup


    subroutine hdf_create_file_with_3Dmesh(filename, x,y,z, x_name, y_name, z_name, n1,n2,n3)

        implicit none
        ! Dataset names
        character(*), intent(in)        :: filename, x_name, y_name, z_name  ! File name
        integer, intent(in)             :: n1, n2, n3

        integer(HID_T)                  :: file_id                          ! File identifier
        integer(HID_T)                  :: dsetx_id, dsety_id, dsetz_id     ! Dataset identifier
        integer(HID_T)                  :: spacex_id, spacey_id, spacez_id  ! Dataspace identifier x
        integer(HSIZE_T), dimension(1)  :: dimsx, dimsy, dimsz              ! Dataset dimensions.
        real*8, dimension(n1)           :: x
        real*8, dimension(n2)           :: y
        real*8, dimension(n3)           :: z

        integer                         :: error

        dimsx = (/n1/)
        dimsy = (/n2/)
        dimsz = (/n3/)

        call h5open_f(error)

        ! Create the file collectively.
        call h5fcreate_f(trim(filename)//".h5", H5F_ACC_TRUNC_F, file_id, error)

        ! Create the data space for the  dataset.
        call h5screate_simple_f(1, dimsx, spacex_id, error)
        call h5screate_simple_f(1, dimsy, spacey_id, error)
        call h5screate_simple_f(1, dimsz, spacez_id, error)

        ! Create the dataset with default properties.
        call h5dcreate_f(file_id, x_name, H5T_NATIVE_DOUBLE, spacex_id, dsetx_id, error)
        call h5dcreate_f(file_id, y_name, H5T_NATIVE_DOUBLE, spacey_id, dsety_id, error)
        call h5dcreate_f(file_id, z_name, H5T_NATIVE_DOUBLE, spacez_id, dsetz_id, error)


        ! Write the dataset independently.
        call h5dwrite_f(dsetx_id, H5T_NATIVE_DOUBLE, x, dimsx, error)
        call h5dwrite_f(dsety_id, H5T_NATIVE_DOUBLE, y, dimsy, error)
        call h5dwrite_f(dsetz_id, H5T_NATIVE_DOUBLE, z, dimsz, error)


        ! Close dataspaces.
        call h5sclose_f(spacex_id, error)
        call h5sclose_f(spacey_id, error)
        call h5sclose_f(spacez_id, error)


        ! Close the dataset and property list.
        call h5dclose_f(dsetx_id, error)
        call h5dclose_f(dsety_id, error)
        call h5dclose_f(dsetz_id, error)


        ! Close the file.
        call h5fclose_f(file_id, error)

        call h5close_f(error)

    end subroutine hdf_create_file_with_3Dmesh
    !
    subroutine hdf_add_2Dfield(filename, field, field_name, n1, n2, n1s, n1e, n2s, n2e)

        implicit none

        character(*), intent(in)        :: filename, field_name  ! File name
        integer, intent(in)             :: n1, n2, n1s, n1e, n2s, n2e
        real*8, dimension(n1s:n1e,n2s:n2e), intent(in)        :: field

        integer(HID_T) :: file_id       ! File identifier
        integer(HID_T) :: dset_id       ! Dataset identifier
        integer(HID_T) :: filespace     ! Dataspace identifier in file
        integer(HID_T) :: memspace      ! Dataspace identifier in memory
        integer(HID_T) :: plist_id, plist_id2      ! Property list identifier

        integer, parameter                  :: rank = 2     ! Dataset rank
        integer(HSIZE_T), dimension(rank)   :: dimsf        ! Dataset dimensions.

        integer(HSIZE_T), dimension(rank)   :: count
        integer(HSSIZE_T), dimension(rank)  :: offset

        integer :: error  ! Error flags

        dimsf = (/n1,n2/)

        call h5open_f(error)


        ! ***************************************************************************************************
        ! ****************************** OPEN FILE *********************************************************
        ! ***************************************************************************************************


        ! Setup file access property list with parallel I/O access.
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

        ! Create the file collectively.
        call h5fopen_f(trim(filename)//".h5", H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
        call h5pclose_f(plist_id, error)

        ! Create the data space for the  dataset.
        call h5screate_simple_f(rank, dimsf, filespace, error)

        ! Create the dataset with default properties.
        call h5dcreate_f(file_id, field_name, H5T_NATIVE_DOUBLE, filespace, dset_id, error)
        call h5sclose_f(filespace, error)

        ! ***************************************************************************************************
        ! ****************************** WRITE DATA IN FILE ************************************************
        ! ***************************************************************************************************

        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        count(1) = min(n1e,n1)-n1s+1 !ysize(1)
        count(2) = min(n2e,n2)-n2s+1 !ysize(2)
        offset(1) = n1s-1
        offset(2) = n2s-1

        call h5screate_simple_f(rank, count, memspace, error)

        ! Select hyperslab in the file.
        call h5dget_space_f(dset_id, filespace, error)
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error)



        ! Create property list for collective dataset write
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ! Write the dataset collectively.
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field, dimsf, error, &
        file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)


        ! ***************************************************************************************************
        ! ****************************** CLOSE HDF5 OBJECTS *************************************************
        ! ***************************************************************************************************

        ! Close dataspaces.
        call h5sclose_f(filespace, error)
        call h5sclose_f(memspace, error)


        ! Close the dataset and property list.
        call h5dclose_f(dset_id, error)
        call h5pclose_f(plist_id, error)


        ! Close the file.
        call h5fclose_f(file_id, error)


        call h5close_f(error)

    end subroutine
    !
    subroutine hdf_read_2Dfield(filename, field, field_name, n1, n2, n1s, n1e, n2s, n2e)

        implicit none

        character(*), intent(in)        :: filename, field_name  ! File name
        integer, intent(in)             :: n1, n2, n1s, n1e, n2s, n2e
        real*8, dimension(n1s:n1e,n2s:n2e), intent(inout)        :: field

        integer(HID_T) :: file_id       ! File identifier
        integer(HID_T) :: dset_id       ! Dataset identifier
        integer(HID_T) :: filespace     ! Dataspace identifier in file
        integer(HID_T) :: memspace      ! Dataspace identifier in memory
        integer(HID_T) :: plist_id, plist_id2      ! Property list identifier

        integer, parameter                  :: rank = 2 ! Dataset rank
        integer(HSIZE_T), dimension(rank)   :: dimsf ! Dataset dimensions.

        integer(HSIZE_T), dimension(rank)   :: count
        integer(HSSIZE_T), dimension(rank)  :: offset

        integer :: error  ! Error flags

        dimsf = (/n1,n2/)

        call h5open_f(error)


        ! ***************************************************************************************************
        ! ****************************** OPEN FILE *********************************************************
        ! ***************************************************************************************************


        ! Setup file access property list with parallel I/O access.
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

        ! Create the file collectively.
        call h5fopen_f(trim(filename)//".h5", H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
        call h5pclose_f(plist_id, error)

        ! Open the dataset
        call h5dopen_f(file_id, field_name, dset_id, error)

        ! ***************************************************************************************************
        ! ****************************** WRITE DATA IN FILE ************************************************
        ! ***************************************************************************************************

        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        count(1) = min(n1e,n1)-n1s+1 !ysize(1)
        count(2) = min(n2e,n2)-n2s+1 !ysize(2)
        offset(1) = n1s-1
        offset(2) = n2s-1

        call h5screate_simple_f(rank, count, memspace, error)

        ! Select hyperslab in the file.
        call h5dget_space_f(dset_id, filespace, error)
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error)



        ! Create property list for collective dataset write
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ! Read the dataset collectively.
        CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, field, count, error, &
        memspace, filespace)



        ! ***************************************************************************************************
        ! ****************************** CLOSE HDF5 OBJECTS *************************************************
        ! ***************************************************************************************************

        ! Close dataspaces.
        call h5sclose_f(filespace, error)
        call h5sclose_f(memspace, error)


        ! Close the dataset and property list.
        call h5dclose_f(dset_id, error)
        call h5pclose_f(plist_id, error)


        ! Close the file.
        call h5fclose_f(file_id, error)


        call h5close_f(error)

    end subroutine


    subroutine hdf_add_3Dfield(filename, field, field_name, n1, n2, n3, n1s, n1e, n2s, n2e, n3s, n3e)

        use MPI
        USE HDF5 ! This module contains all necessary modules

        implicit none

        character(*), intent(in)        :: filename, field_name  ! File name
        integer, intent(in)             :: n1, n2, n3, n1s, n1e, n2s, n2e, n3s, n3e
        real*8, dimension(n1s:n1e,n2s:n2e,n3s:n3e), intent(in)        :: field

        integer(HID_T) :: file_id       ! File identifier
        integer(HID_T) :: dset_id       ! Dataset identifier
        integer(HID_T) :: filespace     ! Dataspace identifier in file
        integer(HID_T) :: memspace      ! Dataspace identifier in memory
        integer(HID_T) :: plist_id, plist_id2      ! Property list identifier

        integer, parameter                  :: rank = 3 ! Dataset rank
        integer(HSIZE_T), dimension(rank)   :: dimsf ! Dataset dimensions.

        integer(HSIZE_T), dimension(rank)   :: count
        integer(HSSIZE_T), dimension(rank)  :: offset

        integer :: error  ! Error flags

        dimsf = (/n1,n2,n3/)

        call h5open_f(error)


        ! ***************************************************************************************************
        ! ****************************** OPEN FILE *********************************************************
        ! ***************************************************************************************************


        ! Setup file access property list with parallel I/O access.
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

        ! Create the file collectively.
        call h5fopen_f(trim(filename)//".h5", H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
        call h5pclose_f(plist_id, error)

        ! Create the data space for the  dataset.
        call h5screate_simple_f(rank, dimsf, filespace, error)

        ! Create the dataset with default properties.
        call h5dcreate_f(file_id, field_name, H5T_NATIVE_DOUBLE, filespace, dset_id, error)
        call h5sclose_f(filespace, error)

        ! ***************************************************************************************************
        ! ****************************** WRITE DATA IN FILE ************************************************
        ! ***************************************************************************************************

        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        count(1) = min(n1e,n1)-n1s+1 !ysize(1)
        count(2) = min(n2e,n2)-n2s+1 !ysize(2)
        count(3) = min(n3e,n3)-n3s+1 !ysize(3)
        offset(1) = n1s-1
        offset(2) = n2s-1
        offset(3) = n3s-1

        call h5screate_simple_f(rank, count, memspace, error)

        ! Select hyperslab in the file.
        call h5dget_space_f(dset_id, filespace, error)
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error)



        ! Create property list for collective dataset write
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ! Write the dataset collectively.
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field, dimsf, error, &
        file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)


        ! ***************************************************************************************************
        ! ****************************** CLOSE HDF5 OBJECTS *************************************************
        ! ***************************************************************************************************

        ! Close dataspaces.
        call h5sclose_f(filespace, error)
        call h5sclose_f(memspace, error)


        ! Close the dataset and property list.
        call h5dclose_f(dset_id, error)
        call h5pclose_f(plist_id, error)


        ! Close the file.
        call h5fclose_f(file_id, error)


        call h5close_f(error)

    end subroutine


    subroutine hdf_write_3Dfield(filepath, field, field_name, n1, n2, n3, n1s, n1e, n2s, n2e, n3s, n3e)

        use MPI
        USE HDF5 ! This module contains all necessary modules

        implicit none

        character(*), intent(in)                                        :: filepath, field_name  ! File name
        integer, intent(in)                                             :: n1, n2, n3, n1s, n1e, n2s, n2e, n3s, n3e
        real*8, dimension(n1s:n1e,n2s:n2e,n3s:n3e), intent(in)          :: field


        integer :: error  ! Error flags
        integer :: rank


        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, error)

        if(rank==0)  call hdf_create_emptyfile(trim(filepath))
        call hdf_add_3Dfield(trim(filepath), field(:,:,:), trim(field_name), n1, n2, n3, n1s, n1e, n2s, n2e, n3s, n3e)

    end subroutine
    !
    subroutine hdf_read_3Dfield(filename, field, field_name, n1, n2, n3, n1s, n1e, n2s, n2e, n3s, n3e)

        implicit none

        character(*), intent(in)        :: filename, field_name  ! File name
        integer, intent(in)             :: n1, n2, n3, n1s, n1e, n2s, n2e, n3s, n3e
        real*8, dimension(n1s:n1e,n2s:n2e,n3s:n3e), intent(inout)        :: field

        integer(HID_T) :: file_id       ! File identifier
        integer(HID_T) :: dset_id       ! Dataset identifier
        integer(HID_T) :: filespace     ! Dataspace identifier in file
        integer(HID_T) :: memspace      ! Dataspace identifier in memory
        integer(HID_T) :: plist_id, plist_id2      ! Property list identifier

        integer, parameter                  :: rank = 3 ! Dataset rank
        integer(HSIZE_T), dimension(rank)   :: dimsf ! Dataset dimensions.

        integer(HSIZE_T), dimension(rank)   :: count
        integer(HSSIZE_T), dimension(rank)  :: offset

        integer :: error  ! Error flags

        dimsf = (/n1,n2,n3/)

        call h5open_f(error)


        ! ***************************************************************************************************
        ! ****************************** OPEN FILE *********************************************************
        ! ***************************************************************************************************


        ! Setup file access property list with parallel I/O access.
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

        ! Create the file collectively.
        call h5fopen_f(trim(filename)//".h5", H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
        call h5pclose_f(plist_id, error)

        ! Open the dataset
        call h5dopen_f(file_id, field_name, dset_id, error)

        ! ***************************************************************************************************
        ! ****************************** WRITE DATA IN FILE ************************************************
        ! ***************************************************************************************************

        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        count(1) = min(n1e,n1)-n1s+1 !ysize(1)
        count(2) = min(n2e,n2)-n2s+1 !ysize(2)
        count(3) = min(n3e,n3)-n3s+1 !ysize(3)
        offset(1) = n1s-1
        offset(2) = n2s-1
        offset(3) = n3s-1

        call h5screate_simple_f(rank, count, memspace, error)

        ! Select hyperslab in the file.
        call h5dget_space_f(dset_id, filespace, error)
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error)



        ! Create property list for collective dataset write
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ! Read the dataset collectively.
        CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, field, count, error, &
        memspace, filespace)



        ! ***************************************************************************************************
        ! ****************************** CLOSE HDF5 OBJECTS *************************************************
        ! ***************************************************************************************************

        ! Close dataspaces.
        call h5sclose_f(filespace, error)
        call h5sclose_f(memspace, error)


        ! Close the dataset and property list.
        call h5dclose_f(dset_id, error)
        call h5pclose_f(plist_id, error)


        ! Close the file.
        call h5fclose_f(file_id, error)


        call h5close_f(error)

    end subroutine

end module HDF5_IO
