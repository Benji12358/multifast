

module HDF5_IO
    USE HDF5
    use MPI
    implicit none

contains


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
