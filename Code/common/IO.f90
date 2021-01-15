module time_writer
    use run_ctxt_data
    use decomp_2d
    implicit none

contains

    subroutine read_timefile(filepath, field_it, field_time)
        implicit none
        character(*)    :: filepath
        integer         :: field_it
        real*8          :: field_time

        open(15,file=trim(filepath))
        read(15,*)field_it
        read(15,*)field_time
        close(15)

    end subroutine read_timefile

    subroutine write_timefile(filepath)

        implicit none
        character(*)    :: filepath

        ! Save the current time in the file advancement.d
        if (nrank==0) then
           write(*,*)"write_timefile directory", filepath

           open(15,file=trim(filepath))

            write(15,*)ntime
            write(15,*)t
            close(15)
        end if

    end subroutine write_timefile

end module time_writer

module snapshot_writer

    use mpi
    use decomp_2d

    implicit none

contains

    subroutine create_snapshot(snaps_dir, snap_dir, field, field_name, pencil)

        use mesh
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
        !if (pencil==1) call hdf_add_3Dfield(file_path, field(:,:,:), trim(field_name), nx_global, ny_global, nz_global, 1,xsize(1),1,xsize(2),1,xsize(3))
        if (pencil==1) call hdf_add_3Dfield(file_path, field(:,:,:), trim(field_name), nx_global, ny_global, nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
        !if (pencil==2) call hdf_add_3Dfield(file_path, field(:,:,:), trim(field_name), nx_global, ny_global, nz_global, 1,ysize(1),1,ysize(2),1,ysize(3))
        if (pencil==2) call hdf_add_3Dfield(file_path, field(:,:,:), trim(field_name), nx_global, ny_global, nz_global, ystart(1),yend(1),ystart(2),yend(2),ystart(3),yend(3))
        !if (pencil==3) call hdf_add_3Dfield(file_path, field(:,:,:), trim(field_name), nx_global, ny_global, nz_global, 1,zsize(1),1,zsize(2),1,zsize(3))
        if (pencil==3) call hdf_add_3Dfield(file_path, field(:,:,:), trim(field_name), nx_global, ny_global, nz_global, zstart(1),zend(1),zstart(2),zend(2),zstart(3),zend(3))

        call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

    end subroutine create_snapshot



    subroutine create_stretch_snapshot(snaps_dir, snap_dir, field, field_name, pencil, X1, X2, X3)

        use mesh
        use HDF5_IO

        implicit none
        real*8, dimension(:,:,:)    :: field
        character(*)                :: snaps_dir, snap_dir, field_name
        integer                     :: pencil

        integer                     :: mpi_err


        character(200)    :: file_path, snap_path
        logical             :: snap_dir_exist

        real*8                              :: X1(n1), X2(n2), X3(n3)

        snap_path=trim(snaps_dir)//"/"//trim(snap_dir)

        file_path=trim(snap_path)//"/"//trim(field_name)

        if(nrank==0)  call hdf_create_file_with_3Dmesh(file_path, X3, X2, X1, "Xaxis", "Yaxis", "Zaxis", n3,n2,n1)
        if (pencil==1) call hdf_add_3Dfield(file_path, field(:,:,:), trim(field_name), nx_global, ny_global, nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
        if (pencil==2) call hdf_add_3Dfield(file_path, field(:,:,:), trim(field_name), nx_global, ny_global, nz_global, ystart(1),yend(1),ystart(2),yend(2),ystart(3),yend(3))
        if (pencil==3) call hdf_add_3Dfield(file_path, field(:,:,:), trim(field_name), nx_global, ny_global, nz_global, zstart(1),zend(1),zstart(2),zend(2),zstart(3),zend(3))

        call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

    end subroutine create_stretch_snapshot



    subroutine create_2D_snapshot(snaps_dir, snap_dir, field, field_name, pencil)

        use HDF5_IO

        implicit none
        real*8, dimension(:,:)    :: field
        character(*)                :: snaps_dir, snap_dir, field_name
        integer                     :: pencil

        integer                     :: mpi_err

        character(200)    :: file_path, snap_path
        logical             :: snap_dir_exist


        snap_path=trim(snaps_dir)//"/"//trim(snap_dir)

        file_path=trim(snap_path)//"/"//trim(field_name)

        if(nrank==0)  call hdf_create_emptyfile(file_path)
        !if (pencil==1) call hdf_add_3Dfield(file_path, field(:,:), trim(field_name), nx_global, ny_global, nz_global, 1,xsize(1),1,xsize(2),1,xsize(3))
        if (pencil==1) call hdf_add_2Dfield(file_path, field(:,:), trim(field_name), ny_global, nz_global, xstart(2),xend(2),xstart(3),xend(3))
        !if (pencil==2) call hdf_add_3Dfield(file_path, field(:,:), trim(field_name), nx_global, ny_global, nz_global, 1,ysize(1),1,ysize(2),1,ysize(3))
        if (pencil==2) call hdf_add_2Dfield(file_path, field(:,:), trim(field_name), nx_global, nz_global, ystart(1),yend(1),ystart(3),yend(3))
        !if (pencil==3) call hdf_add_3Dfield(file_path, field(:,:), trim(field_name), nx_global, ny_global, nz_global, 1,zsize(1),1,zsize(2),1,zsize(3))
        if (pencil==3) call hdf_add_2Dfield(file_path, field(:,:), trim(field_name), nx_global, ny_global, zstart(1),zend(1),zstart(2),zend(2))

        call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

    end subroutine create_2D_snapshot

end module snapshot_writer
