module mod_hdf5IO

  use hdf5
  use mod_myMPI

  implicit none

  private

  public :: create_file,read_file,close_file,write_hdf,read_hdf,datafile

  type :: datafile
      integer(kind=hid_t) :: id
      integer(HID_T) :: current_group
  end type datafile

  interface write_hdf
    module procedure write_int,write_real,write_1dreal,write_2dreal,write_3dreal,write_4dreal
    module procedure write_2dint,write_3dint,write_4dint
  end interface write_hdf

  interface read_hdf
    module procedure read_int,read_real,read_1dreal,read_2dreal,read_3dreal,read_4dreal
    module procedure read_2dint,read_3dint,read_4dint
  end interface read_hdf

    contains

      subroutine create_file(filename,file)
        implicit none
        type(datafile), intent(out) :: file
        character(*), intent(in) :: filename
        integer(kind=hid_t) :: p_id, f_id
    
        call myMPI()

        ! Initialize HDF5 environment
        call h5open_f(ierr)
         if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: create_file'
            write(*,*) 'h5open_f ierr:',ierr
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 

        ! Set up the access properties
        call h5pcreate_f(H5P_FILE_ACCESS_F, p_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: create_file'
            write(*,*) 'h5create_f ierr:',ierr
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if    
        call h5pset_fapl_mpio_f(p_id, comm, mpi_info_null, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: create_file'
            write(*,*) 'h5pset_fapl_mpio_f ierr:',ierr
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
    
        ! Create and Open HDF5 file
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, f_id, ierr, &
                         access_prp = p_id)
        if(ierr/=0) then
          write(*,*) 'Unable to create: ', trim(filename), ': ', ierr
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: create_file'
            write(*,*) 'h5fcreate_f ierr:',ierr
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if  

        call h5pclose_f(p_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: create_file'
            write(*,*) 'h5pclose_f ierr:',ierr
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if

        file%current_group=0
        file%id = f_id

      end subroutine create_file
  
      ! Subroutine to read HDF5 file with collective IO
      subroutine read_file(filename,file)
        implicit none
        type(datafile), intent(out) :: file
        character(*), intent(in) :: filename
        integer(kind=hid_t) :: p_id, f_id
        
        call myMPI()
        
        ! Initialize HDF5 environment
        call h5open_f(ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_file'
            write(*,*) 'h5open_f ierr:',ierr
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if    
        
        ! Set up the access properties
        call h5pcreate_f(H5P_FILE_ACCESS_F, p_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_file'
            write(*,*) 'h5pcreate_f ierr:',ierr
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if    
        
        call h5pset_fapl_mpio_f(p_id, comm, mpi_info_null, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_file'
            write(*,*) 'h5pset_fapl_mpio_f ierr:',ierr
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if  
        
        ! Open HDF5 file
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, f_id, ierr, &
                       access_prp = p_id)
        if(ierr/=0) then
          write(*,*) 'Unable to open: ', trim(filename), ': ', ierr
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_file'
            write(*,*) 'h5fopen_f ierr:',ierr
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if  
                      
        call h5pclose_f(p_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_file'
            write(*,*) 'h5pclose_f ierr:',ierr
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 

        file%current_group=0
        file%id = f_id

      end subroutine read_file  

      ! Subroutine to close HDF5 file
      subroutine close_file(file)
        type(datafile), intent(in) :: file
        
        call myMPI()
        
        ! Close HDF5 file
        call h5fclose_f(file%id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: close_file'
            write(*,*) 'h5fclose_f ierr:',ierr
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        ! Close HDF5 environment
        call h5close_f(ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: close_file'
            write(*,*) 'h5close_f ierr:',ierr
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
      end subroutine close_file
  
      ! Subroutine to write integer attributes
      subroutine write_int(file,data,dname)
        implicit none
        integer(KIND=4), intent(in) :: data
        character(*), intent(in) :: dname
        type(datafile), intent(in) :: file
        integer(HID_T) :: filespace,memspace,dset_id,attr_id,group_id,dspace_id
        integer(hsize_t) :: dims(1)

        call myMPI()

        dims(1) = 1
        group_id = file%current_group

        ! Create scalar dataspace
        call h5screate_f(H5S_SCALAR_F, dspace_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_int'
            write(*,*) 'h5screate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'dspace_id:', dspace_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if

        ! Create the attribute
        call h5acreate_f(group_id, dname, H5T_NATIVE_INTEGER, dspace_id, attr_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_int'
            write(*,*) 'h5acreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'dspace_id:', dspace_id
            write(*,*) 'attr_id:', attr_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
        
        call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, data, dims, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_int'
            write(*,*) 'h5awrite_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'dims:', dims
            write(*,*) 'attr_id:', attr_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if

        call h5aclose_f(attr_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_int'
            write(*,*) 'h5aclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'attr_id:', attr_id        
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
        
        call h5sclose_f(dspace_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_int'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'dspace_id:', dspace_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if

      end subroutine write_int
      
     ! Subroutine to read integer attributes
      subroutine read_int(file,data,dname)
        implicit none
        integer, intent(out) :: data
        character(*), intent(in) :: dname
        type(datafile), intent(in) :: file
        integer(HID_T) :: filespace,memspace,dset_id,attr_id,group_id,dspace_id
        integer(hsize_t) :: ndims(1)

        call myMPI()

        ndims(1) = 1
        group_id = file%current_group   

        call h5aopen_f(group_id, dname, attr_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_int'
            write(*,*) 'h5aopen_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'group_id:', group_id
            write(*,*) 'attr_id:', attr_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if    
        
        call h5aread_f(attr_id, H5T_NATIVE_INTEGER, data, ndims, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_int'
            write(*,*) 'h5aread_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'attr_id:', attr_id
            write(*,*) 'ndims:', ndims
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if  

        call h5aclose_f(attr_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_int'
            write(*,*) 'h5aclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'attr_id:', attr_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 

      end subroutine read_int  
      
      
      ! Subroutine to write real attributes
      subroutine write_real(file,data,dname)
        implicit none
        real*8, intent(in) :: data
        character(*), intent(in) :: dname
        type(datafile), intent(in) :: file
        integer(HID_T) :: filespace,memspace,dset_id,attr_id,group_id,dspace_id
        integer(hsize_t) :: dims(1)

        dims(1) = 1
        group_id = file%current_group

        ! Create scalar dataspace
        call h5screate_f(H5S_SCALAR_F, dspace_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_real'
            write(*,*) 'h5screate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'dspace_id:', dspace_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 

        ! Create the attribute
        call h5acreate_f(group_id, dname, H5T_NATIVE_DOUBLE, dspace_id, attr_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_real'
            write(*,*) 'h5acreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'group_id:', group_id
            write(*,*) 'dspace_id:', dspace_id
            write(*,*) 'attr_id:', attr_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
        call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, data, dims, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_real'
            write(*,*) 'h5awrite_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'attr_id:', attr_id
            write(*,*) 'dims:', dims
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
            
        call h5aclose_f(attr_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_real'
            write(*,*) 'h5aclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'attr_id:', attr_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 

        call h5sclose_f(dspace_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_real'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'dspace_id:', dspace_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 

      end subroutine write_real

      ! Subroutine to read real attributes
      subroutine read_real(file,data,dname)
        implicit none
        real*8, intent(out) :: data
        character(*), intent(in) :: dname
        type(datafile), intent(in) :: file
        integer(HID_T) :: filespace,memspace,dset_id,attr_id,group_id,dspace_id
        integer(hsize_t) :: dims(1)

        dims(1) = 1
        group_id = file%current_group   
        
        call h5aopen_f(group_id, dname, attr_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_real'
            write(*,*) 'h5aopen_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'group_id:', group_id
            write(*,*) 'attr_id:', attr_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
         
        call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, data, dims, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_real'
            write(*,*) 'h5aread_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'attr_id:', attr_id
            write(*,*) 'dims:', dims
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if  
            
        call h5aclose_f(attr_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_real'
            write(*,*) 'h5aclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'attr_id:', attr_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if  

      end subroutine read_real

      
      ! Subroutine to write real 1d arrays
      subroutine write_1dreal(file,data,dname)
        implicit none
        integer, parameter :: rank=1
        real*8, intent(in) :: data(:)
        character(*), intent(in) :: dname
        type(datafile), intent(in) :: file
        integer(HSIZE_T) :: dims(rank)
        integer(HID_T) :: filespace,d_id,group_id

        dims = shape(data)
        group_id = file%current_group

        call h5screate_simple_f(rank, dims, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_1dreal'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if  
        
        call h5dcreate_f(group_id, dname, H5T_NATIVE_DOUBLE, filespace, d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_1dreal'
            write(*,*) 'h5dcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if  
        
        call h5dwrite_f(d_id, H5T_NATIVE_DOUBLE, data, dims, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_1dreal'
            write(*,*) 'h5dwrite_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'dims:', dims
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if  
        
        call h5dclose_f(d_id,ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_1dreal'
            write(*,*) 'h5dclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if  
        
        call h5sclose_f(filespace,ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_1dreal'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if  
        
      end subroutine write_1dreal
      
      ! Subroutine to read real 1d arrays
      subroutine read_1dreal(file,data,dname)
        implicit none
        integer, parameter :: rank=1
        real*8, intent(out) :: data(:)
        character(*), intent(in) :: dname
        type(datafile), intent(in) :: file
        integer(HSIZE_T) :: ndims(rank)
        integer(HID_T) :: filespace,d_id,group_id

        ndims = shape(data)
        group_id = file%current_group

        call h5dopen_f(group_id, dname, d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_1dreal'
            write(*,*) 'h5dopen_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'group_id:', group_id
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if      
        
        call h5dread_f(d_id, H5T_NATIVE_DOUBLE, data, ndims, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_1dreal'
            write(*,*) 'h5dread_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'ndims:', ndims
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if    
        
        call h5dclose_f(d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_1dreal'
            write(*,*) 'h5dclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if    
        
      end subroutine read_1dreal
      
      ! Subroutine to write integer 2d arrays 
      subroutine write_2dint(file,data,dimsf,dname, n1s, n1e, n2s, n2e)
        implicit none
        integer, parameter :: ndims=2
        integer, intent(in) :: data(:,:)
        character(*), intent(in) :: dname
        type(datafile), intent(in) :: file
        integer(HSIZE_T), dimension(ndims), intent(in) :: dimsf
        integer, intent(in) :: n1s, n1e, n2s, n2e
        integer(kind=hid_t) :: x_id, d_id, group_id
        integer(kind=hid_t) :: memspace, filespace
        ! Local hyper slab info
        integer(kind=hsize_t) :: count(ndims), offset(ndims)

        group_id = file%current_group

        ! Create the data space for the dataset.
        call h5screate_simple_f(ndims, dimsf, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dint'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
        call h5dcreate_f(group_id, dname, H5T_NATIVE_INTEGER, filespace, d_id, &
                         ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dint'
            write(*,*) 'h5screate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'group_id:', group_id
            write(*,*) 'filespace:', filespace
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        count(1) = n1e-n1s+1 
        count(2) = n2e-n2s+1 
        offset(1) = n1s-1
        offset(2) = n2s-1

        call h5screate_simple_f(ndims, count, memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dint'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'count:', count
            write(*,*) 'memspace:', memspace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 

        ! Select hyperslab in the file.
        call h5dget_space_f(d_id, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dint'
            write(*,*) 'h5dget_space_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     
        
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dint'
            write(*,*) 'h5sselect_hyperslab_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
            write(*,*) 'offset:', offset
            write(*,*) 'count:', count
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if

        ! Create a data transfer property
        call h5pcreate_f(H5P_DATASET_XFER_F, x_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dint'
            write(*,*) 'h5pcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
        
        call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dint'
            write(*,*) 'h5pset_dxpl_mpio_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        ! Write the data
        call h5dwrite_f(d_id, H5T_NATIVE_INTEGER, data, dimsf, ierr,    &
                        file_space_id=filespace, mem_space_id=memspace, &
                        xfer_prp=x_id)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dint'
            write(*,*) 'h5dwrite_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id', d_id
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace
            write(*,*) 'memspace:', memspace
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if

        ! Close everything and exit
        call h5dclose_f(d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dint'
            write(*,*) 'h5dclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5sclose_f(filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dint'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
        
        call h5sclose_f(memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dint'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'memspace:', memspace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
        
        call h5pclose_f(x_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dint'
            write(*,*) 'h5pclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 

      end subroutine write_2dint
      
      ! Subroutine to read integer 2d arrays 
      subroutine read_2dint(file,data,dimsf,dname, n1s, n1e, n2s, n2e)
        implicit none
        integer, parameter :: ndims=2
        integer, intent(out) :: data(:,:)
        character(*), intent(in) :: dname
        type(datafile), intent(in) :: file
        integer(HSIZE_T), dimension(ndims), intent(in) :: dimsf
        integer, intent(in) :: n1s, n1e, n2s, n2e
        integer(kind=hid_t) :: x_id, d_id, group_id
        integer(kind=hid_t) :: memspace, filespace
        ! Local hyper slab info
        integer(kind=hsize_t) :: count(ndims), offset(ndims)

        group_id = file%current_group

        ! Create the data space for the  dataset.
        call h5screate_simple_f(ndims, dimsf, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dint'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
        call h5dopen_f(group_id, dname, d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dint'
            write(*,*) 'h5dopen_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'group_id:', group_id
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 

        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        count(1) = n1e-n1s+1 
        count(2) = n2e-n2s+1 
        offset(1) = n1s-1
        offset(2) = n2s-1

        call h5screate_simple_f(ndims, count, memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dint'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'count:', count
            write(*,*) 'memspace:', memspace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 

        ! Select hyperslab in the file.
        call h5dget_space_f(d_id, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dint'
            write(*,*) 'h5dget_space_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dint'
            write(*,*) 'h5sselect_hyperslab_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
            write(*,*) 'offset:', offset
            write(*,*) 'count:', count
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
        ! Create a data transfer property
        call h5pcreate_f(H5P_DATASET_XFER_F, x_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dint'
            write(*,*) 'h5pcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
            
        call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dint'
            write(*,*) 'h5pset_dxpl_mpio_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
        ! Write the data
        call h5dread_f(d_id, H5T_NATIVE_INTEGER, data, dimsf, ierr,    &
                       file_space_id=filespace, mem_space_id=memspace, &
                       xfer_prp=x_id)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dint'
            write(*,*) 'h5dread_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id', d_id
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace
            write(*,*) 'memspace:', memspace
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
                      

        ! Close everything and exit
        call h5dclose_f(d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dint'
            write(*,*) 'h5dclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if      
        
        call h5sclose_f(filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dint'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if      
        
        call h5sclose_f(memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dint'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'memspace:', memspace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if      
        
        call h5pclose_f(x_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dint'
            write(*,*) 'h5pclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if      

      end subroutine read_2dint      

      ! Subroutine to write real 2d arrays 
      subroutine write_2dreal(file,data,dimsf,dname, n1s, n1e, n2s, n2e)
        implicit none
        integer, parameter :: ndims=2
        real*8, intent(in) :: data(:,:)
        character(*), intent(in) :: dname
        type(datafile), intent(in) :: file
        integer(HSIZE_T), dimension(ndims), intent(in) :: dimsf
        integer, intent(in) :: n1s, n1e, n2s, n2e
        integer(kind=hid_t) :: x_id, d_id, group_id
        integer(kind=hid_t) :: memspace, filespace
        ! Local hyper slab info
        integer(kind=hsize_t) :: count(ndims), offset(ndims)

        group_id = file%current_group

        ! Create the data space for the  dataset.
        call h5screate_simple_f(ndims, dimsf, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dreal'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5dcreate_f(group_id, dname, H5T_IEEE_F64LE, filespace, d_id, &
                         ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dreal'
            write(*,*) 'h5dcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if

        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        count(1) = n1e-n1s+1 
        count(2) = n2e-n2s+1 
        offset(1) = n1s-1
        offset(2) = n2s-1

        call h5screate_simple_f(ndims, count, memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dreal'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'count:', count
            write(*,*) 'memspace:', memspace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if

        ! Select hyperslab in the file.
        call h5dget_space_f(d_id, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dreal'
            write(*,*) 'h5dget_space_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dreal'
            write(*,*) 'h5sselect_hyperslab_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
            write(*,*) 'offset:', offset
            write(*,*) 'count:', count
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
        
        ! Create a data transfer property
        call h5pcreate_f(H5P_DATASET_XFER_F, x_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dreal'
            write(*,*) 'h5pcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
            
        call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dreal'
            write(*,*) 'h5pset_dxpl_mpio_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
            
        ! Write the data
        call h5dwrite_f(d_id, H5T_IEEE_F64LE, data, dimsf, ierr,        &
                        file_space_id=filespace, mem_space_id=memspace, &
                        xfer_prp=x_id)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dreal'
            write(*,*) 'h5dwrite_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace
            write(*,*) 'memspace:', memspace
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if                    

        ! Close everything and exit
        call h5dclose_f(d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dreal'
            write(*,*) 'h5dclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
            
        call h5sclose_f(filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dreal'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
            
        call h5sclose_f(memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dreal'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'memspace:', memspace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
            
        call h5pclose_f(x_id, ierr) 
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_2dreal'
            write(*,*) 'h5pclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 

      end subroutine write_2dreal


      ! Subroutine to read real 2d arrays 
      subroutine read_2dreal(file,data,dimsf,dname, n1s, n1e, n2s, n2e)
        implicit none
        integer, parameter :: ndims=2
        real*8, intent(out) :: data(:,:)
        character(*), intent(in) :: dname
        type(datafile), intent(in) :: file
        integer(HSIZE_T), dimension(ndims), intent(in) :: dimsf
        integer, intent(in) :: n1s, n1e, n2s, n2e
        integer(kind=hid_t) :: x_id, d_id, group_id
        integer(kind=hid_t) :: memspace, filespace
        ! Local hyper slab info
        integer(kind=hsize_t) :: count(ndims), offset(ndims)

        group_id = file%current_group

        ! Create the data space for the  dataset.
        call h5screate_simple_f(ndims, dimsf, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dreal'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
        call h5dopen_f(group_id, dname, d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dreal'
            write(*,*) 'h5dopen_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'group_id:', group_id
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 

        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        count(1) = n1e-n1s+1 
        count(2) = n2e-n2s+1 
        offset(1) = n1s-1
        offset(2) = n2s-1

        call h5screate_simple_f(ndims, count, memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dreal'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'count:', count
            write(*,*) 'memspace:', memspace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 

        ! Select hyperslab in the file.
        call h5dget_space_f(d_id, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dreal'
            write(*,*) 'h5dget_space_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dreal'
            write(*,*) 'h5sselect_hyperslab_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
            write(*,*) 'offset:', offset
            write(*,*) 'count:', count
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
        ! Create a data transfer property
        call h5pcreate_f(H5P_DATASET_XFER_F, x_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dreal'
            write(*,*) 'h5pcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
            
        call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dreal'
            write(*,*) 'h5pset_dxpl_mpio_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
        ! Write the data
        call h5dread_f(d_id, H5T_IEEE_F64LE, data, dimsf, ierr,    &
                       file_space_id=filespace, mem_space_id=memspace, &
                       xfer_prp=x_id)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dreal'
            write(*,*) 'h5dread_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id', d_id
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace
            write(*,*) 'memspace:', memspace
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
                      

        ! Close everything and exit
        call h5dclose_f(d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dreal'
            write(*,*) 'h5dclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if      
        
        call h5sclose_f(filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dreal'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if      
        
        call h5sclose_f(memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dreal'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'memspace:', memspace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if      
        
        call h5pclose_f(x_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_2dreal'
            write(*,*) 'h5pclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if      

      end subroutine read_2dreal
      
      ! Subroutine to write integer 3d arrays 
      subroutine write_3dint(file,data,dimsf,dname,n1s,n1e,n2s,n2e,n3s,n3e)
        implicit none
        integer, parameter :: ndims=3
        integer, intent(in) :: data(:,:,:)
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
        call h5screate_simple_f(ndims, dimsf, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dint'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5dcreate_f(group_id, dname, H5T_NATIVE_INTEGER, filespace, d_id, &
                         ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dint'
            write(*,*) 'h5dcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'group_id:', group_id
            write(*,*) 'filespace:', filespace
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if

        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        count(1) = n1e-n1s+1
        count(2) = n2e-n2s+1
        count(3) = n3e-n3s+1
        offset(1) = n1s-1
        offset(2) = n2s-1
        offset(3) = n3s-1

        call h5screate_simple_f(ndims, count, memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dint'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'count:', count
            write(*,*) 'memspace:', memspace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
        

        ! Select hyperslab in the file.
        call h5dget_space_f(d_id, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dint'
            write(*,*) 'h5dget_space_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dint'
            write(*,*) 'h5sselect_hyperslab_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
            write(*,*) 'offset:', offset
            write(*,*) 'count:', count
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if

        ! Create a data transfer property
        call h5pcreate_f(H5P_DATASET_XFER_F, x_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dint'
            write(*,*) 'h5pcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dint'
            write(*,*) 'h5pset_dxpl_mpio_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        ! Write the data
        call h5dwrite_f(d_id, H5T_NATIVE_INTEGER, data, dimsf, ierr,    &
                        file_space_id=filespace, mem_space_id=memspace, &
                        xfer_prp=x_id)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dint'
            write(*,*) 'h5dwrite_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace  
            write(*,*) 'memspace:', memspace 
            write(*,*) 'x_id:', x_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if                    

        ! Close everything and exit
        call h5dclose_f(d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dint'
            write(*,*) 'h5dclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5sclose_f(filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dint'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if    
        
        call h5sclose_f(memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dint'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'memspace:', memspace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
            
        call h5pclose_f(x_id, ierr) 
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dint'
            write(*,*) 'h5pclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
      end subroutine write_3dint
      
      ! Subroutine to read integer 3d arrays 
      subroutine read_3dint(file,data,dimsf,dname,n1s,n1e,n2s,n2e,n3s,n3e)
        implicit none
        integer, parameter :: ndims=3
        integer, intent(out) :: data(:,:,:)
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
        call h5screate_simple_f(ndims, dimsf, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dint'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace  
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if       
            
        call h5dopen_f(group_id, dname, d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dint'
            write(*,*) 'h5dopen_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'group_id:', group_id
            write(*,*) 'd_id:', d_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if       
        
        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        count(1) = n1e-n1s+1
        count(2) = n2e-n2s+1
        count(3) = n3e-n3s+1
        offset(1) = n1s-1
        offset(2) = n2s-1
        offset(3) = n3s-1

        call h5screate_simple_f(ndims, count, memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dint'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'count:', count
            write(*,*) 'memspace:', memspace  
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
        ! Select hyperslab in the file.
        call h5dget_space_f(d_id, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dint'
            write(*,*) 'h5dget_space_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'filespace:', filespace  
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     
        
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dint'
            write(*,*) 'h5sselect_hyperslab_f ierr:',ierr
            write(*,*) 'filespace:', filespace 
            write(*,*) 'offset:', offset  
            write(*,*) 'count:', count
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if  
        
        ! Create a data transfer property
        call h5pcreate_f(H5P_DATASET_XFER_F, x_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dint'
            write(*,*) 'h5pcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if    
        
        call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dint'
            write(*,*) 'h5pset_dxpl_mpio_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if        
        
        ! Write the data
        call h5dread_f(d_id, H5T_NATIVE_INTEGER, data, dimsf, ierr,    &
                       file_space_id=filespace, mem_space_id=memspace, &
                       xfer_prp=x_id)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dint'
            write(*,*) 'h5dread_f ierr:',ierr
            write(*,*) 'd_id', d_id
            write(*,*) 'dimsf', dimsf
            write(*,*) 'filesspace:', filespace
            write(*,*) 'memspace:', memspace 
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if                     

        ! Close everything and exit
        call h5dclose_f(d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dint'
            write(*,*) 'h5dclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if       
        
        call h5sclose_f(filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dint'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     
        
        call h5sclose_f(memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dint'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'memspace:', memspace 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     
        
        call h5pclose_f(x_id, ierr) 
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dint'
            write(*,*) 'h5pclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     

      end subroutine read_3dint

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
        call h5screate_simple_f(ndims, dimsf, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dreal'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5dcreate_f(group_id, dname, H5T_IEEE_F64LE, filespace, d_id, &
                         ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dreal'
            write(*,*) 'h5dcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'group_id:', group_id
            write(*,*) 'filespace:', filespace
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if

        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        count(1) = n1e-n1s+1
        count(2) = n2e-n2s+1
        count(3) = n3e-n3s+1
        offset(1) = n1s-1
        offset(2) = n2s-1
        offset(3) = n3s-1

        call h5screate_simple_f(ndims, count, memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dreal'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'count:', count
            write(*,*) 'memspace:', memspace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
        

        ! Select hyperslab in the file.
        call h5dget_space_f(d_id, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dreal'
            write(*,*) 'h5dget_space_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dreal'
            write(*,*) 'h5sselect_hyperslab_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
            write(*,*) 'offset:', offset
            write(*,*) 'count:', count
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if

        ! Create a data transfer property
        call h5pcreate_f(H5P_DATASET_XFER_F, x_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dreal'
            write(*,*) 'h5pcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dreal'
            write(*,*) 'h5pset_dxpl_mpio_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        ! Write the data
        call h5dwrite_f(d_id, H5T_IEEE_F64LE, data, dimsf, ierr,    &
                        file_space_id=filespace, mem_space_id=memspace, &
                        xfer_prp=x_id)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dreal'
            write(*,*) 'h5dwrite_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace  
            write(*,*) 'memspace:', memspace 
            write(*,*) 'x_id:', x_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if                    

        ! Close everything and exit
        call h5dclose_f(d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dreal'
            write(*,*) 'h5dclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5sclose_f(filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dreal'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if    
        
        call h5sclose_f(memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dreal'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'memspace:', memspace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
            
        call h5pclose_f(x_id, ierr) 
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_3dreal'
            write(*,*) 'h5pclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
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
        call h5screate_simple_f(ndims, dimsf, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dreal'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace  
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if       
            
        call h5dopen_f(group_id, dname, d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dreal'
            write(*,*) 'h5dopen_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'group_id:', group_id
            write(*,*) 'd_id:', d_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if       
        
        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        count(1) = n1e-n1s+1
        count(2) = n2e-n2s+1
        count(3) = n3e-n3s+1
        offset(1) = n1s-1
        offset(2) = n2s-1
        offset(3) = n3s-1

        call h5screate_simple_f(ndims, count, memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dreal'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'count:', count
            write(*,*) 'memspace:', memspace  
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
        ! Select hyperslab in the file.
        call h5dget_space_f(d_id, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dreal'
            write(*,*) 'h5dget_space_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'filespace:', filespace  
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     
        
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dreal'
            write(*,*) 'h5sselect_hyperslab_f ierr:',ierr
            write(*,*) 'filespace:', filespace 
            write(*,*) 'offset:', offset  
            write(*,*) 'count:', count
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if  
        
        ! Create a data transfer property
        call h5pcreate_f(H5P_DATASET_XFER_F, x_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dreal'
            write(*,*) 'h5pcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if    
        
        call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dreal'
            write(*,*) 'h5pset_dxpl_mpio_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if        
        
        ! Write the data
        call h5dread_f(d_id, H5T_IEEE_F64LE, data, dimsf, ierr,    &
                       file_space_id=filespace, mem_space_id=memspace, &
                       xfer_prp=x_id)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dreal'
            write(*,*) 'h5dread_f ierr:',ierr
            write(*,*) 'd_id', d_id
            write(*,*) 'dimsf', dimsf
            write(*,*) 'filesspace:', filespace
            write(*,*) 'memspace:', memspace 
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if                     

        ! Close everything and exit
        call h5dclose_f(d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dreal'
            write(*,*) 'h5dclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if       
        
        call h5sclose_f(filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dreal'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     
        
        call h5sclose_f(memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dreal'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'memspace:', memspace 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     
        
        call h5pclose_f(x_id, ierr) 
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_3dreal'
            write(*,*) 'h5pclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     

      end subroutine read_3dreal


      ! Subroutine to write integer 4d arrays 
      subroutine write_4dint(file,data,dimsf,dname,n1s,n1e,n2s,n2e,n3s,n3e,n4s,n4e)
        implicit none
        integer, parameter :: ndims=4
        integer, intent(in) :: data(:,:,:,:)
        character(*), intent(in) :: dname
        type(datafile), intent(in) :: file
        integer(HSIZE_T), dimension(ndims), intent(in) :: dimsf
        integer, intent(in) :: n1s, n1e, n2s, n2e, n3s, n3e, n4s, n4e
        integer(kind=hid_t) :: x_id, d_id, group_id
        integer(kind=hid_t) :: memspace, filespace
        ! Local hyper slab info
        integer(kind=hsize_t) :: count(ndims), offset(ndims)

        group_id = file%current_group

        ! Create the data space for the  dataset.
        call h5screate_simple_f(ndims, dimsf, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dint'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5dcreate_f(group_id, dname, H5T_NATIVE_INTEGER, filespace, d_id, &
                         ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dint'
            write(*,*) 'h5dcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'group_id:', group_id
            write(*,*) 'filespace:', filespace
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if

        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        count(1) = n1e-n1s+1
        count(2) = n2e-n2s+1
        count(3) = n3e-n3s+1
        count(4) = n4e-n4s+1
        offset(1) = n1s-1
        offset(2) = n2s-1
        offset(3) = n3s-1
        offset(4) = n4s-1

        call h5screate_simple_f(ndims, count, memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dint'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'count:', count
            write(*,*) 'memspace:', memspace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
        

        ! Select hyperslab in the file.
        call h5dget_space_f(d_id, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dint'
            write(*,*) 'h5dget_space_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dint'
            write(*,*) 'h5sselect_hyperslab_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
            write(*,*) 'offset:', offset
            write(*,*) 'count:', count
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if

        ! Create a data transfer property
        call h5pcreate_f(H5P_DATASET_XFER_F, x_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dint'
            write(*,*) 'h5pcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dint'
            write(*,*) 'h5pset_dxpl_mpio_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        ! Write the data
        call h5dwrite_f(d_id, H5T_NATIVE_INTEGER, data, dimsf, ierr,    &
                        file_space_id=filespace, mem_space_id=memspace, &
                        xfer_prp=x_id)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dint'
            write(*,*) 'h5dwrite_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace  
            write(*,*) 'memspace:', memspace 
            write(*,*) 'x_id:', x_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if                    

        ! Close everything and exit
        call h5dclose_f(d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dint'
            write(*,*) 'h5dclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5sclose_f(filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dint'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if    
        
        call h5sclose_f(memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dint'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'memspace:', memspace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
            
        call h5pclose_f(x_id, ierr) 
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dint'
            write(*,*) 'h5pclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
      end subroutine write_4dint
      
      ! Subroutine to read integer 4d arrays 
      subroutine read_4dint(file,data,dimsf,dname,n1s,n1e,n2s,n2e,n3s,n3e,n4s,n4e)
        implicit none
        integer, parameter :: ndims=4
        integer, intent(out) :: data(:,:,:,:)
        character(*), intent(in) :: dname
        type(datafile), intent(in) :: file
        integer(HSIZE_T), dimension(ndims), intent(in) :: dimsf
        integer, intent(in) :: n1s, n1e, n2s, n2e, n3s, n3e, n4s, n4e
        integer(kind=hid_t) :: x_id, d_id, group_id
        integer(kind=hid_t) :: memspace, filespace
        ! Local hyper slab info
        integer(kind=hsize_t) :: count(ndims), offset(ndims)

        group_id = file%current_group

        ! Create the data space for the  dataset.
        call h5screate_simple_f(ndims, dimsf, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dint'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace  
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if       
            
        call h5dopen_f(group_id, dname, d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dint'
            write(*,*) 'h5dopen_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'group_id:', group_id
            write(*,*) 'd_id:', d_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if       
        
        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        count(1) = n1e-n1s+1
        count(2) = n2e-n2s+1
        count(3) = n3e-n3s+1
        count(4) = n4e-n4s+1
        offset(1) = n1s-1
        offset(2) = n2s-1
        offset(3) = n3s-1
        offset(4) = n4s-1

        call h5screate_simple_f(ndims, count, memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dint'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'count:', count
            write(*,*) 'memspace:', memspace  
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
        ! Select hyperslab in the file.
        call h5dget_space_f(d_id, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dint'
            write(*,*) 'h5dget_space_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'filespace:', filespace  
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     
        
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dint'
            write(*,*) 'h5sselect_hyperslab_f ierr:',ierr
            write(*,*) 'filespace:', filespace 
            write(*,*) 'offset:', offset  
            write(*,*) 'count:', count
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if  
        
        ! Create a data transfer property
        call h5pcreate_f(H5P_DATASET_XFER_F, x_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dint'
            write(*,*) 'h5pcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if    
        
        call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dint'
            write(*,*) 'h5pset_dxpl_mpio_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if        
        
        ! Write the data
        call h5dread_f(d_id, H5T_NATIVE_INTEGER, data, dimsf, ierr,    &
                       file_space_id=filespace, mem_space_id=memspace, &
                       xfer_prp=x_id)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dint'
            write(*,*) 'h5dread_f ierr:',ierr
            write(*,*) 'd_id', d_id
            write(*,*) 'dimsf', dimsf
            write(*,*) 'filesspace:', filespace
            write(*,*) 'memspace:', memspace 
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if                     

        ! Close everything and exit
        call h5dclose_f(d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dint'
            write(*,*) 'h5dclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if       
        
        call h5sclose_f(filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dint'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     
        
        call h5sclose_f(memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dint'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'memspace:', memspace 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     
        
        call h5pclose_f(x_id, ierr) 
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dint'
            write(*,*) 'h5pclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     

      end subroutine read_4dint

       ! Subroutine to write real 4d arrays 
      subroutine write_4dreal(file,data,dimsf,dname,n1s,n1e,n2s,n2e,n3s,n3e,n4s,n4e)
        implicit none
        integer, parameter :: ndims=4
        real*8, intent(in) :: data(:,:,:,:)
        character(*), intent(in) :: dname
        type(datafile), intent(in) :: file
        integer(HSIZE_T), dimension(ndims), intent(in) :: dimsf
        integer, intent(in) :: n1s, n1e, n2s, n2e, n3s, n3e, n4s, n4e
        integer(kind=hid_t) :: x_id, d_id, group_id
        integer(kind=hid_t) :: memspace, filespace
        ! Local hyper slab info
        integer(kind=hsize_t) :: count(ndims), offset(ndims)

        group_id = file%current_group

        ! Create the data space for the  dataset.
        call h5screate_simple_f(ndims, dimsf, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dreal'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5dcreate_f(group_id, dname, H5T_IEEE_F64LE, filespace, d_id, &
                         ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dreal'
            write(*,*) 'h5dcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'group_id:', group_id
            write(*,*) 'filespace:', filespace
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if

        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        count(1) = n1e-n1s+1
        count(2) = n2e-n2s+1
        count(3) = n3e-n3s+1
        count(4) = n4e-n4s+1
        offset(1) = n1s-1
        offset(2) = n2s-1
        offset(3) = n3s-1
        offset(4) = n4s-1

        call h5screate_simple_f(ndims, count, memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dreal'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'count:', count
            write(*,*) 'memspace:', memspace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
        

        ! Select hyperslab in the file.
        call h5dget_space_f(d_id, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dreal'
            write(*,*) 'h5dget_space_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dreal'
            write(*,*) 'h5sselect_hyperslab_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
            write(*,*) 'offset:', offset
            write(*,*) 'count:', count
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if

        ! Create a data transfer property
        call h5pcreate_f(H5P_DATASET_XFER_F, x_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dreal'
            write(*,*) 'h5pcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dreal'
            write(*,*) 'h5pset_dxpl_mpio_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        ! Write the data
        call h5dwrite_f(d_id, H5T_IEEE_F64LE, data, dimsf, ierr,    &
                        file_space_id=filespace, mem_space_id=memspace, &
                        xfer_prp=x_id)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dreal'
            write(*,*) 'h5dwrite_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace  
            write(*,*) 'memspace:', memspace 
            write(*,*) 'x_id:', x_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if                    

        ! Close everything and exit
        call h5dclose_f(d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dreal'
            write(*,*) 'h5dclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if
            
        call h5sclose_f(filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dreal'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if    
        
        call h5sclose_f(memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dreal'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'memspace:', memspace
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
            
        call h5pclose_f(x_id, ierr) 
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: write_4dreal'
            write(*,*) 'h5pclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
      end subroutine write_4dreal
      
      ! Subroutine to read real 4d arrays 
      subroutine read_4dreal(file,data,dimsf,dname,n1s,n1e,n2s,n2e,n3s,n3e,n4s,n4e)
        implicit none
        integer, parameter :: ndims=4
        real*8, intent(out) :: data(:,:,:,:)
        character(*), intent(in) :: dname
        type(datafile), intent(in) :: file
        integer(HSIZE_T), dimension(ndims), intent(in) :: dimsf
        integer, intent(in) :: n1s, n1e, n2s, n2e, n3s, n3e, n4s, n4e
        integer(kind=hid_t) :: x_id, d_id, group_id
        integer(kind=hid_t) :: memspace, filespace
        ! Local hyper slab info
        integer(kind=hsize_t) :: count(ndims), offset(ndims)

        group_id = file%current_group

        ! Create the data space for the  dataset.
        call h5screate_simple_f(ndims, dimsf, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dreal'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'dimsf:', dimsf
            write(*,*) 'filespace:', filespace  
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if       
            
        call h5dopen_f(group_id, dname, d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dreal'
            write(*,*) 'h5dopen_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'group_id:', group_id
            write(*,*) 'd_id:', d_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if       
        
        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        count(1) = n1e-n1s+1
        count(2) = n2e-n2s+1
        count(3) = n3e-n3s+1
        count(4) = n4e-n4s+1
        offset(1) = n1s-1
        offset(2) = n2s-1
        offset(3) = n3s-1
        offset(4) = n4s-1

        call h5screate_simple_f(ndims, count, memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dreal'
            write(*,*) 'h5screate_simple_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'ndims:', ndims
            write(*,*) 'count:', count
            write(*,*) 'memspace:', memspace  
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if 
        
        ! Select hyperslab in the file.
        call h5dget_space_f(d_id, filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dreal'
            write(*,*) 'h5dget_space_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id
            write(*,*) 'filespace:', filespace  
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     
        
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dreal'
            write(*,*) 'h5sselect_hyperslab_f ierr:',ierr
            write(*,*) 'filespace:', filespace 
            write(*,*) 'offset:', offset  
            write(*,*) 'count:', count
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if  
        
        ! Create a data transfer property
        call h5pcreate_f(H5P_DATASET_XFER_F, x_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dreal'
            write(*,*) 'h5pcreate_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if    
        
        call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dreal'
            write(*,*) 'h5pset_dxpl_mpio_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if        
        
        ! Write the data
        call h5dread_f(d_id, H5T_IEEE_F64LE, data, dimsf, ierr,    &
                       file_space_id=filespace, mem_space_id=memspace, &
                       xfer_prp=x_id)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dreal'
            write(*,*) 'h5dread_f ierr:',ierr
            write(*,*) 'd_id', d_id
            write(*,*) 'dimsf', dimsf
            write(*,*) 'filesspace:', filespace
            write(*,*) 'memspace:', memspace 
            write(*,*) 'x_id:', x_id
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if                     

        ! Close everything and exit
        call h5dclose_f(d_id, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dreal'
            write(*,*) 'h5dclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'd_id:', d_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if       
        
        call h5sclose_f(filespace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dreal'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'filespace:', filespace 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     
        
        call h5sclose_f(memspace, ierr)
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dreal'
            write(*,*) 'h5sclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'memspace:', memspace 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     
        
        call h5pclose_f(x_id, ierr) 
        if(ierr/=0) then
          if (myid==0) then
            write(*,*) '*******ERROR*******'
            write(*,*) 'SUBROUTINE: read_4dreal'
            write(*,*) 'h5pclose_f ierr:',ierr
            write(*,*) 'dname:  ', dname
            write(*,*) 'x_id:', x_id 
          end if
          call mpi_barrier(comm,mpi_er)
          call mpi_abort(comm,1,mpi_er)
        end if     

      end subroutine read_4dreal     

  end module mod_hdf5IO 