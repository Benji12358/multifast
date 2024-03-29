!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! VTK_mod.f90 --- VTK data file format
!! 
!! Auteur          : Jalel Chergui (LIMSI-CNRS) <Jalel.Chergui@limsi.fr>
!! Créé le         : Wed Jul 26 14:36:52 2006
!! Dern. mod. par  : Cyrille Bonamy (LEGI-CNRS) <cyrille.bonamy@legi.cnrs.fr>
!! Dern. mod. le   : Wed Apr 30 11:35:06 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module vtk
  implicit none

  private

  integer, parameter          :: s=selected_real_kind(6)
  integer, parameter          :: si=selected_int_kind(9)
  integer, parameter          :: d=selected_real_kind(12)! ATTENTION AVANT C
  character(len=1), parameter :: newline=achar(10)
  integer                     :: iproc=0, nb_procs=1
  integer           :: n_subdomain=0
  type, public :: VTK_file_handle
   private
   character(len=80) :: prefix
   integer           :: unit
   integer           :: ni, nj, nk
   integer           :: counter=0
   integer           :: restart=0
   logical           :: first=.true.
  end type VTK_file_handle

  private :: handle_error, handle_warning, handle_info

  private :: VTK_write_mesh_2d,   VTK_write_mesh_3d,   &
             VTK_write_mesh_particle_3d, &
             VTK_write_vector_2d, VTK_write_vector_3d, &
             VTK_write_mesh_uns_3d, VTK_write_scalar_uns_3d, &
             VTK_write_vector_uns_3d, &
             VTK_write_scalar_2d, VTK_write_scalar_3d

  public :: VTK_open_file,  &
            VTK_write_mesh, &
            VTK_write_var,  &
            VTK_close_file, &
            VTK_collect_file

  interface VTK_write_mesh
    module procedure VTK_write_mesh_2d,VTK_write_mesh_3d,VTK_write_mesh_uns_3d,&
                     VTK_write_mesh_particle_3d
  end interface

  interface VTK_write_var
    module procedure VTK_write_scalar_2d, VTK_write_scalar_3d, &
                     VTK_write_mesh_particle_3d, &
                     VTK_write_vector_2d, VTK_write_vector_3d, &
                     VTK_write_scalar_uns_3d, VTK_write_vector_uns_3d
  end interface
 
  contains

  subroutine handle_error(name, message)
    implicit none

    character(len=*), intent(in) :: name, message

    print '(/,"   *** Error *** ", A,": ", A,/)',name, message
    stop
  end subroutine handle_error

  subroutine handle_warning(name, message)
    implicit none
    character(len=*), intent(in) :: name, message

    print '(/,"   *** Warning *** ",A,": ", A,/)',name, message
  end subroutine handle_warning

  subroutine handle_info(name, message)
    implicit none
    character(len=*), intent(in) :: name, message

    print '(/,"   *** Info *** ",A,": ", A,/)',name, message
  end subroutine handle_info

    subroutine VTK_open_file(prefix, proc_rank, num_procs, restart, subdomain, fd)
    implicit none

    character(len=*), intent(in)         :: prefix
    integer, optional, intent(in)        :: proc_rank, num_procs, restart, subdomain
    type(VTK_file_handle), intent(inout) :: fd
    character(len=10) :: rank, snapshot, subdom
    character(len=80) :: f
    character(len=256):: MAIN_header
    integer           :: err
    logical           :: file_opened

    !... Looking for a none connected logical file unit.
    fd%prefix=trim(prefix)
    fd%unit=99
    inquire(unit=fd%unit, opened=file_opened)
    do while (file_opened .and. fd%unit /= 0)

      fd%unit = fd%unit - 1
      inquire(unit=fd%unit, opened=file_opened)

    end do
    if (fd%unit == 0 .and. file_opened) then

      call handle_error("VTK_open_file","All file units from 0 to 99 are already connected.")
      stop

    else

      if ( present(proc_rank) .and. present(num_procs) ) then

        iproc    = proc_rank
        nb_procs = num_procs

      else if ( present(proc_rank) ) then

        call handle_error("VTK_open_file","Both PROC_RANK and NUM_PROCS arguments must be present.")

      else if ( present(num_procs) ) then

        call handle_error("VTK_open_file","Both PROC_RANK and NUM_PROCS arguments must be present.")

      end if
      if ((fd%first) .and. (present(restart))) then
         fd%restart=restart
         fd%counter=restart
         fd%first=.false.
      end if
      fd%counter=fd%counter+1
      write(snapshot,'(i8)') fd%counter
      if ( present(proc_rank) ) then
        write(rank,'(i8)') iproc
        f=trim(fd%prefix)//"_"//trim(adjustl(rank))//"_"//trim(adjustl(snapshot))//".vtk"
      else if (present(subdomain)) then
        n_subdomain=max(subdomain,n_subdomain)
        write(subdom,'(i8)') subdomain
        f=trim(fd%prefix)//"_"//trim(adjustl(subdom))//"_"//trim(adjustl(snapshot))//".vtk"
      else
        f=trim(fd%prefix)//"_"//trim(adjustl(snapshot))//".vtk"
      end if
      open(unit=fd%unit, file=trim(adjustl(f)), form="UNFORMATTED", access="STREAM", status="replace", &
           action="write", convert='BIG_ENDIAN', iostat=err)
      if(err /= 0) print '("Problem creating file ",a,".")', trim(f)
    end if
    MAIN_header="# vtk DataFile Version 2.0"//newline//"(c) J.C. April 2009"//newline//"BINARY"//newline
    write(unit=fd%unit) trim(MAIN_header)
  end subroutine VTK_open_file


  subroutine VTK_write_mesh_2d(fd, x, y)
    implicit none

    type(VTK_file_handle), intent(inout)   :: fd
    real(kind=d), intent(in), dimension(:) :: x, y
    character(len=30)  :: buf1, buf2
    character(len=256) :: GRID_header
    integer, parameter :: nk=1

    fd%ni=size(x) ; fd%nj=size(y)
    write(buf1,'(i8," ",i8," ",i8)') fd%ni,fd%nj,nk
    GRID_header="DATASET RECTILINEAR_GRID"//newline//"DIMENSIONS "//trim(adjustl(buf1))//newline
    write(unit=fd%unit) trim(GRID_header)
    write(buf2,'(i8)') fd%ni
    GRID_header="X_COORDINATES "//trim(adjustl(buf2))//" float"//newline
    write(unit=fd%unit) trim(GRID_header),real(x(1:fd%ni),kind=s),newline
    write(buf2,'(i8)') fd%nj
    GRID_header="Y_COORDINATES "//trim(adjustl(buf2))//" float"//newline
    write(unit=fd%unit) trim(GRID_header),real(y(1:fd%nj),kind=s),newline
    !write(buf2,'(i8)') fd%nk
    write(buf2,'(i8)') nk
    GRID_header="Z_COORDINATES "//trim(adjustl(buf2))//" float"//newline
    write(unit=fd%unit) trim(GRID_header),real(0.0_d,kind=s),newline
    write(buf2,'(i8)') fd%ni*fd%nj
    GRID_header="POINT_DATA "//trim(adjustl(buf2))//newline
    write(unit=fd%unit) trim(GRID_header)
  end subroutine VTK_write_mesh_2d

  subroutine VTK_write_mesh_3d(fd, x, y, z)
    implicit none

    type(VTK_file_handle), intent(inout)   :: fd
    real(kind=d), intent(in), dimension(:) :: x, y, z
    character(len=30)  :: buf1, buf2
    character(len=256) :: GRID_header

    fd%ni=size(x) ; fd%nj=size(y) ; fd%nk=size(z)
    write(buf1,'(i8," ",i8," ",i8)') fd%ni,fd%nj,fd%nk
    GRID_header="DATASET RECTILINEAR_GRID"//newline//"DIMENSIONS "//trim(adjustl(buf1))//newline
    write(unit=fd%unit) trim(GRID_header)
    write(buf2,'(i8)') fd%ni
    GRID_header="X_COORDINATES "//trim(adjustl(buf2))//" float"//newline
    write(unit=fd%unit) trim(GRID_header),real(x(1:fd%ni),kind=s),newline
    write(buf2,'(i8)') fd%nj
    GRID_header="Y_COORDINATES "//trim(adjustl(buf2))//" float"//newline
    write(unit=fd%unit) trim(GRID_header),real(y(1:fd%nj),kind=s),newline
    write(buf2,'(i8)') fd%nk
    GRID_header="Z_COORDINATES "//trim(adjustl(buf2))//" float"//newline
    write(unit=fd%unit) trim(GRID_header),real(z(1:fd%nk),kind=s),newline
    write(buf2,'(i8)') fd%ni*fd%nj*fd%nk
    GRID_header="POINT_DATA "//trim(adjustl(buf2))//newline
    write(unit=fd%unit) trim(GRID_header)
  end subroutine VTK_write_mesh_3d

  subroutine VTK_write_mesh_particle_3d(fd, xyz, n)
    implicit none

    type(VTK_file_handle), intent(inout)   :: fd
    real(kind=d), intent(in), dimension(:,:) :: xyz
    integer, intent(in)                     :: n
    character(len=30)  :: buf1,buf2
    character(len=256) :: GRID_header
    integer,dimension(:),allocatable :: cell_type
    integer :: element_num,i

    fd%ni=n ;
    write(buf1,'(i8)') fd%ni !Nombre de POINTS/PARTICLES/BUBLES
    write(buf2,'(i8)') 2*fd%ni !Nombre de POINTS/PARTICLES/BUBLES
    GRID_header="DATASET UNSTRUCTURED_GRID"//newline//"POINTS "//trim(adjustl(buf1))//" float"//newline
    write(unit=fd%unit) trim(GRID_header),real(xyz,kind=s),newline
    GRID_HEADER="CELLS "//trim(adjustl(buf1))//" "//trim(adjustl(buf2))//newline
    allocate(cell_type(2*fd%ni))
    cell_type(1:2*fd%ni)=1;
    do i=1,fd%ni
      cell_type(2*i)=i-1;
    end do
    write(unit=fd%unit) trim(GRID_header),int(cell_type,kind=si),newline
    deallocate(cell_type)

    GRID_HEADER="CELL_TYPES "//trim(adjustl(buf1))//newline
    allocate(cell_type(fd%ni))
    cell_type(1:fd%ni)=1;
    write(unit=fd%unit) trim(GRID_header),int(cell_type,kind=si),newline
    deallocate(cell_type)
    GRID_header="POINT_DATA "//trim(adjustl(buf1))//newline
    write(unit=fd%unit) trim(GRID_header)
  end subroutine VTK_write_mesh_particle_3d

  subroutine VTK_write_mesh_uns_3d(fd, xyz,element_order,element_node,cell_type)
    implicit none

    type(VTK_file_handle), intent(inout)   :: fd
    integer,intent(in) :: element_order
    integer,intent(in),dimension(:,:) :: element_node
    integer,intent(in),dimension(:) :: cell_type
    real(kind=d), intent(in), dimension(:,:) :: xyz
    integer,dimension(:,:),allocatable :: towrite
    character(len=30)  :: buf1, buf2, buf3
    character(len=256) :: GRID_header
    integer :: element_num

    fd%ni=size(xyz,2) ;
    element_num=size(element_node,dim=2)
    write(buf1,'(i8)') fd%ni !Nombre de POINTS/NOEUDS
    write(buf2,'(i8)') element_num !Nombre de CELLULES
    write(buf3,'(i8)') element_num*(element_order+1) ! TAILLE CONNECTIVITE/REPRESENTATION CELLULES
    GRID_header="DATASET UNSTRUCTURED_GRID"//newline//"POINTS "//trim(adjustl(buf1))//" float"//newline
    write(unit=fd%unit) trim(GRID_header),real(xyz,kind=s),newline
    GRID_header="CELLS "//trim(adjustl(buf2))//" "//trim(adjustl(buf3))//newline
    allocate(towrite(element_order+1,element_num))
    towrite(1,1:element_num)=element_order
    towrite(2:element_order+1,1:element_num)=element_node(1:element_order,1:element_num)-1
    write(unit=fd%unit) trim(GRID_header),int(towrite,kind=si),newline
    deallocate(towrite)
    GRID_HEADER="CELL_TYPES "//trim(adjustl(buf2))//newline
    write(unit=fd%unit) trim(GRID_header),int(cell_type,kind=si),newline
    GRID_header="POINT_DATA "//trim(adjustl(buf1))//newline
    write(unit=fd%unit) trim(GRID_header)

  end subroutine VTK_write_mesh_uns_3d

  subroutine VTK_write_vector_2d(fd, name, vx, vy)
    implicit none

    type(VTK_file_handle), intent(in)        :: fd
    character(len=*), intent(in)             :: name
    real(kind=d), intent(in), dimension(:,:) :: vx, vy
    real(kind=d), allocatable, dimension(:,:,:) :: velocity
    integer                                     :: i, j, code=0
    character(len=256)                          :: uname, vname, VAR_header

    if ((size(vx,dim=1) /= fd%ni) .or. &
        (size(vx,dim=2) /= fd%nj)) call handle_warning("VTK_write_var","Incompatible X component and mesh sizes.")

    if ((size(vy,dim=1) /= fd%ni) .or. &
        (size(vy,dim=2) /= fd%nj)) call handle_warning("VTK_write_var","Incompatible Y component and mesh sizes.")

    if (.not.allocated(velocity)) then

      allocate(velocity(3,fd%ni,fd%nj),STAT=code)
      if ( code /= 0 ) &
      call handle_error("VTK_write_var","Not enough memory to allocate VELOCITY array")

    end if

    do j=1, fd%nj
      do i=1, fd%ni

          velocity(1,i,j) = vx(i,j)
          velocity(2,i,j) = vy(i,j)
          velocity(3,i,j) = 0.0_d

      end do
    end do

    VAR_header="VECTORS "//trim(adjustl(name))//" float "//newline
    write(unit=fd%unit) trim(VAR_header),real(velocity(:,1:fd%ni,1:fd%nj),kind=s),newline

    uname="X_"//name
    VAR_header="SCALARS "//trim(adjustl(uname))//" float "//newline//"LOOKUP_TABLE default"//newline
      write(unit=fd%unit) trim(VAR_header),real(velocity(1,1:fd%ni,1:fd%nj),kind=s),newline

    vname="Y_"//name
    VAR_header="SCALARS "//trim(adjustl(vname))//" float "//newline//"LOOKUP_TABLE default"//newline
      write(unit=fd%unit) trim(VAR_header),real(velocity(2,1:fd%ni,1:fd%nj),kind=s),newline

    if (allocated(velocity)) deallocate(velocity)
  end subroutine VTK_write_vector_2d

  subroutine VTK_write_vector_3d(fd, name, vx, vy, vz)
    implicit none

    type(VTK_file_handle), intent(in)          :: fd
    character(len=*), intent(in)               :: name
    real(kind=d), intent(in), dimension(:,:,:) :: vx, vy, vz
    real(kind=d), allocatable, dimension(:,:,:,:) :: velocity
    integer                                       :: i, j, k, code=0
    character(len=256)                            :: uname, vname, wname, VAR_header

    if ((size(vx,dim=1) /= fd%ni) .or. &
        (size(vx,dim=2) /= fd%nj) .or. &
        (size(vx,dim=3) /= fd%nk)) call handle_warning("VTK_write_var","Incompatible X component and mesh sizes.")

    if ((size(vy,dim=1) /= fd%ni) .or. &
        (size(vy,dim=2) /= fd%nj) .or. &
        (size(vy,dim=3) /= fd%nk)) call handle_warning("VTK_write_var","Incompatible Y component and mesh sizes.")

    if ((size(vz,dim=1) /= fd%ni) .or. &
        (size(vz,dim=2) /= fd%nj) .or. &
        (size(vz,dim=3) /= fd%nk)) call handle_warning("VTK_write_var","Incompatible Z component and mesh sizes.")

    if (.not.allocated(velocity)) then
      allocate(velocity(3,fd%ni,fd%nj,fd%nk),STAT=code)
      if ( code /= 0 ) &
      call handle_error("VTK_write_var","Not enough memory to allocate VELOCITY array")
    end if

    do k=1, fd%nk
      do j=1, fd%nj
        do i=1, fd%ni

          velocity(1,i,j,k) = vx(i,j,k)
          velocity(2,i,j,k) = vy(i,j,k)
          velocity(3,i,j,k) = vz(i,j,k)

        end do
      end do
    end do

    VAR_header="VECTORS "//trim(adjustl(name))//" float "//newline
    write(unit=fd%unit) trim(VAR_header),real(velocity(:,1:fd%ni,1:fd%nj,1:fd%nk),kind=s),newline

    uname="X_"//name
    VAR_header="SCALARS "//trim(adjustl(uname))//" float "//newline//"LOOKUP_TABLE default"//newline
      write(unit=fd%unit) trim(VAR_header),real(velocity(1,1:fd%ni,1:fd%nj,1:fd%nk),kind=s),newline

    vname="Y_"//name
    VAR_header="SCALARS "//trim(adjustl(vname))//" float "//newline//"LOOKUP_TABLE default"//newline
      write(unit=fd%unit) trim(VAR_header),real(velocity(2,1:fd%ni,1:fd%nj,1:fd%nk),kind=s),newline

    wname="Z_"//name
    VAR_header="SCALARS "//trim(adjustl(wname))//" float "//newline//"LOOKUP_TABLE default"//newline
      write(unit=fd%unit) trim(VAR_header),real(velocity(3,1:fd%ni,1:fd%nj,1:fd%nk),kind=s),newline

    if (allocated(velocity)) deallocate(velocity)
  end subroutine VTK_write_vector_3d

  subroutine VTK_write_vector_uns_3d(fd, name, vx, vy, vz, n)
    implicit none

    type(VTK_file_handle), intent(in)          :: fd
    character(len=*), intent(in)               :: name
    integer, intent(in)               :: n
    real(kind=d), intent(in), dimension(:) :: vx, vy, vz
    real(kind=d), allocatable, dimension(:,:) :: velocity
    integer                                       :: code=0
    character(len=256)                            :: uname, vname, wname, VAR_header

    write(*,*) 'VTK_write_vector_uns_3d'

    if ((n /= fd%ni)) call handle_warning("VTK_write_var","Incompatible component and mesh sizes.")

    if (.not.allocated(velocity)) then
      allocate(velocity(3,fd%ni),STAT=code)
      if ( code /= 0 ) &
      call handle_error("VTK_write_var","Not enough memory to allocate VELOCITY array")
    end if

    velocity(1,:) = vx
    velocity(2,:) = vy
    velocity(3,:) = vz

    VAR_header="VECTORS "//trim(adjustl(name))//" float "//newline
    write(unit=fd%unit) trim(VAR_header),real(velocity,kind=s),newline

    uname="X_"//name
    VAR_header="SCALARS "//trim(adjustl(uname))//" float "//newline//"LOOKUP_TABLE default"//newline
      write(unit=fd%unit) trim(VAR_header),real(velocity(1,:),kind=s),newline

    vname="Y_"//name
    VAR_header="SCALARS "//trim(adjustl(vname))//" float "//newline//"LOOKUP_TABLE default"//newline
      write(unit=fd%unit) trim(VAR_header),real(velocity(2,:),kind=s),newline

    wname="Z_"//name
    VAR_header="SCALARS "//trim(adjustl(wname))//" float "//newline//"LOOKUP_TABLE default"//newline
      write(unit=fd%unit) trim(VAR_header),real(velocity(3,:),kind=s),newline

    if (allocated(velocity)) deallocate(velocity)
  end subroutine VTK_write_vector_uns_3d

  subroutine VTK_write_scalar_2d(fd, name, field)
    implicit none

    type(VTK_file_handle), intent(in)        :: fd
    character(len=*), intent(in)             :: name
    real(kind=d), intent(in), dimension(:,:) :: field
    character(len=256) :: VAR_header

    if ((size(field,dim=1) /= fd%ni) .or. (size(field,dim=2) /= fd%nj)) &
       call handle_warning("VTK_write_var","Incompatible FIELD and MESH sizes.")

    VAR_header="SCALARS "//trim(adjustl(name))//" float "//newline//"LOOKUP_TABLE default"//newline
    write(unit=fd%unit) trim(VAR_header),real(field(1:fd%ni,1:fd%nj),kind=s),newline
  end subroutine VTK_write_scalar_2d

  subroutine VTK_write_scalar_3d(fd, name, field)
    implicit none

    type(VTK_file_handle), intent(in)          :: fd
    character(len=*), intent(in)               :: name
    real(kind=d), intent(in), dimension(:,:,:) :: field
    character(len=256) :: VAR_header

    if ((size(field,dim=1) /= fd%ni) .or. &
        (size(field,dim=2) /= fd%nj) .or. &
        (size(field,dim=3) /= fd%nk)) call handle_warning("VTK_write_var","Incompatible FIELD and MESH sizes.")

    VAR_header="SCALARS "//trim(adjustl(name))//" float "//newline//"LOOKUP_TABLE default"//newline
    write(unit=fd%unit) trim(VAR_header),real(field(1:fd%ni,1:fd%nj,1:fd%nk),kind=s),newline
  end subroutine VTK_write_scalar_3d

  subroutine VTK_write_scalar_uns_3d(fd, name, field)
    implicit none

    type(VTK_file_handle), intent(in)          :: fd
    character(len=*), intent(in)               :: name
    real(kind=d), intent(in), dimension(:) :: field
    character(len=256) :: VAR_header

    if ((size(field,dim=1) /= fd%ni)) call handle_warning("VTK_write_var","Incompatible FIELD and MESH sizes.")

    VAR_header="SCALARS "//trim(adjustl(name))//" float "//newline//"LOOKUP_TABLE default"//newline
    write(unit=fd%unit) trim(VAR_header),real(field(1:fd%ni),kind=s),newline
  end subroutine VTK_write_scalar_uns_3d

  subroutine VTK_close_file(fd)
   implicit none

   type(VTK_file_handle), intent(in) :: fd
   logical :: file_opened

   inquire(unit=fd%unit, opened=file_opened)
   if (file_opened) then
     close(unit=fd%unit)
   else
     call handle_warning("VTK_close_file","No such file to close. Please, check file descriptor.")
   end if
  end subroutine VTK_close_file

  subroutine VTK_collect_file(fd)
    implicit none

    type(VTK_file_handle), intent(inout) :: fd
    character(len=10) :: rank, snapshot, subdom
    character(len=80) :: f, vtrfile
    integer           :: shot, err, nt, np, k
    logical           :: file_opened

    !... Looking for a none connected logical file unit.
    if (iproc == 0) then
      fd%unit=99
      inquire(unit=fd%unit, opened=file_opened)
      do while (file_opened .and. fd%unit /= 0)
        fd%unit = fd%unit - 1
        inquire(unit=fd%unit, opened=file_opened)
      end do
      if (fd%unit == 0 .and. file_opened) then
        call handle_warning("VTK_open_file","warning, all file units from 0 to 99 are already connected.")
      else
        f=trim(adjustl(fd%prefix))//".pvd"
        open(unit=fd%unit, file=trim(adjustl(f)), form="FORMATTED", status="replace", &
           action="write", iostat=err)
        if(err /= 0) print '("VTK_collect_file: Error, problem creating file ",a,".")', trim(f)
        write(unit=fd%unit,fmt='(100A)')   '<?xml version="1.0"?>'
        write(unit=fd%unit,fmt='(100A)')   '<VTKFile type="Collection" version="0.1" &
          byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
        write(unit=fd%unit,fmt='(100A)')   '  <Collection>'
        nt=len_trim(fd%prefix)
        np=scan(STRING=fd%prefix, SET="/", BACK=.true.)
        vtrfile=fd%prefix(np+1:nt)
        if ( nb_procs == 1 ) then

          if ( n_subdomain == 0 ) then

          do shot = 1, fd%counter
            write(snapshot, '(i6)') shot
            write(unit=fd%unit,fmt='(100A)') '    <DataSet timestep="'//trim(adjustl(snapshot))//&
                                      &'" part="0'//'" file="'//trim(adjustl(vtrfile))//&
                                      &"_"//trim(adjustl(snapshot))//'.vtr"/>'
          end do

          else
          do k = 0, n_subdomain
            write(subdom, '(i6)') k
            do shot = 1, fd%counter
              write(snapshot, '(i6)') shot
              write(unit=fd%unit,fmt='(100A)') '    <DataSet timestep="'//trim(adjustl(snapshot))//&
                                        &'" part="'//trim(adjustl(subdom))//'" file="'//&
                                        &trim(adjustl(vtrfile))//"_"//trim(adjustl(subdom))//&
                                        &"_"//trim(adjustl(snapshot))//'.vtr"/>'
            end do
          end do

          endif

        else

          do k = 0, nb_procs-1
            write(rank, '(i6)') k
            do shot = 1, fd%counter
              write(snapshot, '(i6)') shot
              write(unit=fd%unit,fmt='(100A)') '    <DataSet timestep="'//trim(adjustl(snapshot))//&
                                        &'" part="'//trim(adjustl(rank))//'" file="'//&
                                        &trim(adjustl(vtrfile))//"_"//trim(adjustl(rank))//&
                                        &"_"//trim(adjustl(snapshot))//'.vtr"/>'
            end do
          end do

        end if
        write(unit=fd%unit,fmt='(100A)')    '  </Collection>'
        write(unit=fd%unit,fmt='(100A)')    '</VTKFile>'
        close(unit=fd%unit)
      end if
    end if
    fd%counter=0 ; fd%restart=0 ; fd%first=.true. ; iproc=0 ; nb_procs=1
  end subroutine VTK_collect_file
end module vtk
