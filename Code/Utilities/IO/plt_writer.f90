
module plt_file_writer

    implicit none

    CHARACTER(LEN=20),              PARAMETER :: &
    DP_FMT_DEFAULT = "(f10.4)",         &
    DP_FMT_HIGH_PRECISION="(f15.6)",    &
    INT_FMT = "(i5)",                   &
    DP_FMT_SCIENTIFIC="(ES21.14)"

    integer, parameter          :: PLT_INTEGER=0, PLT_DOUBLE_PRECISION=1
    integer, parameter          :: NB_ZONE_MAX=1000, NVAR_MAX=20, MAX_TITLE_SIZE=40, NB_FILES_MAX=10


    TYPE :: variable
        integer                     :: type
        character(20)               :: fmt
        character(MAX_TITLE_SIZE)   :: name
    END TYPE variable


    TYPE :: zone
        integer, dimension(3)           :: z_shape
        integer                         :: ndims
        integer                         :: nb_pts
        character(MAX_TITLE_SIZE)       :: title

    END TYPE zone


    TYPE :: set_of_int_values
        integer, dimension(:), pointer    :: values
        logical                             :: isNaN=.true.

    END TYPE set_of_int_values


    TYPE :: set_of_dp_values
        real*8, dimension(:), pointer    :: values
        logical                             :: isNaN=.true.

    END TYPE set_of_dp_values



    TYPE :: plt_file_view
        integer ::nb_zones=0, nb_vars=0, id=0

        type(set_of_int_values), dimension(NB_ZONE_MAX, NVAR_MAX)     :: int_values
        type(set_of_dp_values), dimension(NB_ZONE_MAX, NVAR_MAX) :: dp_values

        type(zone), dimension(NB_ZONE_MAX)     :: zones
        type(variable), dimension(NVAR_MAX)     :: variables

    END TYPE plt_file_view


    type(plt_file_view), dimension(NB_FILES_MAX), save, target         :: plt_files


    integer :: nb_of_files=0, next_id_index=0, UNIT_SHIFT=40
    integer, dimension(NB_FILES_MAX)    :: free_ids

    interface add_values
        module procedure add_int_values, add_dp_values
    end interface add_values


contains

    subroutine print_variable(var)
        implicit none
        type(variable) :: var

        write(*,*)

        select case (var%type)

            case (PLT_DOUBLE_PRECISION)

                write(*,*) 'Type: REAL*8'
                write(*,*)'Format :', var%fmt
                write(*,*)'Title :', var%name

            case (PLT_INTEGER)
                write(*,*) 'Type: INTEGER'
                write(*,*)'Title :', var%name

            case default

        end select

    end subroutine print_variable

    subroutine print_zone(z)
        implicit none
        type(zone) :: z
        integer     :: i

        write(*,*)
        write(*,*) 'Affichage dune zone'
        write(*,*) 'Nb pts par dim  : ', z%nb_pts
        write(*,*)'Title :', z%title

    end subroutine print_zone



    subroutine print_plt_file(plt_file)
        implicit none
        type(plt_file_view) :: plt_file

        write(*,*)
        write(*,*) 'Nb de zones  : ', plt_file%nb_zones
        write(*,*) 'Nb de variables  : ', plt_file%nb_vars
        write(*,*) 'ID  : ', plt_file%id

    end subroutine print_plt_file



    subroutine print_set_of_dp_values(SoDPv)
        implicit none
        type(set_of_dp_values) :: SoDPv

        write(*,*)
        write(*,*) 'Affichage dun ensemble de valeurs réelles'

        write(*,*) 'is NaN ?:', SoDPv%isNaN
        write(*,*) 'Valeurs:',SoDPv%values

    end subroutine print_set_of_dp_values



    subroutine print_set_of_int_values(SoIV)
        implicit none
        type(set_of_int_values) :: SoIV

        write(*,*)
        write(*,*) 'Affichage dun ensemble de valeurs entieres'

        write(*,*) 'is NaN ?:', SoIv%isNaN
        write(*,*) 'Valeurs:',SoIV%values

    end subroutine print_set_of_int_values



    subroutine get_string_for_value(plt_file, zone_key_arg, var_key_arg, value_index, str)
        implicit none
        character(*)    ::str
        integer             :: zone_key_arg, var_key_arg, value_index
        type(set_of_int_values) :: SoIV
        type(set_of_dp_values) :: SoDPV
        type(variable)          :: var
        type(plt_file_view)     :: plt_file

        var=plt_file%variables(var_key_arg)
        str="   NaN"

        if (var%type==PLT_INTEGER) then
            SoIV=plt_file%int_values(zone_key_arg, var_key_arg)

            if (.not. SoIV%isNaN) then
                write(str,var%fmt) SoIV%values(value_index)
            else
                str="    NaN"
            end if

        end if

        if (var%type==PLT_DOUBLE_PRECISION) then

            SoDPV=plt_file%dp_values(zone_key_arg, var_key_arg)

            if (.not. SoDPV%isNaN) then
                write(str,var%fmt) SoDPV%values(value_index)
            else
                str="    NaN"
            end if

        end if

    end subroutine get_string_for_value




    subroutine add_zone(plt_file_id, zone_id, zone_title, n1)
        implicit none
        character(*)                :: zone_title
        integer                     :: n1, zone_id, plt_file_id
        integer, dimension(3)       :: shape_1D
        type(variable) :: ds_null
        type(plt_file_view), pointer :: pf

        shape_1D=(/n1,0,0/)

        pf=>plt_files(plt_file_id)
        pf%nb_zones=pf%nb_zones+1
        zone_id=pf%nb_zones

        pf%zones(pf%nb_zones)=zone(z_shape=shape_1D, ndims=1, nb_pts=n1, title=zone_title)

    end subroutine add_zone


    subroutine add_2Dzone(plt_file_id, zone_id, zone_title, n1, n2)
        implicit none
        character(*)                :: zone_title
        integer                     :: n1, n2, zone_id, plt_file_id
        integer, dimension(3)       :: shape_2D
        type(variable) :: ds_null
        type(plt_file_view), pointer :: pf

        shape_2D=(/n1,n2,0/)

        pf=>plt_files(plt_file_id)
        pf%nb_zones=pf%nb_zones+1
        zone_id=pf%nb_zones

        pf%zones(pf%nb_zones)=zone(z_shape=shape_2D, ndims=2, nb_pts=n1*n2, title=zone_title)

    end subroutine add_2Dzone


    subroutine add_3Dzone(plt_file_id, zone_id, zone_title, n1, n2, n3)
        implicit none
        character(*)                :: zone_title
        integer                     :: n1, n2, n3, zone_id, plt_file_id
        integer, dimension(3)       :: shape_3D
        type(variable) :: ds_null
        type(plt_file_view), pointer :: pf

        shape_3D=(/n1,n2,n3/)

        pf=>plt_files(plt_file_id)
        pf%nb_zones=pf%nb_zones+1
        zone_id=pf%nb_zones

        pf%zones(pf%nb_zones)=zone(z_shape=shape_3D, ndims=3, nb_pts=n1*n2*n3, title=zone_title)

    end subroutine add_3Dzone




    subroutine add_variable(plt_file_id, var_id, var_title, values_type, values_format)
        implicit none
        character(*)                :: var_title, values_format
        integer                     :: var_id, values_type, plt_file_id
        type(variable) :: new_var
        type(plt_file_view), pointer :: pf

        pf=>plt_files(plt_file_id)


        pf%nb_vars=pf%nb_vars+1
        var_id=pf%nb_vars
        pf%variables(var_id)=variable(type=values_type, fmt=values_format, name=var_title)

    end subroutine add_variable




    subroutine add_dp_values(plt_file_id, zone_key_arg, var_key_arg, values_arg, nb_values)
        implicit none
        integer                     :: zone_key_arg, var_key_arg, nb_values, plt_file_id
        real*8, dimension(:)       :: values_arg

        real*8, dimension(:), pointer       :: values_ptr

        type(zone)                  :: z
        type(variable)                  :: var
        character(30)                   :: error
        type(plt_file_view), pointer :: pf

        pf=>plt_files(plt_file_id)

        z=pf%zones(zone_key_arg)
        var=pf%variables(var_key_arg)

        error=""

        if (nb_values.ne.z%nb_pts) then
            error="NOMBRE DE PTS INCORRECT"
        end if

        if (var%type.ne.PLT_DOUBLE_PRECISION) then
            error="TYPE INCORRECT"
        end if

        if ((nb_values.eq.z%nb_pts).and.(var%type.eq.PLT_DOUBLE_PRECISION)) then

            allocate(values_ptr(nb_values))
            values_ptr(1:nb_values)=values_arg(1:nb_values)

            pf%dp_values(zone_key_arg, var_key_arg)=    &
            set_of_dp_values(values=values_ptr, isNaN=.false.)

        else

            write(*,*) 'ERREUR:'//trim(error)

        end if

    end subroutine add_dp_values

    subroutine add_int_values(plt_file_id, zone_key_arg, var_key_arg, values_arg, nb_values)
        implicit none
        integer                     :: zone_key_arg, var_key_arg, nb_values, plt_file_id
        integer, dimension(:)       :: values_arg
        integer, dimension(:), pointer       :: values_ptr
        type(zone)                  :: z
        type(variable)                  :: var
        character(30)                   :: error
        type(plt_file_view), pointer :: pf

        pf=>plt_files(plt_file_id)

        z=pf%zones(zone_key_arg)
        var=pf%variables(var_key_arg)

        error=""

        if (nb_values.ne.z%nb_pts) then
            error="NOMBRE DE PTS INCORRECT"
        end if

        if (var%type.ne.PLT_INTEGER) then
            error="TYPE INCORRECT"
        end if

        if ((nb_values.eq.z%nb_pts).and.(var%type.eq.PLT_INTEGER)) then

            allocate(values_ptr(nb_values))
            values_ptr(1:nb_values)=values_arg(1:nb_values)

            pf%int_values(zone_key_arg, var_key_arg)=   &
            set_of_int_values(values=values_ptr, isNaN=.false.)

        else

            write(*,*) 'ERREUR:'//trim(error)

        end if

    end subroutine add_int_values



    subroutine open_plt_file(id, path, file_name)
        implicit none
        character(*) :: path, file_name
        character(200)  :: file_url
        integer     :: id

        call alloc_plt_file(id)
        plt_files(id)%id=id

        file_url=trim(path)//trim(file_name)//".plt"
        open (unit = id+UNIT_SHIFT, file = file_url)

    contains

        subroutine alloc_plt_file(id)
            implicit none
            integer id

            if (next_id_index.gt.0) then
                id=free_ids(next_id_index)
                free_ids(next_id_index)=0
                next_id_index=next_id_index-1
            else
                nb_of_files=nb_of_files+1
                id=nb_of_files
            end if

        end subroutine alloc_plt_file


    end subroutine open_plt_file



    subroutine close_plt_file(id)
        implicit none
        integer     :: id, z, j
        type(zone)  :: zn
        type(plt_file_view), pointer :: pf

        pf=>plt_files(id)


        do z = 1, pf%nb_zones
            zn=pf%zones(z)

            do j = 1, pf%nb_vars
                call free_values(pf, zone_key_arg=z, var_key_arg=j)
            end do
        end do


        pf%nb_zones=0
        pf%nb_vars=0
        pf%id=0

        next_id_index=next_id_index+1
        free_ids(next_id_index)=id

        close(id)

    contains





        subroutine free_values(plt_file, zone_key_arg, var_key_arg)
            implicit none
            integer             :: zone_key_arg, var_key_arg, value_index
            type(set_of_int_values) :: SoIV
            type(set_of_dp_values) :: SoDPV
            type(variable)          :: var
            type(plt_file_view)     :: plt_file

            var=plt_file%variables(var_key_arg)

            if (var%type==PLT_INTEGER) then
                SoIV=plt_file%int_values(zone_key_arg, var_key_arg)
                plt_file%int_values(zone_key_arg, var_key_arg)%isNaN=.true.
                if (associated(plt_file%int_values(zone_key_arg, var_key_arg)%values)) then
                    deallocate(plt_file%int_values(zone_key_arg, var_key_arg)%values)
                end if

            end if

            if (var%type==PLT_DOUBLE_PRECISION) then
                SoDPV=plt_file%dp_values(zone_key_arg, var_key_arg)
                plt_file%dp_values(zone_key_arg, var_key_arg)%isNaN=.true.
                if (associated(plt_file%dp_values(zone_key_arg, var_key_arg)%values)) then
                    deallocate(plt_file%dp_values(zone_key_arg, var_key_arg)%values)
                end if

            end if

        end subroutine free_values

    end subroutine close_plt_file



    subroutine write_in_file(plt_file_id, title)
        implicit none
        character(*) :: title
        character(40), dimension(NVAR_MAX)  :: line_values
        integer     :: file_id, i, j, z, plt_file_id
        type(zone)  :: zn
        type(plt_file_view), pointer :: pf

        pf=>plt_files(plt_file_id)


        call create_plt_heading(pf%id+UNIT_SHIFT, title, pf%nb_vars, pf%variables(1:pf%nb_vars)%name, pf%zones(1)%title)

        do z = 1, pf%nb_zones
            zn=pf%zones(z)

            call write_current_zone

            do i = 1, zn%nb_pts
                do j = 1, pf%nb_vars
                    call get_string_for_value(pf, zone_key_arg=z, var_key_arg=j, value_index=i, str=line_values(j))
                end do

                call add_line_to_plt(pf%id+UNIT_SHIFT, pf%nb_vars, line_values)

            end do

        end do

        close(pf%id+UNIT_SHIFT)

    contains

        subroutine create_plt_heading(file_id, title, n_columns,labels, zone_name)

            implicit none
            character(*) :: zone_name, title
            integer :: n_columns, i, j, file_id
            character(*), dimension(*)::labels
            character(1000) :: line

            write(file_id,*)'TITLE = "'//title//'"'
            write(file_id,*)'VARIABLES = "'//trim(labels(1))//'"'

            do i = 2, n_columns
                write(file_id,*)'"'//trim(labels(i))//'"'
            end do

            return

        end subroutine


        subroutine write_current_zone()
            implicit none

1           format('ZONE T="', a,'"')
2           format('ZONE T="', a,'", I=', i4, ', J=', i4, ', F=POINT')
3           format('ZONE T="', a,'", I=', i4, ', J=', i4, ', K=', i4, ', F=POINT')

            select case (zn%ndims)

                case (1)
                    write(pf%id+UNIT_SHIFT,1)trim(zn%title)

                case (2)
                    write(pf%id+UNIT_SHIFT,2)trim(zn%title), zn%z_shape(1), zn%z_shape(2)

                case (3)
                    write(pf%id+UNIT_SHIFT,3)trim(zn%title), zn%z_shape(1), zn%z_shape(2), zn%z_shape(3)

                case default

            end select


        end subroutine write_current_zone

        subroutine add_line_to_plt(file_id, n_columns, line_data)

            implicit none
            integer :: file_id, n_lines, n_columns, i, j, line_cursor, column_size
            character(*), dimension(*):: line_data
            character(1000) :: line

            line_cursor=0
            column_size=30

            do i = 1, n_columns
                line=line(1:line_cursor)//line_data(i)
                line_cursor=line_cursor+column_size
            end do

            write(file_id,*) line(1:line_cursor)

            return

        end subroutine


    end subroutine write_in_file



end module plt_file_writer
