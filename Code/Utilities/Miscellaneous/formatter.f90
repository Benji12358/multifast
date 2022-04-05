module formatter
    implicit none

contains

    subroutine format_integer(number, str, str_size)
        implicit none
        character(*)    :: str
        integer         :: str_size, number

        character(str_size) :: tmp_str
        integer     :: i

        write(str, "(i10.10)")number

        do i = 1, str_size
            if (ichar(str(i:i))==ichar('0')) then
                str(i:i)=' '
            end if
        end do

        do i = 1, str_size
            write(*,*)"str", i, "=",str(i:i)
        end do
    end subroutine format_integer

end module formatter
