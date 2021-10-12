module file_copy
    implicit none

contains

    ! This fonction copy an ascii file. To use instead of the system calls
    ! Problem: some spaces are added in the copied file
    subroutine copy_ascii_file(file_dest, file_src,nline)
        implicit none
        integer, parameter  :: LINE_MAX_SIZE=500
        character(*)    :: file_src, file_dest
        character(LINE_MAX_SIZE)    :: currentLine
        integer :: isrc, idest, io
        integer, intent(out), optional ::  nline

        isrc=120
        idest=121

        open(unit=isrc, file=file_src)
        open(unit=idest, file=file_dest)

        read(isrc,'(a)', iostat=io)currentLine
        if (present(nline)) nline=0

        do while (io==0)
            if (present(nline)) nline=nline+1
            write(idest,'(a)') trim(currentLine)
            read(isrc,'(a)', iostat=io)currentLine
        end do

        close(isrc)
        close(idest)


    end subroutine copy_ascii_file

end module file_copy
