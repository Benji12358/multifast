


module object_file_reader
    implicit none

contains

    subroutine read_object_size(file_path, nb_vertex, nb_faces)
        implicit none
        character(*)    :: file_path
        integer             :: nb_vertex, nb_faces

        character(200)  :: current_line
        integer             :: ivertex, iface

        open(15, file=trim(file_path))

        ivertex=0
        iface=0


        do while (.true.)

            read(15, '(a)', end=10) current_line

            !write(*,*)
            !write(*,*) current_line

            if(current_line(1:1)=="v") then
                ivertex=ivertex+1
            endif

            if(current_line(1:1)=="f") then
                iface=iface+1
            endif

        end do

10      close(15)

        nb_vertex=ivertex
        nb_faces=iface


    end subroutine read_object_size

    subroutine read_object(file_path, vertex, faces, nb_vertex, nb_faces, verbose)
        implicit none
        character(*)    :: file_path
        integer             :: nb_vertex, nb_faces
        character(200)  :: current_line
        character(200)      :: point_str, face_str
        integer, optional             :: verbose

        real*8, dimension(:,:), allocatable              :: vertex
        integer             :: ivertex, iface, i
        integer, dimension(:,:), allocatable             :: faces
        real*8  :: center(3)

        if (allocated(vertex)) then
            deallocate(vertex)
        endif

        if (allocated(faces)) then
            deallocate(faces)
        endif

        allocate(vertex(nb_vertex,3))
        allocate(faces(nb_faces,3))

        open(15, file=trim(file_path))

        ivertex=0
        iface=0


        do while (.true.)

            read(15, '(a)', end=10) current_line

            if (present(verbose)) write(*,*)
            if (present(verbose)) write(*,*) current_line

            if(current_line(1:1)=="v") then
                ivertex=ivertex+1
                point_str=current_line(3:LEN_TRIM(current_line))
                read(point_str,*) vertex(ivertex, 1), vertex(ivertex, 2), vertex(ivertex, 3)
                if (present(verbose)) write(*,*)'POINT:', vertex(ivertex, :)
            endif

            if(current_line(1:1)=="f") then
                iface=iface+1
                face_str=current_line(3:LEN_TRIM(current_line))
                read(face_str,*) faces(iface, 1), faces(iface, 2), faces(iface, 3)
                if (present(verbose)) then
                    write(*,*)'Face ', iface, " composee des sommets :"
                    write(*,*)"------", faces(iface, 1)
                    write(*,*)"------", faces(iface, 2)
                    write(*,*)"------", faces(iface, 3)
                endif
            endif

        end do

10      close(15)


        ! Recenter the object
        center=0.d0
        do i = 1, nb_vertex
            center=center+vertex(i, :)
        end do
        center=center/nb_vertex

        do i = 1, nb_vertex
            vertex(i, :)=vertex(i, :)-center
        end do

    end subroutine read_object


end module object_file_reader
