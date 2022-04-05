! This module provide buffered IO for 2d fields. The 2d fields are stored in memory (buffer)
! until the max number of fields is reached. The buffers are then stored in a file in order to be used thereafter.
! The module can be used in two ways: READ or WRITE:

!   - WRITE: Store 2d fields in buffer then, after the maximum storage limits is reached, save 2d fields in files.
!   - READ: Read 2d fields from buffer. The fields in buffer are previously loaded from files. When the buffers are full,
!   they are updated (overwritten) by reading next files.

! Use :
! call init(mode, n, recovery_dir_arg, working_dir_arg) :
!               Open a session in READ or WRIDE mode (mode input argument).
!               n is the size of the buffers.
!               recovery_dir_arg is the directory that will be used to restore buffer data if continuing from previous run.
!               working_arg is the directory were the buffers will be stored or read when buffers limit size is exceeded.

! call open_buffer(id, "fieldname"): open buffer identified by id. The buffer array are allocated and the output argument
! id will be used to access to the buffer from the client code.

!call add_2dfield(id, field): Store the input 2d field in the buffer identified by "id".
!call get_2dfield(id, field): Fill the input 2d field from the buffer identified by "id".




module flow_buffer_data
    implicit none

    type buffer_wrap
        real*8, dimension(:,:,:), pointer   :: data
        character(20)                       :: name
    end type buffer_wrap

    integer, parameter  :: MAX_BUFFER=20
    integer             :: nb_buffer=0

    integer             :: cursor=1
    integer             :: cursor_mod=1

    type(buffer_wrap), dimension(MAX_BUFFER)    :: buffer_list
    integer                                     :: nb_frames

    integer, parameter                          :: READ_MODE=0, WRITE_MODE=1
    integer                                     :: access_mode

    logical :: verbose=.true.

    character(200)  :: recovery_dir, working_dir

contains

end module flow_buffer_data


module flow_buffer_dao

    use flow_buffer_data
contains


    subroutine load_buffer(id, buffer_file)

        use decomp_2d

        use HDF5_IO

        implicit none
        integer, intent(in)         :: id
        character(200), intent(in)  :: buffer_file

        logical     :: buffer_file_exists
        character(20)                                                           :: cur_buffer_name

        cur_buffer_name=trim(buffer_list(id)%name)

        inquire( file=trim(buffer_file)//".h5", exist=buffer_file_exists)

        if(buffer_file_exists) then

            if (verbose) then
                open(15,file="inoutflow2", position="append")
                if(nrank==0) write(15,*)"Load buffer ", trim(cur_buffer_name), "from file:", trim(buffer_file)
                close(15)
            endif

            call hdf_read_3Dfield(buffer_file, buffer_list(id)%data(:,:,:), cur_buffer_name, nb_frames, ny_global, nz_global, 1,nb_frames, xstart(2),xend(2), xstart(3),xend(3))

        endif

    end subroutine load_buffer

    subroutine create_buffer_file(buffer_file)

        use decomp_2d

        use HDF5_IO

        implicit none
        character(200), intent(in)  :: buffer_file

        if(nrank==0)  call hdf_create_emptyfile(buffer_file)


    end subroutine create_buffer_file


    subroutine save_buffer(id, buffer_file)

        use decomp_2d

        use HDF5_IO

        implicit none
        integer, intent(in)         :: id
        character(200), intent(in)  :: buffer_file

        logical     :: buffer_file_exists
        real*8, dimension(nb_frames, xstart(2):xend(2), xstart(3):xend(3))    :: cur_buffer
        character(20)                                                           :: cur_buffer_name

        cur_buffer=buffer_list(id)%data(:,:,:)
        cur_buffer_name=trim(buffer_list(id)%name)

        inquire( file=trim(buffer_file)//".h5", exist=buffer_file_exists)


        if (verbose) then
            open(15,file="inoutflow2", position="append")
            if(nrank==0) write(15,*)"export buffer ", trim(cur_buffer_name), "in file", trim(buffer_file)
            close(15)
        endif


        call hdf_add_3Dfield(buffer_file, cur_buffer, cur_buffer_name, nb_frames, ny_global, nz_global, 1,nb_frames, xstart(2),xend(2), xstart(3),xend(3))

    end subroutine save_buffer


    subroutine flush_buffer(id, buffer_file)

        use decomp_2d

        use HDF5_IO

        implicit none
        integer, intent(in)         :: id
        character(200), intent(in)  :: buffer_file


        call save_buffer(id, buffer_file)
        buffer_list(id)%data(:,:,:)=0.d0

    end subroutine flush_buffer


    subroutine write_cursorfile()

        use decomp_2d

        use run_ctxt_data

        implicit none
        character(200)  :: cursor_file
        character(20)   :: cursor_filename
        character(10)   :: tmp_str

        write(tmp_str, "(i10)")ntime
        cursor_filename="cursor_position"//trim(adjustl(tmp_str))
        cursor_file=trim(working_dir)//"/"//trim(cursor_filename)

        if (verbose) then
            open(15,file="inoutflow2", position="append")
            if (nrank==0) write(15,*) "write cursor file at", trim(cursor_file)
            close(15)
        endif


        open(15, file=trim(cursor_file))
        write(15,*) cursor_mod
        write(15,*) cursor
        close(15)

    end subroutine write_cursorfile


    subroutine read_cursorfile(cursor_file)

        implicit none
        character(200)  :: cursor_file

        open(15, file=trim(cursor_file))
        read(15,*) cursor_mod
        read(15,*) cursor
        close(15)

    end subroutine read_cursorfile

end module flow_buffer_dao


module flow_buffer_handler

use flow_buffer_data
use flow_buffer_dao

contains

!
!   Open a session in READ or WRIDE mode (mode input argument).
!               n is the size of the buffers.
!               recovery_dir_arg is the directory that will be used to restore buffer data if continuing from previous run.
!               working_arg is the directory were the buffers will be stored or read when buffers limit size is exceeded.

    subroutine init(mode, n, recovery_dir_arg, working_dir_arg)

        use run_ctxt_data
        use start_settings

        implicit none
        integer, intent(in)     :: mode, n
        character(*), intent(in)    :: recovery_dir_arg, working_dir_arg

        character(200)  :: cursor_file
        character(20)   :: cursor_filename
        character(10)   :: tmp_str

        access_mode=mode
        nb_frames=n
        recovery_dir=trim(recovery_dir_arg)
        working_dir=trim(working_dir_arg)

        ! Load the cursor position  ---------------------------------------------------
        if ((run_ctxt==CONTINUE_FROM_PREVIOUS_RUN).or.(run_ctxt==RECOVERY_A_RUN)) then

            cursor_file=trim(recovery_dir)//"/cursor_position"

            call read_cursorfile(cursor_file)


        elseif (start_source_type==HDF5_FILE) then

            write(tmp_str, "(i10)")start_it
            cursor_filename="cursor_position"//trim(adjustl(tmp_str))
            cursor_file=trim(working_dir)//"/"//trim(cursor_filename)

            call read_cursorfile(cursor_file)

        else

            cursor_mod=1
            cursor=1

        end if


    end subroutine init

    subroutine display_buffer(id)

        use decomp_2d
        use COMMON_fieldstools

        implicit none
        integer, intent(in)        :: id
        real*8, dimension(nb_frames, xstart(2):xend(2), xstart(3):xend(3))    :: cur_buffer

        integer :: i
        real*8  :: cs


        open(15,file="inoutflow2", position="append")
        if (nrank==0)   write(15,*)"Affichage du buffer:", trim(buffer_list(id)%name)


        cur_buffer=buffer_list(id)%data(:,:,:)
        do i = 1, nb_frames

            call perform_checksum2D_x(cur_buffer(i, :,:), cs)

            if(i==cursor)then
                if (nrank==0)   write(15,*)"_", i, cs
            endif
            if(i/=cursor)then
                if (nrank==0)   write(15,*)"+", i, cs
            endif
        end do

        close(15)


    end subroutine display_buffer

    subroutine open_newbuffer(id, name)

        use decomp_2d
        use HDF5_IO

        use run_ctxt_data
        use start_settings

        implicit none
        integer, intent(out)        :: id
        character(*)                :: name


        character(200)  :: buffers_file
        real*8, dimension(nb_frames, xstart(2):xend(2), xstart(3):xend(3))    :: cur_buffer
        character(20)                                                           :: cur_buffer_name
        character(20)                                                           :: buffer_filename
        character(10)   :: tmp_str

        nb_buffer=nb_buffer+1

        allocate(buffer_list(nb_buffer)%data(nb_frames, xstart(2):xend(2), xstart(3):xend(3)) )
        buffer_list(nb_buffer)%name=trim(name)

        buffer_list(nb_buffer)%data=0.d0


        id=nb_buffer




        ! Load the cursor position  ---------------------------------------------------

        cur_buffer_name=trim(buffer_list(id)%name)

        if ((run_ctxt==CONTINUE_FROM_PREVIOUS_RUN).or.(run_ctxt==RECOVERY_A_RUN)) then

            buffers_file=trim(recovery_dir)//"/buffers"

            call load_buffer(id, buffers_file)



        elseif (start_source_type==HDF5_FILE) then

            write(tmp_str, "(i10)")cursor_mod
            buffer_filename="buffers_"//trim(adjustl(tmp_str))
            buffers_file=trim(working_dir)//"/"//trim(buffer_filename)


            if (verbose) then
                open(15,file="inoutflow2", position="append")
                if(nrank==0) write(15,*)"HDF5_FILE reading buffer from file ", trim(buffers_file)
                close(15)
            endif

            call load_buffer(id, buffers_file)


        end if

    end subroutine open_newbuffer

    subroutine close_allbuffers()
        implicit none
        integer     :: i

        do i = 1, nb_buffer
            deallocate(buffer_list(i)%data)
        end do

        nb_buffer=0
        cursor_mod=1
        cursor=1

        nb_frames=0

    end subroutine close_allbuffers

    subroutine flush_all_buffers(buffers_file)

        use decomp_2d

        implicit none
        character(200), intent(in)  :: buffers_file
        integer i

        call create_buffer_file(buffers_file)

        do i = 1, nb_buffer
            call flush_buffer(i, buffers_file)
        end do

    end subroutine flush_all_buffers

    subroutine save_all_buffers(buffers_file)

        use decomp_2d

        implicit none
        character(200), intent(in)  :: buffers_file
        integer i

        call create_buffer_file(buffers_file)

        do i = 1, nb_buffer
            call save_buffer(i, buffers_file)
        end do

    end subroutine save_all_buffers

    subroutine load_all_buffers(buffers_file)

        use decomp_2d

        implicit none
        character(200), intent(in)  :: buffers_file
        integer i

        do i = 1, nb_buffer
            call load_buffer(i, buffers_file)
        end do

    end subroutine load_all_buffers








    subroutine update_buffer(id, flow2d, action)

        use decomp_2d

        use HDF5_IO
        use run_ctxt_data

        implicit none
        integer, intent(in)                                                     :: action
        integer, intent(in)                                                     :: id
        real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)), intent(inout)  :: flow2d

        character(200)                                                          :: buffer_file
        character(10)                                                           :: tmp_str

        logical, save, dimension(MAX_BUFFER)                                    :: busy_buffers=.false.
        integer, save                                                           :: nb_busy_buffers=0




        ! Le cas où le buffer est enregistré comme occupé mais que les autres buffers n'ont pas tous été remplis
        ! pour la frame courante est un cas anormal qui correspond à une mauvaise utilisation du module. Ce cas
        ! provoque donc l'arret du programme.
        if ((nb_busy_buffers/=nb_buffer).and.(busy_buffers(id))) then

            call exit(1)

        endif


        ! On remplit le buffer identifié par "id" avec le champ donné en paramètre. A la fin de cette opération,
        ! on enregistre le buffer comme étant occupé et on incrémente le nombre de buffers occupés.
        if(action==WRITE_MODE) then
            buffer_list(id)%data(cursor,:,:)=flow2d

            if (verbose) then
                open(15,file="inoutflow2", position="append")
                if(nrank==0) write(15,*)ntime, 'w ', trim(buffer_list(id)%name), &
                            sum(abs(buffer_list(id)%data(cursor,:,:)))/(xsize(2)*xsize(3)), cursor, cursor_mod
                close(15)
            endif


        endif

        if(action==READ_MODE) then
            flow2d=buffer_list(id)%data(cursor,:,:)


            if (verbose) then
                open(15,file="inoutflow2", position="append")
                if(nrank==0) write(15,*)ntime, 'r ', trim(buffer_list(id)%name), &
                            sum(abs(flow2d))/(xsize(2)*xsize(3)), cursor, cursor_mod
                close(15)
            endif


        endif


        busy_buffers(id)=.true.
        nb_busy_buffers=nb_busy_buffers+1



        ! Si tous les buffers ont déjà été remplis pour la frame courante, on passe à la frame suivante.
        if (nb_busy_buffers==nb_buffer) then
            cursor=cursor+1

            busy_buffers=.false.
            nb_busy_buffers=0

            if (verbose) then
                open(15,file="inoutflow2", position="append")
                if(nrank==0) write(15,*)"NEW FRAME"
                close(15)
            endif

        endif



        ! Si le curseur dépasse le nombre de frames autorisées, on exporte les frames précédentes dans un fichier
        ! et on remet le curseur à zéro. Les anciennes frames sont donc effacées.
        if (cursor>nb_frames) then


            if (verbose) then
                open(15,file="inoutflow2", position="append")
                if(nrank==0) write(15,*)"Depassement: ", cursor,"/",nb_frames
                close(15)
            endif


            if(action==WRITE_MODE) then

                write(tmp_str, "(i10)")cursor_mod
                buffer_file=trim(working_dir)//'buffers_'//trim(adjustl(tmp_str))
                call flush_all_buffers(buffer_file)

            endif

            if(action==READ_MODE) then

                write(tmp_str, "(i10)")cursor_mod+1
                buffer_file=trim(working_dir)//'/buffers_'//trim(adjustl(tmp_str))

                if (verbose) then
                    open(15,file="inoutflow2", position="append")
                    if(nrank==0) write(15,*)"READ BUFFER FILE:  ", trim(buffer_file)
                    close(15)
                endif

                !buffer_list(id)%data(:,:,:)=-1.d0
                call load_all_buffers(buffer_file)
!
!                buffer_list(1)%data(:,:,:)=-1.d0
!                buffer_list(2)%data(:,:,:)=-2.d0
!                buffer_list(3)%data(:,:,:)=-3.d0
                if (verbose) call display_buffer(1)
            endif


            cursor_mod=cursor_mod+1
            cursor=1

        end if



    end subroutine update_buffer


!   Store the input 2d field in the buffer identified by "id".
    subroutine add_2dfield(id, flow2d)

        use decomp_2d

        implicit none
        integer, intent(in)                                                     :: id
        real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)), intent(inout)  :: flow2d

        call update_buffer(id, flow2d, WRITE_MODE)

    end subroutine add_2dfield

    ! Fill the input 2d field from the buffer identified by "id".
    subroutine get_2dfield(id, flow2d)

        use decomp_2d

        implicit none
        integer, intent(in)                                                     :: id
        real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)), intent(inout)  :: flow2d

        call update_buffer(id, flow2d, READ_MODE)

    end subroutine get_2dfield

    subroutine save_state()
        implicit none
        character(200)      :: buffers_file
        character(200)      :: cursor_file

        buffers_file=trim(recovery_dir)//"/buffers"

        call save_all_buffers(buffers_file)

        cursor_file=trim(recovery_dir)//"/cursor_position"

        open(15, file=trim(cursor_file))
        write(15,*) cursor_mod
        write(15,*) cursor
        close(15)

    end subroutine save_state

end module flow_buffer_handler

module flow_buffer_handler_test

    use COMMON_workspace_view

contains

    ! Test if the buffer are correctly restored in case where the simulation is performed from a previous run.
    subroutine recovery_test()


        use run_ctxt_data

        use mpi
        use decomp_2d

        use flow_buffer_handler
        use flow_buffer_data, only: verbose
        implicit none


        real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: linear, linear2, linear3
        integer :: id_lin, id_lin2, id_lin3
        integer :: outflow_buff
        integer::sb
        integer:: sbmax=3

        verbose=.true.

        outflow_buff=10

        call close_allbuffers
        call init(WRITE_MODE, outflow_buff, COMMON_recovery_outflow_dir, COMMON_outflow_path)

        call open_newbuffer(id_lin, "linear1")
        call open_newbuffer(id_lin2, "linear2")
        call open_newbuffer(id_lin3, "linear3")

        do ntime = 1, 23

            linear=ntime*1.d0
            linear2=ntime*2.d0
            linear3=ntime*3.d0

            do sb = 1, sbmax
                call add_2dfield(id_lin, linear)
                call add_2dfield(id_lin2, linear2)
                call add_2dfield(id_lin3, linear3)
            end do

        end do



        open(15,file="inoutflow2", position="append")
        if(nrank==0) write(15,*) "Affichage des buffers JUSTE AVANT save_state"
        close(15)
        call display_buffer(id_lin)
        call display_buffer(id_lin2)
        call display_buffer(id_lin3)

        call save_state
        call close_allbuffers




        run_ctxt=CONTINUE_FROM_PREVIOUS_RUN

        call init(WRITE_MODE, outflow_buff, COMMON_recovery_outflow_dir, COMMON_outflow_path)

        call open_newbuffer(id_lin, "linear1")
        call open_newbuffer(id_lin2, "linear2")
        call open_newbuffer(id_lin3, "linear3")

        open(15,file="inoutflow2", position="append")
        if(nrank==0) write(15,*) "Affichage des buffers JUSTE APRES ouverture depuis recovery"
        close(15)
        call display_buffer(id_lin)
        call display_buffer(id_lin2)
        call display_buffer(id_lin3)

!
!        do ntime = 24, 27
!
!            linear=ntime*1.d0
!            linear2=ntime*2.d0
!            linear3=ntime*3.d0
!
!
!            do sb = 1, sbmax
!                call add_2dfield(id_lin, linear)
!                call add_2dfield(id_lin2, linear2)
!                call add_2dfield(id_lin3, linear3)
!            end do
!
!        end do
!
!        call display_buffer(id_lin)
!        call display_buffer(id_lin2)
!        call display_buffer(id_lin3)
!
!
!        call close_allbuffers


    end subroutine recovery_test

    ! 2 step test:
    ! - 1st step: use of the buffer in the write mode to store 2D fields. At the end of this step, the buffers are
    !   expected to be stored in outflow directory
    ! - 2nd step: use of the buffer in the read mode to read 2d fields. The buffer are filled from previously generated
    !   buffers files and the 2d fields are obtained by reading the resulting buffer.

    ! This test check the correct communication between the write and read mode.
    subroutine inflow_outflow_communication_test()


        use run_ctxt_data
        use start_settings

        use mpi
        use decomp_2d

        use flow_buffer_handler
        use flow_buffer_data, only: verbose

        implicit none


        real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: linear, linear2, linear3
        integer :: id_lin, id_lin2, id_lin3
        integer :: outflow_buff
        integer::sb
        integer:: sbmax=3

        character(400)  :: command

        integer :: mpi_err


        verbose=.true.

        open(15,file="inoutflow2", position="append")
        if(nrank==0) write(15,*) "***********************************************************"
        if(nrank==0) write(15,*) "***********************************************************"
        if(nrank==0) write(15,*) "***********************************************************"
        if(nrank==0) write(15,*) "***********************************************************"
        if(nrank==0) write(15,*) "***********************************************************"
        if(nrank==0) write(15,*) "***********************************************************"
        close(15)

        outflow_buff=10

        run_ctxt=NEW_SIMULATION


        call close_allbuffers
        call init(WRITE_MODE, outflow_buff, COMMON_recovery_outflow_dir, COMMON_outflow_path)

        call open_newbuffer(id_lin, "linear1")
        call open_newbuffer(id_lin2, "linear2")
        call open_newbuffer(id_lin3, "linear3")

        do ntime = 1, 50

            linear=ntime*1.d0
            linear2=ntime*2.d0
            linear3=ntime*3.d0

            do sb = 1, sbmax
                call add_2dfield(id_lin, linear)
                call add_2dfield(id_lin2, linear2)
                call add_2dfield(id_lin3, linear3)
            end do

            if (mod(ntime, 8).eq.0) then


                open(15,file="inoutflow2", position="append")
                if(nrank==0) write(15,*) "Affichage des buffers JUSTE AVANT write_cursorfile"
                close(15)

                call display_buffer(id_lin)
                call display_buffer(id_lin2)
                call display_buffer(id_lin3)

                call write_cursorfile
            endif

        end do

        call close_allbuffers


        command="cp -rf "//trim(COMMON_outflow_path)//" "//trim(COMMON_results_path)//"/Inflow"
        if (nrank==0) write(*,*) command
        if (nrank==0) call system("cp -rf "//trim(COMMON_outflow_path)//" "//trim(COMMON_results_path)//"/Inflow")

        call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

        open(15,file="inoutflow2", position="append")
        if(nrank==0) write(15,*) "___________________________________________________________"
        if(nrank==0) write(15,*) "___________________________________________________________"
        if(nrank==0) write(15,*) "___________________________________________________________"
        close(15)



        run_ctxt=NEW_SIMULATION
        COMMON_inflow_path=trim(COMMON_results_path)//"Inflow"
        start_source_type=HDF5_FILE
        start_it=8

        call init(WRITE_MODE, outflow_buff, COMMON_recovery_outflow_dir, COMMON_inflow_path)

        call open_newbuffer(id_lin, "linear1")
        call open_newbuffer(id_lin2, "linear2")
        call open_newbuffer(id_lin3, "linear3")


        open(15,file="inoutflow2", position="append")
        if(nrank==0) write(15,*) "Affichage des buffers JUSTE APRES ouverture des buffers"
        close(15)

        call display_buffer(id_lin)
        call display_buffer(id_lin2)
        call display_buffer(id_lin3)


        start_it=start_it+1
        do ntime = start_it, 50

            linear=ntime*1.d0
            linear2=ntime*2.d0
            linear3=ntime*3.d0

            do sb = 1, sbmax
                call get_2dfield(id_lin, linear)
                call get_2dfield(id_lin2, linear2)
                call get_2dfield(id_lin3, linear3)
            end do

        end do

        open(15,file="inoutflow2", position="append")
        if(nrank==0) write(15,*) "***********************************************************"
        close(15)

    end subroutine inflow_outflow_communication_test

end module flow_buffer_handler_test
