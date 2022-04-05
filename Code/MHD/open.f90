!
!! TO TEST !!! NOT EVER USED !!!
!module SCALAR_inout_flow
!    use SCALAR_data
!    use DNS_settings
!    use inflow_settings
!    use mesh
!    use decomp_2d
!    implicit none
!
!    integer :: id_sca, id_lin
!    logical :: inout_newversion=.true.
!
!    real*8, dimension(:,:), allocatable:: linear
!
!contains
!
!    subroutine inoutflow_init()
!
!        use flow_buffer_handler
!
!        allocate(linear(xstart(2):xend(2), xstart(3):xend(3)))
!
!        call open_newbuffer(id_sca, "scalar")
!
!        ! For tests purpose !!
!        call open_newbuffer(id_lin, "linear")
!
!    end subroutine
!
!    subroutine add_outflow()
!
!        use flow_buffer_handler
!        use run_ctxt_data
!        use COMMON_fieldstools
!
!        implicit none
!        real*8  :: cs1, cs2
!
!        linear=ntime*1.d0
!        call add_2dfield(id_sca, sca_x(n1-1, :,:))
!        call add_2dfield(id_lin, linear)
!
!        call perform_checksum2D_x(linear(:,:), cs1)
!        call perform_checksum2D_x(sca_x(n1-1, :,:), cs2)
!
!        if (nrank==0) then
!            open(15, file="debugoutflow_sca", position="append")
!            write(15,*) ntime, cs1, cs2
!            close(15)
!        endif
!
!    end subroutine add_outflow
!
!    subroutine get_inflow()
!
!        use inflow_settings
!        use flow_buffer_handler
!
!        implicit none
!        integer, save   :: inflow_ptr=1
!        integer, save   :: inflow_nb=0
!        integer         :: inflow_nb2
!        character(200)       :: current_inflow_path
!
!        if (inflow_buff>0) call inflow_from_file
!        if (inflow_buff==0) call default_inflow
!
!
!        return
!
!        contains
!
!            subroutine default_inflow()
!                implicit none
!                integer j
!
!            end subroutine default_inflow
!
!
!
!            subroutine inflow_from_file()
!
!                use flow_buffer_handler
!                use run_ctxt_data
!                use COMMON_fieldstools
!
!                implicit none
!                real*8  :: cs1, cs2
!
!                call get_2dfield(id_sca, sca_x(1, :,:))
!                call get_2dfield(id_lin, linear)
!
!                call perform_checksum2D_x(linear(:,:), cs1)
!                call perform_checksum2D_x(sca_x(1, :,:), cs2)
!
!                if (nrank==0) then
!                    open(15, file="debugoutflow_sca", position="append")
!                    write(15,*) ntime, cs1, cs2
!                    close(15)
!                endif
!
!            end subroutine inflow_from_file
!
!
!    end subroutine get_inflow
!
!end module SCALAR_inout_flow
!
!
!module SCALAR_inout_flow_old
!    use SCALAR_data
!    use DNS_settings
!    use inflow_settings
!    use mesh
!    use decomp_2d
!    implicit none
!
!    integer :: ptr_size
!
!    real*8,dimension (:,:,:,:), allocatable       :: q1_outflow, q2_outflow, q3_outflow
!    real*8,dimension (:,:,:,:), allocatable       :: q1_inflow, q2_inflow, q3_inflow
!
!    logical :: inout_newversion=.false.
!
!contains
!
!    subroutine inoutflow_init(nb_substep, out_buffer_sz, in_buffer_sz)
!        implicit none
!        integer, intent(in) :: nb_substep, out_buffer_sz, in_buffer_sz
!
!        ptr_size=nb_substep
!
!        if (outflow_buff>0) then
!            allocate(sca_outflow(nb_substep, out_buffer_sz, xstart(2):xend(2), xstart(3):xend(3)))
!        end if
!
!        if (inflow_buff>0) then
!            allocate(sca_inflow(nb_substep, in_buffer_sz, xstart(2):xend(2), xstart(3):xend(3)))
!        end if
!
!    end subroutine
!
!    subroutine add_outflow(ntime, ns)
!        use HDF5_IO
!        use COMMON_workspace_view
!
!        implicit none
!        integer, optional   :: ntime, ns
!        integer, save       :: outflow_ptr=0
!        character(200)      :: current_outflow_path
!
!        real*8, save       :: flowmark=0.d0
!
!        character*10 tmp_str, tmp_str2
!
!        flowmark=flowmark+1.d0
!
!        if (ns==1) outflow_ptr=outflow_ptr+1
!
!        sca_outflow(ns,outflow_ptr, :,:)=sca_x(n1-1, :,:)
!
!
!!        open(15,file="inoutflow", position="append")
!!        if(nrank==0) write(15,*)ntime, sum(q1_outflow(ns,outflow_ptr, :,:))/(xsize(2)*xsize(3))
!!        close(15)
!
!        if (nrank==0) write(*,*)
!        if (nrank==0) write(*,*) "Add outflow at slice n=", outflow_ptr, "(ntime=", ntime,", ns=", ns, ")"
!
!
!        if (mod(ntime+1, outflow_buff)==0) then
!
!            if (ns==ptr_size) outflow_ptr=0
!
!
!            write(tmp_str, "(i10)")ntime-outflow_buff+1
!            write(tmp_str2, "(i10)")ns
!            current_outflow_path=trim(COMMON_outflow_path)//'outflow_'//trim(adjustl(tmp_str))//"_"//trim(adjustl(tmp_str2))
!
!            if (nrank==0) write(*,*)
!            if (nrank==0) write(*,*) "Writting outflow in file: ", trim(current_outflow_path)//".h5"
!
!            call hdf_add_3Dfield(current_outflow_path, sca_outflow(ns,:,:,:), "sca_out", outflow_buff, ny_global, nz_global, 1,outflow_buff, xstart(2),xend(2), xstart(3),xend(3))
!
!        endif
!
!    end subroutine add_outflow
!
!    subroutine get_inflow(ntime, ns)
!        use HDF5_IO
!        use inflow_settings
!
!        implicit none
!        integer         :: ntime, ns
!        integer, save   :: inflow_ptr=1
!        integer, save   :: inflow_nb=0
!        integer         :: inflow_nb2
!        character(200)       :: current_inflow_path
!
!        if (inflow_buff>0) call inflow_from_file
!        if (inflow_buff==0) call default_inflow
!
!
!        return
!
!        contains
!
!            subroutine default_inflow()
!                implicit none
!                integer j
!
!            end subroutine default_inflow
!
!            subroutine inflow_from_file()
!
!                use COMMON_workspace_view
!                use start_settings, only:start_it
!
!                implicit none
!
!                character*10 tmp_str, tmp_str2
!
!                open(15, file="debugoutflow", position="append")
!                write(15,*)"OLD OUTFLOW VELOCITY"
!                close(15)
!
!                if (inflow_ptr==1) then
!
!                    inflow_nb2=inflow_nb*inflow_buff+start_it
!                    write(tmp_str, "(i10)")inflow_nb2
!                    write(tmp_str2, "(i10)")ns
!                    write(*,*)'ns:', ns, tmp_str2
!                    current_inflow_path=trim(COMMON_inflow_path)//'outflow_'//trim(adjustl(tmp_str))//"_"//trim(adjustl(tmp_str2))
!
!                    if (nrank==0) write(*,*)
!                    if (nrank==0) write(*,*) "Reading inflow from file:", trim(current_inflow_path)//".h5"
!
!                    call hdf_read_3Dfield(current_inflow_path, sca_inflow(ns,:,:,:), "sca_out", inflow_buff, ny_global, nz_global, 1,inflow_buff, xstart(2),xend(2), xstart(3),xend(3))
!
!                    if (ns==ptr_size) inflow_nb=inflow_nb+1
!
!                endif
!
!                if (inflow_ptr<=inflow_buff) then
!
!                    if (nrank==0) write(*,*)
!                    if (nrank==0) write(*,*) "Get inflow n=", inflow_ptr, "(ntime=", ntime,", ns=", ns, ")"
!
!                    sca_x(1, :,:)=sca_inflow(ns, inflow_ptr, :,:)
!
!                    if (ns==ptr_size) inflow_ptr=inflow_ptr+1
!
!
!                    !q1_x(1, :,:)=1.d0
!!                    open(15,file="inoutflow", position="append")
!!                    if(nrank==0) write(15,*)ntime, sum(abs(q1_x(1,:,:)))/(xsize(2)*xsize(3))
!!                    close(15)
!!                    call default_inflow
!
!                endif
!
!                if (inflow_ptr>inflow_buff) inflow_ptr=1
!
!            end subroutine inflow_from_file
!
!
!    end subroutine get_inflow
!
!end module SCALAR_inout_flow_old
!
!
