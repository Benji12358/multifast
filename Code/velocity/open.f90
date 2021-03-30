! WARNING !!!!!!!!!
! CHANGE IMPLEMENTATION USE THE COMMON/UTILS/buffers2d MODULE !!!
! HERE WE ARE IN VELOCITY MODULE, AND WE SHOULD TREAT SCALAR AS IT IS DONE HERE !!!!!!!!!!!!!!!

module VELOCITY_inout_flow
    use physical_fields
    use DNS_settings
    use inflow_settings
    use mesh
    use decomp_2d
    use boundaries
    implicit none

    integer :: id_q1, id_q2, id_q3

contains

    subroutine inoutflow_init()

        use flow_buffer_handler

        call open_newbuffer(id_q1, "q1")
        call open_newbuffer(id_q2, "q2")
        call open_newbuffer(id_q3, "q3")

        call display_buffer(id_q1)


    end subroutine

    subroutine add_outflow()

        use flow_buffer_handler
        use run_ctxt_data
        use COMMON_fieldstools

        implicit none
        real*8  :: cs1, cs2, cs3

        call add_2dfield(id_q1, q1_x(n1-1, :,:))
        call add_2dfield(id_q2, q2_x(n1-1, :,:))
        call add_2dfield(id_q3, q3_x(n1-1, :,:))

        call perform_checksum2D_x(q1_x(n1-1, :,:), cs1)
        call perform_checksum2D_x(q2_x(n1-1, :,:), cs2)
        call perform_checksum2D_x(q3_x(n1-1, :,:), cs3)

        if (nrank==0) then
            open(15, file="debugoutflow", position="append")
            write(15,*) ntime, cs1, cs2, cs3
            close(15)
        endif

    end subroutine add_outflow

    subroutine get_inflow()
        use HDF5_IO
        use inflow_settings

        implicit none
        integer, save   :: inflow_ptr=1
        integer, save   :: inflow_nb=0
        integer         :: inflow_nb2
        character(200)       :: current_inflow_path

        if (inflow_buff>0) call inflow_from_file
        if (inflow_buff==0) call default_inflow


        return

        contains

            subroutine default_inflow()
                implicit none
                integer j

                if (nrank==0) write(*,*)'PERFORMING default poiseuille inflow'

                !if (nrank==0) open(15, file="Poiseuille.csv")
                ! do j = xstart(2), min(xend(2), n2-1)
                !     q1_x(1,j,:)=1.d0!-(Yc(j)-1.d0)**2
                !     !if (nrank==0) write(15,*) Yc(j), 1.d0-(Yc(j)-1.d0)**2
                ! end do

                !if (nrank==0) close(15)

                if (inflow_mode==BOUNDARY_LAYER) then
                    q1_x(1,:,:)=0.d0
                    q2_x(1,:,:)=0.d0
                    q3_x(1,:,:)=0.d0
                    if (streamwise==1) then
                        call perform_boundary_layer_1(q1_x(1,:,:), delta_BL)
                        q1_wall10(:,:)=q1_x(1,:,:)
                    endif

                    if (streamwise==3) then
                        call perform_boundary_layer_3(q3_x(1,:,:), delta_BL)
                        q3_wall30(:,:)=q3_x(1,:,:)
                    endif

                endif

                if (inflow_mode==CHANNEL_FLOW) then
                    q1_x(1,:,:)=0.d0
                    q2_x(1,:,:)=0.d0
                    q3_x(1,:,:)=0.d0
                    if (streamwise==1) then
                        call perform_stream1(q1_x(1,:,:), BC2, BC3)
                        q1_wall10(:,:)=q1_x(1,:,:)
                    endif

                    if (streamwise==3) then
                        call perform_stream3(q3_x(1,:,:), BC2, BC1)
                        q3_wall30(:,:)=q3_x(1,:,:)
                    endif

                endif                

                if (inflow_mode==CONSTANT_FLOW) then
                    q1_x(1,:,:)=0.d0
                    q2_x(1,:,:)=0.d0
                    q3_x(1,:,:)=0.d0
                    if (streamwise==1) then
                        q1_x(1,:,:)=1.d0
                        q1_wall10(:,:)=q1_x(1,:,:)
                    endif

                    if (streamwise==3) then
                        q3_x(1,:,:)=1.d0
                        q3_wall30(:,:)=q3_x(1,:,:)
                    endif
                endif

                ! call perform_stream1(q1_x(1,:,:), BC2, BC3)


            end subroutine default_inflow

            subroutine inflow_from_file()
                use run_ctxt_data
                use flow_buffer_handler
                use COMMON_fieldstools

                implicit none
                real*8  :: cs1, cs2, cs3

                call get_2dfield(id_q1, q1_x(1, :,:))
                call get_2dfield(id_q2, q2_x(1, :,:))
                call get_2dfield(id_q3, q3_x(1, :,:))

                call perform_checksum2D_x(q1_x(1, :,:), cs1)
                call perform_checksum2D_x(q2_x(1, :,:), cs2)
                call perform_checksum2D_x(q3_x(1, :,:), cs3)

                if (nrank==0) then
                    open(15, file="debugoutflow", position="append")
                    write(15,*) ntime, cs1, cs2, cs3
                    close(15)
                endif


            end subroutine inflow_from_file

            subroutine perform_stream1(stream1, BC2, BC3)
              implicit none

              real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: stream1
              integer                                                 :: BC2, BC3

              real*8                                                  :: f2(n2), f3(n3), c2, c3
              integer                                                 :: j, k

              f3=1.d0

              if (mod(n2, 2)==0) then
                  c2=Yc(n2/2)
              else
                  c2=Y((n2-1)/2+1)
              end if
              

              if (BC2==NOSLIP) then

                do j=1,n2-1
                    f2(j)=(1.d0-((Yc(j)-c2)/(0.5d0*L2))**2)
                enddo
                f2(n2)=f2(n2-1)

                  if (BC3==NOSLIP) then

                      do k=1,n3

                          if (mod(n3, 2)==0) then
                              c3=(n3/2-0.5d0)*dx3
                          else
                              c3=((n3-1)/2+1-0.5d0)*dx3
                          end if

                          f3(k)=(1.d0-(((k-0.5d0)*dx3-c3)/(0.5d0*L3))**2)

                      enddo

                  end if

                  do k = xstart(3),xend(3)
                      do j= xstart(2),xend(2)
                          stream1(j,k)=f2(j)*f3(k)
                      end do
                  end do

              else
                  stream1=0.d0
              end if

            end subroutine perform_stream1

            subroutine perform_stream3(stream3, BC1, BC2)
                implicit none

                real*8, dimension(zstart(1):zend(1), zstart(2):zend(2)) :: stream3
                integer                                                 :: BC1, BC2

                real*8                                                  :: f1(n1), f2(n2), c1, c2
                integer                                                 :: i, j,temp_test

                f1=1.d0


                if (BC2==NOSLIP) then

                    if (mod(n2, 2)==0) then
                       c2=Yc(n2/2)
                    else
                        c2=Y((n2-1)/2+1)
                    end if

                    do j=1,n2-1

                        f2(j)=(1.d0-((Yc(j)-c2)/(0.5d0*L2))**2)

                    enddo

                   if (BC1==NOSLIP) then

                        if (mod(n1, 2)==0) then
                           c1=(n1/2-0.5d0)*dx1
                        else
                          c1=((n1-1)/2+1-0.5d0)*dx1
                        end if

                        f1(i)=(1.d0-(((i-0.5d0)*dx1-c1)/(0.5d0*L1))**2)
                        
                     else !

                        do j = zstart(2),zend(2)
                         do i= zstart(1),zend(1)
                            stream3(i,j)=f2(j)*f1(i)
                         end do
                       end do

                   endif !BC1

                else !
                    stream3=0.d0
                end if !BC2

            end subroutine perform_stream3

            subroutine perform_boundary_layer_1(stream1, delta_BL)
            implicit none

            real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: stream1
            real*8, dimension(1:n2)                                 :: f1

            real*8                                                  :: delta_BL
            integer                                                 :: j
            logical                                                 :: pair_n2

            pair_n2 = (mod(n2, 2)==0)

            do j=1,n2/2
                if (Yc(j)<delta_BL) then
                    f1(j) = 1.5d0 * (Yc(j)/delta_BL) - 0.5d0 * (Yc(j)/delta_BL)**2
                    f1(n2-j) = 1.5d0 * (Yc(j)/delta_BL) - 0.5d0 * (Yc(j)/delta_BL)**2
                else
                    f1(j) = 1.d0
                    f1(n2-j) = 1.d0
                endif
            enddo

            if (pair_n2) then
                f1(n2/2+1) = 1.d0
            endif

            do j = xstart(2),xend(2)
                stream1(j,:)=f1(j)
            end do

            end subroutine perform_boundary_layer_1

            subroutine perform_boundary_layer_3(stream3, delta_BL)
                implicit none

                real*8, dimension(zstart(1):zend(1), zstart(2):zend(2)) :: stream3
                real*8, dimension(1:n2)                                 :: f3

                real*8                                                  :: delta_BL
                integer                                                 :: j
                logical                                                 :: pair_n2

                pair_n2 = (mod(n2, 2)==0)

                do j=1,n2/2
                    if (Yc(j)<delta_BL) then
                        f3(j) = (3/2) * (Yc(j)/delta_BL) - (1/2) * (Yc(j)/delta_BL)**2
                        f3(n2-j) = (3/2) * (Yc(j)/delta_BL) - (1/2) * (Yc(j)/delta_BL)**2
                    else
                        f3(j) = 1.d0
                        f3(n2-j) = 1.d0
                    endif
                enddo

                if (pair_n2) then
                    f3(n2/2+1) = 1.d0
                endif

                do j = zstart(2),zend(2)
                    stream3(:,j)=f3(j)
                end do

            end subroutine perform_boundary_layer_3


    end subroutine get_inflow

end module VELOCITY_inout_flow


! WARNING !!!!!!!!!
! CHANGE IMPLEMENTATION USE THE COMMON/UTILS/buffers2d MODULE !!!
! HERE WE ARE IN VELOCITY MODULE, AND WE SHOULD TREAT SCALAR AS IT IS DONE HERE !!!!!!!!!!!!!!!

module VELOCITY_inout_flow_old
    use physical_fields
    use DNS_settings
    use inflow_settings
    use mesh
    use decomp_2d
    implicit none

    integer :: ptr_size

    real*8,dimension (:,:,:,:), allocatable       :: q1_outflow, q2_outflow, q3_outflow
    real*8,dimension (:,:,:,:), allocatable       :: q1_inflow, q2_inflow, q3_inflow

contains

    subroutine inoutflow_init(nb_substep, out_buffer_sz, in_buffer_sz)
        implicit none
        integer, intent(in) :: nb_substep, out_buffer_sz, in_buffer_sz

        ptr_size=nb_substep

        if (outflow_buff>0) then
            allocate(q1_outflow(nb_substep, out_buffer_sz, xstart(2):xend(2), xstart(3):xend(3)))
            allocate(q2_outflow(nb_substep, out_buffer_sz, xstart(2):xend(2), xstart(3):xend(3)))
            allocate(q3_outflow(nb_substep, out_buffer_sz, xstart(2):xend(2), xstart(3):xend(3)))
!            allocate(sca_outflow(nb_substep, out_buffer_sz, xstart(2):xend(2), xstart(3):xend(3)))
        end if

        if (inflow_buff>0) then
            allocate(q1_inflow(nb_substep, in_buffer_sz, xstart(2):xend(2), xstart(3):xend(3)))
            allocate(q2_inflow(nb_substep, in_buffer_sz, xstart(2):xend(2), xstart(3):xend(3)))
            allocate(q3_inflow(nb_substep, in_buffer_sz, xstart(2):xend(2), xstart(3):xend(3)))
!            allocate(sca_inflow(nb_substep, in_buffer_sz, xstart(2):xend(2), xstart(3):xend(3)))
        end if

    end subroutine

    subroutine add_outflow(ntime, ns)
        use HDF5_IO
        use COMMON_workspace_view

        implicit none
        integer, optional   :: ntime, ns
        integer, save       :: outflow_ptr=0
        character(200)      :: current_outflow_path

        real*8, save       :: flowmark=0.d0

        character*10 tmp_str, tmp_str2

        flowmark=flowmark+1.d0

        if (ns==1) outflow_ptr=outflow_ptr+1

        q1_outflow(ns,outflow_ptr, :,:)=q1_x(n1-1, :,:)
!        q1_outflow(ns,outflow_ptr, :,:)=flowmark!q1_x(n1-1, :,:)
        q2_outflow(ns,outflow_ptr, :,:)=q2_x(n1-1, :,:)
        q3_outflow(ns,outflow_ptr, :,:)=q3_x(n1-1, :,:)
!        sca_outflow(ns,outflow_ptr, :,:)=sca_x(n1-1, :,:)


!        open(15,file="inoutflow", position="append")
!        if(nrank==0) write(15,*)ntime, sum(q1_outflow(ns,outflow_ptr, :,:))/(xsize(2)*xsize(3))
!        close(15)

        if (nrank==0) write(*,*)
        if (nrank==0) write(*,*) "Add outflow at slice n=", outflow_ptr, "(ntime=", ntime,", ns=", ns, ")"


        if (mod(ntime+1, outflow_buff)==0) then

            if (ns==ptr_size) outflow_ptr=0


            write(tmp_str, "(i10)")ntime-outflow_buff+1
            write(tmp_str2, "(i10)")ns
            current_outflow_path=trim(COMMON_outflow_path)//'outflow_'//trim(adjustl(tmp_str))//"_"//trim(adjustl(tmp_str2))

            if (nrank==0) write(*,*)
            if (nrank==0) write(*,*) "Writting outflow in file: ", trim(current_outflow_path)//".h5"

            if(nrank==0)  call hdf_create_emptyfile(current_outflow_path)
            call hdf_add_3Dfield(current_outflow_path, q1_outflow(ns,:,:,:), "q1_out", outflow_buff, ny_global, nz_global, 1,outflow_buff, xstart(2),xend(2), xstart(3),xend(3))
            call hdf_add_3Dfield(current_outflow_path, q2_outflow(ns,:,:,:), "q2_out", outflow_buff, ny_global, nz_global, 1,outflow_buff, xstart(2),xend(2), xstart(3),xend(3))
            call hdf_add_3Dfield(current_outflow_path, q3_outflow(ns,:,:,:), "q3_out", outflow_buff, ny_global, nz_global, 1,outflow_buff, xstart(2),xend(2), xstart(3),xend(3))
!            call hdf_add_3Dfield(current_outflow_path, sca_outflow(ns,:,:,:), "sca_out", outflow_buff, ny_global, nz_global, 1,outflow_buff, xstart(2),xend(2), xstart(3),xend(3))

        endif

    end subroutine add_outflow

    subroutine get_inflow(ntime, ns)
        use HDF5_IO
        use inflow_settings

        implicit none
        integer         :: ntime, ns
        integer, save   :: inflow_ptr=1
        integer, save   :: inflow_nb=0
        integer         :: inflow_nb2
        character(200)       :: current_inflow_path

        if (inflow_buff>0) call inflow_from_file
        if (inflow_buff==0) call default_inflow


        return

        contains

            subroutine default_inflow()
                implicit none
                integer j

                if (nrank==0) write(*,*)'PERFORMING default inflow'

                !if (nrank==0) open(15, file="Poiseuille.csv")
                do j = xstart(2), min(xend(2), n2-1)
                    q1_x(1,j,:)=1.d0!-(Yc(j)-1.d0)**2
                    !if (nrank==0) write(15,*) Yc(j), 1.d0-(Yc(j)-1.d0)**2
                end do

                !if (nrank==0) close(15)



                q1_wall10(:,:)=q1_x(1, :,:)
                q2_x(1, :,:)=0.d0
                q3_x(1, :,:)=0.d0
                !sca_x(1, :,:)=0.d0


            end subroutine default_inflow

            subroutine inflow_from_file()

                use COMMON_workspace_view
                use start_settings, only:start_it

                implicit none

                character*10 tmp_str, tmp_str2

                open(15, file="debugoutflow", position="append")
                write(15,*)"OLD OUTFLOW VELOCITY"
                close(15)

                if (inflow_ptr==1) then

                    inflow_nb2=inflow_nb*inflow_buff+start_it
                    write(tmp_str, "(i10)")inflow_nb2
                    write(tmp_str2, "(i10)")ns
                    write(*,*)'ns:', ns, tmp_str2
                    current_inflow_path=trim(COMMON_inflow_path)//'outflow_'//trim(adjustl(tmp_str))//"_"//trim(adjustl(tmp_str2))

                    if (nrank==0) write(*,*)
                    if (nrank==0) write(*,*) "Reading inflow from file:", trim(current_inflow_path)//".h5"

                    call hdf_read_3Dfield(current_inflow_path, q1_inflow(ns,:,:,:), "q1_out", inflow_buff, ny_global, nz_global, 1,inflow_buff, xstart(2),xend(2), xstart(3),xend(3))
                    call hdf_read_3Dfield(current_inflow_path, q2_inflow(ns,:,:,:), "q2_out", inflow_buff, ny_global, nz_global, 1,inflow_buff, xstart(2),xend(2), xstart(3),xend(3))
                    call hdf_read_3Dfield(current_inflow_path, q3_inflow(ns,:,:,:), "q3_out", inflow_buff, ny_global, nz_global, 1,inflow_buff, xstart(2),xend(2), xstart(3),xend(3))
!                    call hdf_read_3Dfield(current_inflow_path, sca_inflow(ns,:,:,:), "sca_out", inflow_buff, ny_global, nz_global, 1,inflow_buff, xstart(2),xend(2), xstart(3),xend(3))

                    if (ns==ptr_size) inflow_nb=inflow_nb+1

                endif

                if (inflow_ptr<=inflow_buff) then

                    if (nrank==0) write(*,*)
                    if (nrank==0) write(*,*) "Get inflow n=", inflow_ptr, "(ntime=", ntime,", ns=", ns, ")"

                    q1_x(1, :,:)=q1_inflow(ns, inflow_ptr, :,:)
                    q2_x(1, :,:)=q2_inflow(ns, inflow_ptr, :,:)
                    q3_x(1, :,:)=q3_inflow(ns, inflow_ptr, :,:)
!                    sca_x(1,1, :,:)=sca_inflow(ns, inflow_ptr, :,:)

                    q1_wall10(:,:)=q1_x(1, :,:)

                    if (ns==ptr_size) inflow_ptr=inflow_ptr+1








                    !q1_x(1, :,:)=1.d0
!                    open(15,file="inoutflow", position="append")
!                    if(nrank==0) write(15,*)ntime, sum(abs(q1_x(1,:,:)))/(xsize(2)*xsize(3))
!                    close(15)
!                    call default_inflow

                endif

                if (inflow_ptr>inflow_buff) inflow_ptr=1

            end subroutine inflow_from_file


    end subroutine get_inflow

end module VELOCITY_inout_flow_old





