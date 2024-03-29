module CORE_IO_settings

    use COMMON_workspace_view

    implicit none

contains

        subroutine read_start_settings()


            use start_settings

            implicit none

            integer             :: reading_format
            integer start_from_coarse_file_int, IBM_activated_int, generic_poisson_flag

            ! Start **********************************************************************************
            open(15,file=trim(COMMON_settings_path)//'start.d')

            read(15,*) start_source_type
            read(15,*) start_from_coarse_file_int, n1c, n2c, n3c !
            read(15,*)
            read(15,'(a)') external_fields_path  ! path of the read field  file (cha.rea)
            read(15,*) start_it
            read(15,*) vper
            read(15,*) reading_format
            close(15)

            if (start_from_coarse_file_int==1) then

                if (start_source_type==HDF5_FILE) then
                    start_from_coarse_file=.true.
                else
                    write(*,*) "La reprise a partir d'un maillage grossier n'est possible qu'a partir d'un fichier 2Decomp"
                    call exit
                end if

            else
                start_from_coarse_file=.false.
            end if

        end subroutine read_start_settings


        subroutine read_discretization_settings()

            use numerical_methods_settings
            use DNS_settings

            implicit none
            integer     :: generic_poisson_flag

            open(15,file=trim(COMMON_settings_path)//'discretization.d')

            read(15,*)schemes_configuration
            read(15,*) time_advancement, dt
            read(15,*) generic_poisson_flag

            use_generic_poisson=(generic_poisson_flag==1)

            close(15)

        end subroutine read_discretization_settings

        subroutine read_domain_settings()

            use mathematical_constants

            use mesh

            implicit none

            ! Computational domain *******************************************************************
            open(15,file=trim(COMMON_settings_path)//'computational_domain.d')

            read(15,*) n1,n2,n3
            read(15,*) L3, L2, L1
            read(15,*) stretch_Y, mesh_type

            close(15)

            n1m=n1-1                !number of spanwise cells
            n2m=n2-1                !number of normal cells
            n3m=n3-1                !number of streamwise cells

            L1=L1*pi       !spanwise size of the whole calcul box
            L3=L3*pi       !streamwise size of the whole calcul box

        end subroutine read_domain_settings

        subroutine read_global_settings(nb_iteration, save_3D_frequency, save_recovery_frequency)

            use boundaries
            use DNS_settings
            use bubble_data
            use inflow_settings

            implicit none

            integer,intent(out)         :: nb_iteration, save_3D_frequency, save_recovery_frequency

            open(15,file=trim(COMMON_settings_path)//'global.d')

            read(15,*) nb_iteration, save_3D_frequency, save_recovery_frequency
            read(15,*) BC1, BC2, BC3
            read(15,*) streamwise
            read(15,*) flow_type
            read(15,*) outflow_buff
            read(15,*) inflow_buff
            read(15,*) inflow_mode
            read(15,'(a)')COMMON_inflow_path  ! path for inflow
            read(15,*) divro
            read(15,*) d ! with physical unity (mm)
            read(15,*) g ! with physical unity (mm/s/s)
            read(15,*) ren
            read(15,*) nu_viscosity !(mm*mm/s) characteristic velocity
            read(15,*) h_height !(mm) characteristic length
            read(15,*) Uc !(mm/s) characteristic velocity
            read(15,*) save_gradP_frequency ! mean_gradP export frequency


!            open(1544, file="global_parameters", position="append")
!
!            write(1544,*)'*********************************'
!            write(1544,*)'The Re number in the simulation is:  ',ren
!            write(1544,*)'*********************************'
            g_modified = g*h_height/(Uc**2) !modification of expression : h_height instead of L2
            d_modified = d/h_height
            rad_modified = d_modified*0.5d0

!            write(1544,*)'*********************************'
!            write(1544,*)'The gamma number in the simulation is:  ',g_modified
!            write(1544,*)'*********************************'
            Re_bubble = abs(g)*(d**3)/(12*nu_viscosity)

!            if (Re_bubble <= 0.1) then
!                write(1544,*)'The limit Reynolds of the Bubble is:  ',Re_bubble
!            else if (Re_bubble > 0.1) then
!
!                write(1544,*)'The Reynolds of the Bubble exeed 0.1:  ',Re_bubble
!                write(1544,*)'WARNING:  Change the Bubble Diameter'
!                write(*,*)'The Diameter is:  ', d
!     !*************OLD case cf Mauricio********
!           !     write(1544,*)'The Reynolds of the Bubble exeed 0.1:  ',Re_bubble
!           !     write(1544,*)'Changing the Bubble Diameter'
!           !      d = (1.2d0*(nu_viscosity)/abs(g))**(1/3)
!           !     write(*,*)'The new Diameter is:  ', d
!            end if
!            close(1544)

            close(15)


        end subroutine read_global_settings


!        subroutine read_MHD_settings()
!
!            use DNS_settings
!            use MHD_data
!
!            implicit none
!
!            open(15,file=trim(settings_path)//'MHD.d')
!
!             read(15,*) MHD_activated                          ! 0:OFF 1:ON
!             read(15,*) Hartmann_number                          !=> Leads also to the Hartmann number by combinng with the Reynolds number
!!            read(15,*) Monodirec_Magnetic_field               ! 0:OFF 1:ON
!             read(15,*) Magnetic_field_unit_vector_x
!             read(15,*) Magnetic_field_unit_vector_y
!             read(15,*) Magnetic_field_unit_vector_z
!             ! Warning the B0 'applied magnetic field' vector should be unit vector due to non dimensionalization
!             ! One direction only for the moment (eiher x , y or z)
!            close(15)
!
!          Stuart_number = Hartmann_number*Hartmann_number/ren
!
!         write(*,*)'Stuart_number =',  Stuart_number
!
!        end subroutine read_MHD_settings

        subroutine read_IBM_settings()

            use IBM_settings

            implicit none
            integer :: IBM_activated_int

            ! IBM **********************************************************************************
            open(15,file=trim(COMMON_settings_path)//'IBM.d')

            read(15,*) IBM_activated_int
            read(15,'(a)')obj_file_path  ! path of the read field  file (cha.rea)
            read(15,*) body_x1, body_x2, body_x3
            read(15,*) body_scale_x1, body_scale_x2, body_scale_x3
            close(15)

            IBM_activated=(IBM_activated_int==1)

        end subroutine read_IBM_settings


        subroutine read_animation_settings()

            use anim2D
            use animBubble

            implicit none
            integer                     :: n
            integer, dimension(1:7)     :: anim_flag
            character(200)              :: current_line

            open(15,file=trim(COMMON_settings_path)//'anim2D.d')

            read(15,*) anim2D_1, anim2D_2, anim2D_3

            do while (.true.)

                read(15, '(a)', end=10) current_line
                if ((anim2D_1/=0).and.(current_line(1:2)=="X1")) then
!                    write(*,*)'lol', current_line

                    read(15,*) param_anim2D_1%nb_slices, (param_anim2D_1%slices(n),n=1,param_anim2D_1%nb_slices)
                    read(15,*) param_anim2D_1%nb_steps, param_anim2D_1%step_size
                    read(15,*) anim_flag(1), anim_flag(2), anim_flag(3), anim_flag(4)
                    read(15,*) anim_flag(5)

                    param_anim2D_1%export_q1=(anim_flag(1)==1)
                    param_anim2D_1%export_q2=(anim_flag(2)==1)
                    param_anim2D_1%export_q3=(anim_flag(3)==1)
                    param_anim2D_1%export_pr=(anim_flag(4)==1)
                    param_anim2D_1%export_sca=(anim_flag(5)==1)
                endif

                if ((anim2D_2/=0).and.(current_line(1:2)=="X2")) then
!                    write(*,*)'lol', current_line

                    read(15,*) param_anim2D_2%nb_slices, (param_anim2D_2%slices(n),n=1,param_anim2D_2%nb_slices)
                    read(15,*) param_anim2D_2%nb_steps, param_anim2D_2%step_size
                    read(15,*) anim_flag(1), anim_flag(2), anim_flag(3), anim_flag(4)

                    param_anim2D_2%export_q1=(anim_flag(1)==1)
                    param_anim2D_2%export_q2=(anim_flag(2)==1)
                    param_anim2D_2%export_q3=(anim_flag(3)==1)
                    param_anim2D_2%export_pr=(anim_flag(4)==1)
                endif

                if ((anim2D_3/=0).and.(current_line(1:2)=="X3")) then
!                    write(*,*)'lol', current_line

                    read(15,*) param_anim2D_3%nb_slices, (param_anim2D_3%slices(n),n=1,param_anim2D_3%nb_slices)
                    read(15,*) param_anim2D_3%nb_steps, param_anim2D_3%step_size
                    read(15,*) anim_flag(1), anim_flag(2), anim_flag(3), anim_flag(4)
                    read(15,*) anim_flag(5)

                    param_anim2D_3%export_q1=(anim_flag(1)==1)
                    param_anim2D_3%export_q2=(anim_flag(2)==1)
                    param_anim2D_3%export_q3=(anim_flag(3)==1)
                    param_anim2D_3%export_pr=(anim_flag(4)==1)
                    param_anim2D_3%export_sca=(anim_flag(5)==1)
                endif

                if (current_line(1:2)=="Bu") then
!                    write(*,*)'lol', current_line

                    read(15,*) slice_cpt
                    read(15,*) (slice_tbl(n)%p1,n=1,slice_cpt)
                    read(15,*) (slice_tbl(n)%p2,n=1,slice_cpt)
                    read(15,*) (slice_tbl(n)%p3,n=1,slice_cpt)
                    read(15,*) (slice_tbl(n)%s1,n=1,slice_cpt)
                    read(15,*) (slice_tbl(n)%s2,n=1,slice_cpt)
                    read(15,*) (slice_tbl(n)%s3,n=1,slice_cpt)
                    read(15,*) (slice_tbl(n)%step_size, n=1,slice_cpt)

                endif

            end do

    10      close(15)
        end subroutine read_animation_settings


        ! MOVE IN BUBBLE MODULE !!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine read_bubble_settings()

            use bubble_data

            implicit none
            integer         :: n

            ! Bubbles ********************************************************************************
            open(15,file=trim(COMMON_settings_path)//'bubbles.d')

            read(15,*) sub_step
            read(15,*) bublgen_nb
            read(15,*) (bubl_gen(n)%n1,n=1,bublgen_nb)
            read(15,*) (bubl_gen(n)%n2,n=1,bublgen_nb)
            read(15,*) (bubl_gen(n)%n3,n=1,bublgen_nb)
            read(15,*) (bubl_gen(n)%s1,n=1,bublgen_nb)
            read(15,*) (bubl_gen(n)%s2,n=1,bublgen_nb)
            read(15,*) (bubl_gen(n)%s3,n=1,bublgen_nb)
            read(15,*) (bubl_gen(n)%nb_start,n=1,bublgen_nb)
            read(15,*) (bubl_gen(n)%nb_regen,n=1,bublgen_nb)
            close(15)
        end subroutine read_bubble_settings

        subroutine read_vortices_settings()

            use DNS_settings

            implicit none

            open(68,file=trim(COMMON_settings_path)//'vortices.d')

            read(68,*) vort_dir ! (Orientation)
            read(68,*) x1vs, hvort, x3vs !(starting point)
            read(68,*) Lvort !(vortex length)
            read(68,*) dvort !(vortex diameter)
            read(68,*) Urot !(rotation velocity)

            close(68)

        end subroutine read_vortices_settings


        subroutine resume_settings()

            use decomp_2d

            use DNS_settings
            use IBM_settings
            use start_settings
            use inflow_settings
            use boundaries
            use numerical_methods_settings
            use anim2D

            use schemes_loader

            implicit none
            integer :: s

            if (nrank==0) then

                write(*,*)'*****************************************************************'
                write(*,*)'*****************************************************************'
                write(*,*)'*                                                               *'
                write(*,*)'*                Fully 3D Turbulent Channel Flow                *'
                write(*,*)'*                         (Double precision)                    *'
                write(*,*)'*                      Parallal Version (v 5.0)                 *'
                write(*,*)'*                                                               *'
                write(*,*)'*****************************************************************'
                write(*,*)'*****************************************************************'
                write(*,*)'  '
                write(*,*) '  '



                write(*,*)'Simulation starting ---------------------------------------'

                select case (start_source_type)
                    case (NO_SOURCE)
                        write(*,*)'Begin from generated turbulent fields'

                    case (HDF5_FILE)
                        write(*,*)'Begin from HDF5 file'

                        if (start_from_coarse_file) then
                            write(*,*) '====> The field contains a reduced set of data'
                        end if


                    case default

                end select

                write(*,*)'Computational domain ---------------------------------------'
                write(*,*) '...... Resolution : n1=', n1, 'n2', n2, 'n3', n3
                write(*,*) '...... Size       : L1=', L1, 'L2', L2, 'L3', L3
                if (mesh_type==ORLANDI_MESH) write(*,*) '...... Mesh type  : ORLANDI (beta=', stretch_Y,')'
                if (mesh_type==LAMBALLAIS_MESH) write(*,*) '...... Mesh type  : LAMBALLAIS (beta=', stretch_Y,')'

                write(*,*)
                write(*,*)'2D domain decomposition ------------------------------------'
                write(*,*) '...... x(1)', xstart(1), 'to',xend(1)
                write(*,*) '...... x(2)', xstart(2), 'to',xend(2)
                write(*,*) '...... x(3)', xstart(3), 'to',xend(3)

                write(*,*)'ANIM2D -----------------------------------------------------'
                if (anim2D_1/=0) then

                    if (anim2D_1==1 )write(*,*) 'Anim x1 ? -YES , export in single .h5 file_________________'
                    if (anim2D_1==2 )write(*,*) 'Anim x1 ? -YES , export in multiple .h5 files______________'
                    write(*,*)
                    write(*,*) 'Export:'
                    if (param_anim2D_1%export_q1) write(*,*) '......Q1 : YES'
                    if (param_anim2D_1%export_q2) write(*,*) '......Q2 : YES'
                    if (param_anim2D_1%export_q3) write(*,*) '......Q3 : YES'
                    if (param_anim2D_1%export_pr) write(*,*) '......Pr : YES'
                    if (param_anim2D_1%export_sca) write(*,*) '......Scalar : YES'

                    write(*,*)
                    write(*,*)'At slices :',param_anim2D_1%slices(1:param_anim2D_1%nb_slices)

                    write(*,*)
                    write(*,*)'Save ', param_anim2D_1%nb_steps,'steps every :',param_anim2D_1%step_size,'iterations.'

                else
                    write(*,*) 'Anim x1 ? -NO'
                end if

                if (anim2D_2/=0) then

                    write(*,*) 'Anim x2 ? -YES _________________________'
                    write(*,*)
                    write(*,*) 'Export:'
                    if (param_anim2D_2%export_q1) write(*,*) '......Q1 : YES'
                    if (param_anim2D_2%export_q2) write(*,*) '......Q2 : YES'
                    if (param_anim2D_2%export_q3) write(*,*) '......Q3 : YES'
                    if (param_anim2D_2%export_pr) write(*,*) '......Pr : YES'

                    write(*,*)
                    write(*,*)'At slices :',param_anim2D_2%slices(1:param_anim2D_2%nb_slices)

                    write(*,*)
                    write(*,*)'Save ', param_anim2D_2%nb_steps,'steps every :',param_anim2D_2%step_size,'iterations.'
                else
                    write(*,*) 'Anim x2 ? -NO'
                end if

                if (anim2D_3/=0) then

                    if (anim2D_3==1 )write(*,*) 'Anim x3 ? -YES , export in single .h5 file_________________'
                    if (anim2D_3==2 )write(*,*) 'Anim x3 ? -YES , export in multiple .h5 files______________'
                    write(*,*)
                    write(*,*) 'Export:'
                    if (param_anim2D_3%export_q1) write(*,*) '......Q1 : YES'
                    if (param_anim2D_3%export_q2) write(*,*) '......Q2 : YES'
                    if (param_anim2D_3%export_q3) write(*,*) '......Q3 : YES'
                    if (param_anim2D_3%export_pr) write(*,*) '......Pr : YES'
                    if (param_anim2D_3%export_sca) write(*,*) '......Scalar : YES'

                    write(*,*)
                    write(*,*)'At slices :',param_anim2D_3%slices(1:param_anim2D_3%nb_slices)

                    write(*,*)
                    write(*,*)'Save ', param_anim2D_3%nb_steps,'steps every :',param_anim2D_3%step_size,'iterations.'

                else
                    write(*,*) 'Anim x3 ? -NO'
                end if

                write(*,*)'IBM --------------------------------------------------------'
                if (IBM_activated) then
                    write(*,*) 'IBM activated ? -YES'
                    write(*,*) '...... OBJECT file   :', trim(obj_file_path)
                    write(*,*) '...... Body position :', 'X1=', body_x1, "X2=",body_x2, "X3=", body_x3
                    write(*,*) '...... Scale         :', 'X1=', body_scale_x1, "X2=",body_scale_x2, "X3=", body_scale_x3
                else
                    write(*,*) 'IBM activated ? -NO'
                end if


                write(*,*)
                write(*,*)'Flow configuration -------------------------------------------'

                select case (BC1)

                    case (UNBOUNDED)
                        write(*,*)'BC1: Periodic'
                    case (FREESLIP)
                        write(*,*)'BC1: Freeslip'
                    case (NOSLIP)
                        write(*,*)'BC1: Noslip'
                    case (OPEN)
                        write(*,*)'BC1: In/Outflow'
                        if (flow_type==FLOW_FROM_INFLOW)write(*,*)'First field from outflow'
                        if (inflow_buff>0)write(*,*)'Inflow files path :', trim(COMMON_inflow_path)
                        if (inflow_buff>0)write(*,*)'Inflow buffer size:', inflow_buff
                    case default

                end select

                select case (BC2)

                    case (UNBOUNDED)
                        write(*,*)'BC2: Periodic'
                    case (FREESLIP)
                        write(*,*)'BC2: Freeslip'
                    case (NOSLIP)
                        write(*,*)'BC2: Noslip'
                    case (OPEN)
                        write(*,*)'BC2: In/Outflow'

                    case default

                end select

                select case (BC3)

                    case (UNBOUNDED)
                        write(*,*)'BC3: Periodic'
                    case (FREESLIP)
                        write(*,*)'BC3: Freeslip'
                    case (NOSLIP)
                        write(*,*)'BC3: Noslip'
                    case (OPEN)
                        write(*,*)'BC3: In/Outflow'

                    case default

                end select

                if (outflow_buff>0)write(*,*)'Outflow files path:', trim(COMMON_outflow_path)



                write(*,*)
                write(*,*)'Slows configuration ----------------------------------------'
                write(*,*)'Number of slows', nb_slows

                do s = 1, nb_slows
                    write(*,*)'Configuration of slow ', s, ':'
                    write(*,*)'    -->begin at : ', slows(s)%xst, ', ', slows(s)%zst
                    write(*,*)'    -->end at   : ', slows(s)%xen, ', ', slows(s)%zen
                    write(*,*)'    -->blowing  : ', slows(s)%blowing
                end do

                write(*,*)
                write(*,*)'Poisson equation resolution ----------------------------------'
                if (use_generic_poisson) then
                    write(*,*)'Use of spectral Poisson solver of Lamballais'
                else
                    write(*,*)'Use of XZ-spectral and Y-physical solver'
                end if

                write(*,*)
                write(*,*)'Schemes configuration ----------------------------------------'
                select case (schemes_configuration)

                    case (O2_SCHEMES)
                        write(*,*)'O2'

                    case (OPTIMIZED_SCHEMES)
                        write(*,*)'DRP'

                    case (OPTIMIZED_SCHEMES_WITH_CPT)
                        write(*,*)'Hybrid'

                    case (COMPACT_SCHEMES)
                        write(*,*)'Compact'

                    case default

                end select


            end if

        end subroutine resume_settings

end module CORE_IO_settings

module CORE_log_writers

    use physical_fields
    use SCALAR_data

    implicit none

contains

    subroutine write_sensorfile(i, j_v, k)

        use mpi
        use decomp_2d

        use COMMON_fieldstools
        use run_ctxt_data

        use IVELOCITY_workspace_view
        use VELOCITY_properties

        use SCALAR_data, only: SCA_state

        implicit none
        integer :: i, j_v, k
        real*8  :: u_s, v_s, w_s, cs1, cs2, cs3, cs4
        real*8  :: q1_cs, q2_cs, q3_cs, pr_cs
        real*8  :: q1w10_cs, q2w10_cs, q3w10_cs
        real*8  :: q1w20_cs, q2w20_cs, q3w20_cs
        real*8  :: q1w30_cs, q2w30_cs, q3w30_cs
        real*8  :: q1w11_cs, q2w11_cs, q3w11_cs
        real*8  :: q1w21_cs, q2w21_cs, q3w21_cs
        real*8  :: q1w31_cs, q2w31_cs, q3w31_cs
        real*8  :: q1_csg, q2_csg, q3_csg, pr_csg
        real*8  :: q1w10_csg, q2w10_csg, q3w10_csg
        real*8  :: q1w20_csg, q2w20_csg, q3w20_csg
        real*8  :: q1w30_csg, q2w30_csg, q3w30_csg
        real*8  :: q1w11_csg, q2w11_csg, q3w11_csg
        real*8  :: q1w21_csg, q2w21_csg, q3w21_csg
        real*8  :: q1w31_csg, q2w31_csg, q3w31_csg

        integer :: sensor_file_id=21, mpi_err

11      format(f18.6,',',x,f18.16,',',x,f18.16,',',x,f18.16,',',x,f18.16,',',5x,i6)

        call perform_checksum_x(q1_x(:,:,:), cs1)
        call perform_checksum_x(q2_x(:,:,:), cs2)
        call perform_checksum_x(q3_x(:,:,:), cs3)

        cs4=0.d0
        if(SCA_state/=0) call perform_checksum_x(sca_x(:,:,:), cs4)

        if ((zstart(2)<=j_v).and.(zend(2)>=j_v).and.(zstart(1)<=i).and.(zend(1)>=i).and.(zstart(3)<=k).and.(zend(3)>=k)) then

            open(sensor_file_id, file=VELOCITY_sensor_file, position="append")

            if (streamwise==1) write(sensor_file_id,*)t,cs3, cs2, cs1, cs4, ntime
            if (streamwise==3) write(sensor_file_id,*)t,cs1, cs2, cs3, cs4, ntime

            close(sensor_file_id)

            write(6,*)'_____________________________________________________________'
            write(6,*)
            write(6,*) 'Iteration n°', ntime
            write(6,*) '...... U:', cs3
            write(6,*) '...... V:', cs2
            write(6,*) '...... W:', cs1

            close(sensor_file_id)

        end if

        q1_cs=sum(abs(q1_x))
        q2_cs=sum(abs(q2_y))
        q3_cs=sum(abs(q3_y))
        pr_cs=sum(abs(dp_x))

        q1w10_cs=sum(abs(q1_wall10)); q2w10_cs=sum(abs(q2_wall10)); q3w10_cs=sum(abs(q3_wall10))
        q1w11_cs=sum(abs(q1_wall11)); q2w11_cs=sum(abs(q2_wall11)); q3w11_cs=sum(abs(q3_wall11))

        q1w20_cs=sum(abs(q1_wall20)); q2w20_cs=sum(abs(q2_wall20)); q3w20_cs=sum(abs(q3_wall20))
        q1w21_cs=sum(abs(q1_wall21)); q2w21_cs=sum(abs(q2_wall21)); q3w21_cs=sum(abs(q3_wall21))

        q1w30_cs=sum(abs(q1_wall30)); q2w30_cs=sum(abs(q2_wall30)); q3w30_cs=sum(abs(q3_wall30))
        q1w31_cs=sum(abs(q1_wall31)); q2w31_cs=sum(abs(q2_wall31)); q3w31_cs=sum(abs(q3_wall31))


        call MPI_ALLREDUCE (q1_cs, q1_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (q2_cs, q2_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (q3_cs, q3_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (pr_cs, pr_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        call MPI_ALLREDUCE (q1w10_cs, q1w10_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (q2w10_cs, q2w10_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (q3w10_cs, q3w10_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        call MPI_ALLREDUCE (q1w11_cs, q1w11_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (q2w11_cs, q2w11_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (q3w11_cs, q3w11_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        call MPI_ALLREDUCE (q1w20_cs, q1w20_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (q2w20_cs, q2w20_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (q3w20_cs, q3w20_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        call MPI_ALLREDUCE (q1w21_cs, q1w21_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (q2w21_cs, q2w21_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (q3w21_cs, q3w21_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        call MPI_ALLREDUCE (q1w20_cs, q1w20_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (q2w20_cs, q2w20_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (q3w20_cs, q3w20_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        call MPI_ALLREDUCE (q1w21_cs, q1w21_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (q2w21_cs, q2w21_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (q3w21_cs, q3w21_csg, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        if (nrank==-1) then

            write(*,*)
            write(*,*)'q1', q1_csg
            write(*,*)'q2', q2_csg
            write(*,*)'q3', q3_csg
            write(*,*)'pr', pr_csg


            write(*,*)'qi_wall10', q1w10_csg, q2w10_csg, q3w10_csg
            write(*,*)'qi_wall11', q1w11_csg, q2w11_csg, q3w11_csg

            write(*,*)'qi_wall20', q1w20_csg, q2w20_csg, q3w20_csg
            write(*,*)'qi_wall21', q1w21_csg, q2w21_csg, q3w21_csg

            write(*,*)'qi_wall30', q1w30_csg, q2w30_csg, q3w30_csg
            write(*,*)'qi_wall31', q1w31_csg, q2w31_csg, q3w31_csg

        endif

    end subroutine write_sensorfile

end module CORE_log_writers

module animation

    use decomp_2d

    implicit none

    integer, parameter          :: MAXNAME=20, MAXPATH=200, MAXFIELDS=20

    type anim2D
        character(MAXPATH)                                  :: path
        character(MAXNAME)                                  :: name
        integer, dimension(2)                               :: sz
        integer                                             :: max_nbframes
        integer                                             :: nbframes

        character(MAXNAME), dimension(MAXFIELDS)            :: global_fields
        integer                                             :: nb_global_fields
        character(MAXNAME), dimension(MAXFIELDS)            :: frames_fields
        integer                                             :: nb_frame_fields

    end type anim2D

    type(anim2D),dimension(20)  :: animations
    integer                     :: anim_ptr=0


contains

    subroutine open(path, name, nbframes, sz, id)
        implicit none
        character(*)          :: path
        character(*)           :: name
        integer                 :: nbframes
        integer, dimension(2)   :: sz
        integer                 :: id

        anim_ptr=anim_ptr+1
        id=anim_ptr

        animations(id)%path=path
        animations(id)%name=trim(name)
        animations(id)%max_nbframes=nbframes
        animations(id)%sz=sz

    end subroutine open

    subroutine declare_globalfields(id, global_fields, N)
        implicit none
        integer                         :: N, id
        character(20), dimension(:)     :: global_fields

        integer                         :: i

        do i = 1, N
            animations(id)%global_fields(i)=global_fields(i)
        end do
        animations(id)%nb_global_fields=N

    end subroutine declare_globalfields

    subroutine declare_framefields(id, frames_fields, N)
        implicit none
        integer                         :: N, id
        character(20), dimension(:)     :: frames_fields

        integer                         :: i

        do i = 1, N
            animations(id)%frames_fields(i)=frames_fields(i)
        end do
        animations(id)%nb_frame_fields=N

    end subroutine declare_framefields

    subroutine add_frame(id)
        implicit none
        integer             :: id

        if (animations(id)%nbframes < animations(id)%max_nbframes) then
            animations(id)%nbframes=animations(id)%nbframes+1
        else
            animations(id)%nbframes=0
            call export(id)
        endif



    end subroutine add_frame

    subroutine see(id)
        implicit none
        integer                 :: id, i

        if (nrank==0) then

            write(*,*)'ANIMATION numero', id
            write(*,*)'PATH:      ', trim(animations(id)%path)
            write(*,*)'NAME:      ', trim(animations(id)%name)
            write(*,*)'MAXFRAME:  ', animations(id)%max_nbframes
            write(*,*)'DIMENSION: ', animations(id)%sz

            write(*,*)'NBFRAME:            ', animations(id)%nbframes

            write(*,*)'GLOBAL FIELDS:      '
            do i = 1, animations(id)%nb_global_fields
                write(*,*)trim(animations(id)%global_fields(i))
            end do
            write(*,*)'NB GLOBAL FIELDS:   ', animations(id)%nb_global_fields

            write(*,*)'FRAME FIELDS:       '
            do i = 1, animations(id)%nb_frame_fields
                write(*,*)trim(animations(id)%frames_fields(i))
            end do

            write(*,*)'NB FRAME FIELDS:   ', animations(id)%nb_frame_fields

        endif

    end subroutine see


    subroutine export(id)
        implicit none
        integer                 :: id

        write(*,*)'EXPORT ANIMATION', id

    end subroutine export

end module animation

module anim2D_writer

    use mpi
    use decomp_2d
    use bubble_data
    use physical_fields

    implicit none

contains

    ! Export bubbles or particle located in subdomains in VTK file (can be viewed with paraview)
    ! The subdomain and export frequency for each subdomain is given by slice_tbl

    ! AFAIRE ajouter l'export des bulles dans une tranche
    subroutine anim2D_bubble(nt) ! Exports xb & vb as a vtk file, (format used by Paraview)

        use DNS_settings
        use mesh
        use subdomains_view
        use bubble_data
        use bubble_parallel_data
        use animBubble
        use VTK

        implicit none

        integer, intent(in)                         :: nt

        real*8                                      :: x1s, x2s, x3s, x1e, x2e, x3e
        integer                                     :: i, target

        target=0

        ! For each subdomain
        do i = 1, slice_cpt
            ! Export at the frequency defined by slice_tbl%step_size
            if (mod(nt, slice_tbl(i)%step_size)==0) then

                x1s=Z(max(1, (slice_tbl(i)%p1-slice_tbl(i)%s1)))
                x2s=Y(max(1, (slice_tbl(i)%p2-slice_tbl(i)%s2)))
                x3s=X(max(1, (slice_tbl(i)%p3-slice_tbl(i)%s3)))

                x1e=Z(min((slice_tbl(i)%p1+slice_tbl(i)%s1), n1))
                x2e=Y(min((slice_tbl(i)%p2+slice_tbl(i)%s2), n2))
                x3e=X(min((slice_tbl(i)%p3+slice_tbl(i)%s3), n3))

                call export_bubble(1, bubl, MAX_BUBBLES, bubble_cpt, x1s, x2s, x3s, x1e, x2e, x3e, i)
                call export_bubble(2, bubl_tmp, MAX_BUBBLES_TMP, bubble_tmp_cpt, x1s, x2s, x3s, x1e, x2e, x3e, i)

            endif

        end do


    contains

        subroutine export_bubble(file_nb, import, imp_size, imp_cpt, x1s, x2s, x3s, x1e, x2e, x3e, slice)

            use COMMON_workspace_view

            implicit none

            real*8, intent(in)                                              :: x1s, x2s, x3s, x1e, x2e, x3e
            integer, intent(in)                                             :: file_nb, imp_size, imp_cpt, slice
            type(bubble), dimension(imp_size), intent(in)                   :: import

            real*8, dimension(:,:), allocatable, save                       :: positions, vitesses, term1, term2, term3
            real*8, dimension(:), allocatable, save                         :: rank, term4
            type(bubble), dimension(imp_size)                               :: bubbles_to_export
            character(len=100)                                              :: vtk_out
            type(VTK_file_handle), dimension(MAX_EXPORT, MAX_SLICES), save  :: fd
            integer                                                         :: i, n
            character*10                                                    :: tmp_str, tmp_str2

            if (nrank/=target) then
                do i = 1, imp_cpt
                    if ((import(i)%x3<=x3e).and.(import(i)%x3>=x3s).and.(import(i)%x2<=x2e).and.(import(i)%x2>=x2s).and.(import(i)%x1<=x1e).and.(import(i)%x1>=x1s)) call post_bubble(target, import(i))
                end do
            end if

            n=0
            if (nrank==target) then
                do i = 1, imp_cpt
                    if ((import(i)%x3<=x3e).and.(import(i)%x3>=x3s).and.(import(i)%x2<=x2e).and.(import(i)%x2>=x2s).and.(import(i)%x1<=x1e).and.(import(i)%x1>=x1s)) call add_particule(bubbles_to_export, n, import(i)%x1,import(i)%x2,import(i)%x3, import(i)%v1,import(i)%v2,import(i)%v3, import(i)%u01,import(i)%u02,import(i)%u03, import(i)%u11,import(i)%u12,import(i)%u13, import(i)%om01, import(i)%om02, import(i)%om03, import(i)%om11, import(i)%om12, import(i)%om13, import(i)%a1, import(i)%a2, import(i)%a3, import(i)%b1, import(i)%b2, import(i)%b3, import(i)%c1, import(i)%c2, import(i)%c3, import(i)%d, import(i)%proc)
                end do
            end if

            call send_bubble(bubbles_to_export, n)

            if (allocated(positions)) deallocate(positions)
            if (allocated(vitesses)) deallocate(vitesses)
            if (allocated(rank)) deallocate(rank)
            if (allocated(term1)) deallocate(term1)
            if (allocated(term2)) deallocate(term2)
            if (allocated(term3)) deallocate(term3)
            if (allocated(term4)) deallocate(term4)

            if (nrank==target) then

                allocate(positions(3, n))
                allocate(vitesses(3, n))
                allocate(term1(3, n))
                allocate(term2(3, n))
                allocate(term3(3, n))
                allocate(rank(n))
                allocate(term4(n))

                positions(1,:)=bubbles_to_export(1:n)%x1
                positions(2,:)=bubbles_to_export(1:n)%x2
                positions(3,:)=bubbles_to_export(1:n)%x3

                vitesses(1,:)=bubbles_to_export(1:n)%v1
                vitesses(2,:)=bubbles_to_export(1:n)%v2
                vitesses(3,:)=bubbles_to_export(1:n)%v3

                term1(1,:)=bubbles_to_export(1:n)%a1
                term1(2,:)=bubbles_to_export(1:n)%a2
                term1(3,:)=bubbles_to_export(1:n)%a3

                term2(1,:)=bubbles_to_export(1:n)%b1
                term2(2,:)=bubbles_to_export(1:n)%b2
                term2(3,:)=bubbles_to_export(1:n)%b3

                term3(1,:)=bubbles_to_export(1:n)%c1
                term3(2,:)=bubbles_to_export(1:n)%c2
                term3(3,:)=bubbles_to_export(1:n)%c3

                term4=bubbles_to_export(1:n)%d

                rank=bubbles_to_export(1:n)%proc

                !                write(*,*) 'POSITION1:', maxval(positions(1,:))
                !                write(*,*) 'POSITION2:', maxval(positions(2,:))
                !                write(*,*) 'POSITION3:', maxval(positions(3,:))
                !
                !                write(*,*) 'VITESSE1:', maxval(vitesses(1,:))
                !                write(*,*) 'VITESSE2:', maxval(vitesses(2,:))
                !                write(*,*) 'VITESSE3:', maxval(vitesses(3,:))

                write(tmp_str, "(i10)")slice
                write(tmp_str2, "(i10)")file_nb
                vtk_out=trim(adjustl(COMMON_anim2D_bubble_path))//'bubbles_out_'//trim(adjustl(tmp_str2))//'_slice_'//trim(adjustl(tmp_str))

                call VTK_open_file(PREFIX=vtk_out, FD=fd(file_nb, slice))
                call VTK_write_mesh(FD=fd(file_nb, slice), XYZ=positions(:,:), n=n)
                call VTK_write_var(FD=fd(file_nb, slice),name="Vit_long",FIELD=rank)
                call VTK_write_var(FD=fd(file_nb, slice),name="v1",FIELD=vitesses(1,:))
                call VTK_write_var(FD=fd(file_nb, slice),name="Vitesse",VX=vitesses(1,:),VY=vitesses(2,:),VZ=vitesses(3,:), n=n)
                call VTK_write_var(FD=fd(file_nb, slice),name="Vorticite",VX=term1(1,:),VY=term1(2,:),VZ=term1(3,:), n=n)
                call VTK_write_var(FD=fd(file_nb, slice),name="Drag_force",VX=term2(1,:),VY=term2(2,:),VZ=term2(3,:), n=n)
                call VTK_write_var(FD=fd(file_nb, slice),name="Term3",VX=term3(1,:),VY=term3(2,:),VZ=term3(3,:), n=n)
                call VTK_write_var(FD=fd(file_nb, slice),name="Arkhimedes",FIELD=term4)

                call VTK_close_file(FD=fd(file_nb, slice))

            end if

        end subroutine export_bubble

    end subroutine anim2D_bubble


    subroutine anim2D_addframe1(anim_dir)
        use anim2D
        use HDF5_IO
        use IBM_data
        use physical_fields
        use mesh
        use mpi
        use IBM_settings
        use SCALAR_data

        implicit none
        character(*)            :: anim_dir
        integer, save           :: file_nb=1
        integer, save           :: nbsteps=1
        character(200), save    :: file_path
        character(20), save     :: file_name
!        character(200)          :: file_path
!        character(20)           :: file_name
        character(20)           :: frame_name, slice_name
        character(10)           :: tmp_str

        integer         :: i,s, xdmf_id, ierr

        xdmf_id=511

        !write(*,*)'subroutine anim2D_addframe1'

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! ACTIVATE this part if you want to export animation to various .h5 files !!!!
        !!!! AVOIDS problems while computing on various NODES !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! USE the appropriate write_xdmf2 subroutine !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (anim2D_1==2) then
            write(tmp_str, "(i10)")nbsteps
            file_name='Frame'//trim(adjustl(tmp_str))
            frame_name=''
            file_path=trim(anim_dir)//'X1/'//trim(file_name)
            if(nrank==0)  call hdf_create_file_with_2Dmesh(trim(file_path), Y, X, "x2", "x3", n2, n3)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if ((file_nb==1).and.(nbsteps==1)) then

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! ACTIVATE this part if you want to export animation to a single .h5 file !!!!
        !!!! You also need to activate the creation of the Frame group !!!!!!!!!!!!!!!!!!
        !!!! MIGHT cause problems while computing on various NODES !!!!!!!!!!!!!!!!!!!!!!
        !!!! USE the appropriate write_xdmf2 subroutine !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (anim2D_1==1) then
                file_name="anim1"
                file_path=trim(anim_dir)//'X1/'//trim(file_name)
                if(nrank==0)  call hdf_create_file_with_2Dmesh(trim(file_path), Y, X, "x2", "x3", n2, n3)
            endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (IBM_activated) then

                if(nrank==0)  call hdf_addgroup(file_path, "MASKS")

                do s = 1, param_anim2D_1%nb_slices

                    i=param_anim2D_1%slices(s)
                    write(tmp_str, "(i10)")i
                    slice_name='MASKS/'//'slice_'//trim(adjustl(tmp_str))
                    if (nrank==0) write(*,*)'ADD GROUPE 1', slice_name
                    if(nrank==0)  call hdf_addgroup(file_path, slice_name)

                    if (param_anim2D_1%export_q1) then
                        call hdf_add_2Dfield(trim(file_path), IBM_mask1(i, :,:), trim(slice_name)//"/q1", n2, n3, xstart(2),xend(2), xstart(3),xend(3))
                    end if

                    if (param_anim2D_1%export_q2) then
                        call hdf_add_2Dfield(trim(file_path), IBM_mask2(i, :,:), trim(slice_name)//"/q2", n2, n3, xstart(2),xend(2), xstart(3),xend(3))
                    end if

                    if (param_anim2D_1%export_q3) then
                        call hdf_add_2Dfield(trim(file_path), IBM_mask2(i, :,:), trim(slice_name)//"/q3", n2, n3, xstart(2),xend(2), xstart(3),xend(3))
                    end if

                end do

            endif

        end if


        ! ADDING THE CURRENT STEP
        if (nrank==0) write(*,*) 'ADDING FRAME n°', nbsteps, 'in file : ', trim(file_path)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! ACTIVATE this part if you want to export animation to a single .h5 file !!!!
        !!!! MIGHT cause problems while computing on various NODES !!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (anim2D_1==1) then
            write(tmp_str, "(i10)")nbsteps
            frame_name='Frame'//trim(adjustl(tmp_str))//'/'
            if(nrank==0) write(*,*)'ADD_GROUPE : ', frame_name
            if(nrank==0)  call hdf_addgroup(file_path, frame_name)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        do s = 1, param_anim2D_1%nb_slices
            i=param_anim2D_1%slices(s)
            write(tmp_str, "(i10)")i
            slice_name=trim(frame_name)//'slice_'//trim(adjustl(tmp_str))
            !if (nrank==0) write(*,*)'ADD_GROUPE', slice_name
            if(nrank==0)  call hdf_addgroup(file_path, slice_name)

            if (param_anim2D_1%export_q1) then
                call hdf_add_2Dfield(trim(file_path), q1_x(i, :,:), trim(slice_name)//"/q1", n2, n3, xstart(2),xend(2), xstart(3),xend(3))
            end if

            if (param_anim2D_1%export_q2) then
                call hdf_add_2Dfield(trim(file_path), q2_x(i, :,:), trim(slice_name)//"/q2", n2, n3, xstart(2),xend(2), xstart(3),xend(3))
            end if

            if (param_anim2D_1%export_q3) then
                call hdf_add_2Dfield(trim(file_path), q3_x(i, :,:), trim(slice_name)//"/q3", n2, n3, xstart(2),xend(2), xstart(3),xend(3))
            end if

            if (param_anim2D_1%export_pr) then
                call hdf_add_2Dfield(trim(file_path), dp_x(i, :,:), trim(slice_name)//"/pr", n2, n3, xstart(2),xend(2), xstart(3),xend(3))
            end if

            if ((param_anim2D_1%export_sca).and.(SCA_state/=0)) then
                call hdf_add_2Dfield(trim(file_path), sca_x(i, :,:), trim(slice_name)//"/sca", n2, n3, xstart(2),xend(2), xstart(3),xend(3))
            end if


        end do

        if (anim2D_1==1) call write_xdmf2_single_HDF5
        if (anim2D_1==2) call write_xdmf2_multi_HDF5

        nbsteps=nbsteps+1

!        if (nbsteps>param_anim2D_1%nb_steps) then
!
!            nbsteps=1
!            file_nb=file_nb+1
!
!            write(tmp_str, "(i10)")file_nb
!            file_name="anim"//trim(adjustl(tmp_str))
!            file_path=trim(anim_dir)//'X1/'//trim(file_name)
!            if(nrank==0)  call hdf_create_file_with_2Dmesh(trim(file_path), Y, X, "x2", "x3", n2, n3)
!
!
!            if (IBM_activated) then
!
!                if(nrank==0)  call hdf_addgroup(file_path, "MASKS")
!
!                do s = 1, param_anim2D_1%nb_slices
!
!                    i=param_anim2D_1%slices(s)
!                    write(tmp_str, "(i10)")i
!                    slice_name='MASKS/'//'slice_'//trim(adjustl(tmp_str))
!                    if (nrank==0) write(*,*)'ADD GROUPE 1', slice_name
!                    if(nrank==0)  call hdf_addgroup(file_path, slice_name)
!
!                    if (param_anim2D_1%export_q1) then
!                        call hdf_add_2Dfield(trim(file_path), IBM_mask1(i, :,:), trim(slice_name)//"/q1", n2, n3, xstart(2),xend(2), xstart(3),xend(3))
!                    end if
!
!                    if (param_anim2D_1%export_q2) then
!                        call hdf_add_2Dfield(trim(file_path), IBM_mask2(i, :,:), trim(slice_name)//"/q2", n2, n3, xstart(2),xend(2), xstart(3),xend(3))
!                    end if
!
!                    if (param_anim2D_1%export_q3) then
!                        call hdf_add_2Dfield(trim(file_path), IBM_mask2(i, :,:), trim(slice_name)//"/q3", n2, n3, xstart(2),xend(2), xstart(3),xend(3))
!                    end if
!
!                end do
!
!            end if
!
!        end if

    contains

        subroutine write_xdmf2_single_HDF5()
            use file_copy
            use DNS_settings, only:dt
            use run_ctxt_data, only:tinit
            implicit none
            integer     :: f
            character(200)  :: xdmf_path, attribute_path
            character(20)   :: frame_label, attribute_name
            character(80)   :: tmp_str, tmp_str2, tmp_str3, tmp_str4

            write(tmp_str, "(i10)")file_nb
            xdmf_path=trim(anim_dir)//'X1/anim'//trim(adjustl(tmp_str))//".xdmf"

            if(nrank==0)  call copy_ascii_file(file_dest=trim(xdmf_path), file_src=trim(anim_dir)//"X1/modele2.xdmf")

            if(nrank==0)  open(xdmf_id, file=trim(xdmf_path), position="append")
            if(nrank==0)  write(xdmf_id,'(F10.5,F10.5,I4)')tinit, dt*param_anim2D_1%step_size, param_anim2D_1%nb_steps
            if(nrank==0)  write(xdmf_id,*)"</DataItem>"
            if(nrank==0)  write(xdmf_id,*)"</Time>"
            if(nrank==0)  write(xdmf_id,*)


            do f = 1, nbsteps

                write(tmp_str, "(i10)")f
                frame_label='Frame'//trim(adjustl(tmp_str))
                if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Grid Name="', trim(frame_label), '" GridType="Uniform">'
                if(nrank==0)  write(xdmf_id,*)

                if(nrank==0)  write(xdmf_id,'(a,I4,I4,a)')'<Topology TopologyType="2DRectMesh" NumberOfElements="', n3, n2,'"/>'
                if(nrank==0)  write(xdmf_id,*)'<Geometry GeometryType="VXVY">'

                if(nrank==0)  write(xdmf_id,'(a,I4,a)')'<DataItem Dimensions="',n3,'" NumberType="Float" Precision="8" Format="HDF">'
                if(nrank==0)  write(xdmf_id,*) trim(file_name)//".h5:/x3"
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"

                if(nrank==0)  write(xdmf_id,'(a,I4,a)')'<DataItem Dimensions="',n2,'" NumberType="Float" Precision="8" Format="HDF">'
                if(nrank==0)  write(xdmf_id,*)trim(file_name)//".h5:/x2"
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"


                if(nrank==0)  write(xdmf_id,*)"</Geometry>"

                if(nrank==0)  write(xdmf_id,*)

                do s = 1, param_anim2D_1%nb_slices

                    i=param_anim2D_1%slices(s)
                    write(tmp_str, "(i10)")i
                    slice_name=trim(frame_label)//'/slice_'//trim(adjustl(tmp_str))

                    if(nrank==0)  write(xdmf_id,*)

                    if (param_anim2D_1%export_q1) then
                        attribute_name='q1_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q1'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
                    end if

                    if (param_anim2D_1%export_q2) then
                        attribute_name='q2_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q2'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"

                    end if

                    if (param_anim2D_1%export_q3) then
                        attribute_name='q3_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q3'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
                    end if

                    if (param_anim2D_1%export_pr) then
                        attribute_name='pr_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/pr'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
                    end if

                    if ((param_anim2D_1%export_sca).and.(SCA_state/=0)) then
                        attribute_name='sca_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/sca'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
                    end if

!                    if (.true.) then
!                        attribute_name='fb1_scalar_'//trim(adjustl(tmp_str))
!                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/fb_sca_2'
!
!                        tmp_str2='<DataItem Dimensions="'
!                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
!                        tmp_str4='</DataItem>'
!                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
!                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
!                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
!                    end if

                end do


                if(nrank==0)  write(xdmf_id,*)
                if(nrank==0)  write(xdmf_id,*)
                if(nrank==0)  write(xdmf_id,*)"</Grid>"
                if(nrank==0)  write(xdmf_id,*)
            end do


            if(nrank==0)  write(xdmf_id,*)"</Grid>"
            if(nrank==0)  write(xdmf_id,*)


            ! WRITTING MASKS SECTION


            if (IBM_activated) then

                frame_label='MASKS'
                if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Grid Name="', trim(frame_label), '" GridType="Uniform">'
                if(nrank==0)  write(xdmf_id,*)

                if(nrank==0)  write(xdmf_id,'(a,I4,I4,a)')'<Topology TopologyType="2DRectMesh" NumberOfElements="', n3, n2,'"/>'
                if(nrank==0)  write(xdmf_id,*)'<Geometry GeometryType="VXVY">'

                if(nrank==0)  write(xdmf_id,'(a,I4,a)')'<DataItem Dimensions="',n3,'" NumberType="Float" Precision="8" Format="HDF">'
                if(nrank==0)  write(xdmf_id,*) trim(file_name)//".h5:/x3"
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"

                if(nrank==0)  write(xdmf_id,'(a,I4,a)')'<DataItem Dimensions="',n2,'" NumberType="Float" Precision="8" Format="HDF">'
                if(nrank==0)  write(xdmf_id,*)trim(file_name)//".h5:/x2"
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"


                if(nrank==0)  write(xdmf_id,*)"</Geometry>"

                if(nrank==0)  write(xdmf_id,*)

                do s = 1, param_anim2D_1%nb_slices

                    i=param_anim2D_1%slices(s)
                    write(tmp_str, "(i10)")i
                    slice_name=trim(frame_label)//'/slice_'//trim(adjustl(tmp_str))

                    if(nrank==0)  write(xdmf_id,*)

                    if (param_anim2D_1%export_q1) then
                        attribute_name='q1_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q1'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
                    end if

                    if (param_anim2D_1%export_q2) then
                        attribute_name='q2_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q2'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"

                    end if

                    if (param_anim2D_1%export_q3) then
                        attribute_name='q3_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q3'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
                    end if

                end do


                if(nrank==0)  write(xdmf_id,*)
                if(nrank==0)  write(xdmf_id,*)
                if(nrank==0)  write(xdmf_id,*)"</Grid>"
                if(nrank==0)  write(xdmf_id,*)

            endif


            if(nrank==0)  write(xdmf_id,*)"</Domain>"
            if(nrank==0)  write(xdmf_id,*)"</Xdmf>"
            if(nrank==0)  close(xdmf_id)

        end subroutine write_xdmf2_single_HDF5

        subroutine write_xdmf2_multi_HDF5()
            use file_copy
            use DNS_settings, only:dt
            use run_ctxt_data, only:tinit
            implicit none
!            integer     :: f
            character(200)  :: xdmf_path, attribute_path
            character(20)   :: frame_label, attribute_name
            character(80)   :: tmp_str, tmp_str2, tmp_str3, tmp_str4
            integer,save    :: nline_file

            write(tmp_str, "(i10)")file_nb
            xdmf_path=trim(anim_dir)//'X1/anim'//trim(adjustl(tmp_str))//".xdmf"

            if ((nbsteps==1).and.(nrank==0)) call copy_ascii_file(file_dest=trim(xdmf_path), &
                                                                  file_src=trim(anim_dir)//"X1/modele2.xdmf",nline=nline_file)

            if(nrank==0)  open(xdmf_id, file=trim(xdmf_path), position="rewind")

            if (nbsteps==1) then
                !!!! SKIPPING the already existing text !!!!
                do i=1,nline_file
                    if(nrank==0) read(xdmf_id,*)
                enddo

                if(nrank==0)  write(xdmf_id,'(F10.5,F10.5,I4)')tinit, dt*param_anim2D_1%step_size, param_anim2D_1%nb_steps
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"
                if(nrank==0)  write(xdmf_id,*)"</Time>"
                if(nrank==0)  write(xdmf_id,*)

                nline_file=nline_file+4
            else
                !!!! SKIPPING the already existing text, excepting the last 2 lines (Domain and Xdmf) !!!!
                do i=1,nline_file
                    if(nrank==0) read(xdmf_id,*)
                enddo
            endif


!            do f = 1, nbsteps

                write(tmp_str, "(i10)")nbsteps
                frame_label='Frame'//trim(adjustl(tmp_str))
                if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Grid Name="', trim(frame_label), '" GridType="Uniform">'
                if(nrank==0)  write(xdmf_id,*)

                if(nrank==0)  write(xdmf_id,'(a,I4,I4,a)')'<Topology TopologyType="2DRectMesh" NumberOfElements="', n3, n2,'"/>'
                if(nrank==0)  write(xdmf_id,*)'<Geometry GeometryType="VXVY">'

                if(nrank==0)  write(xdmf_id,'(a,I4,a)')'<DataItem Dimensions="',n3,'" NumberType="Float" Precision="8" Format="HDF">'
                if(nrank==0)  write(xdmf_id,*) trim(file_name)//".h5:/x3"
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"

                if(nrank==0)  write(xdmf_id,'(a,I4,a)')'<DataItem Dimensions="',n2,'" NumberType="Float" Precision="8" Format="HDF">'
                if(nrank==0)  write(xdmf_id,*)trim(file_name)//".h5:/x2"
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"


                if(nrank==0)  write(xdmf_id,*)"</Geometry>"

                if(nrank==0)  write(xdmf_id,*)

                nline_file=nline_file+12

                do s = 1, param_anim2D_1%nb_slices

                    i=param_anim2D_1%slices(s)
                    write(tmp_str, "(i10)")i
                    slice_name='slice_'//trim(adjustl(tmp_str))

                    if(nrank==0)  write(xdmf_id,*)

                    nline_file=nline_file+1

                    if (param_anim2D_1%export_q1) then
                        attribute_name='q1_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q1'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"

                        nline_file=nline_file+3
                    end if

                    if (param_anim2D_1%export_q2) then
                        attribute_name='q2_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q2'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"

                        nline_file=nline_file+3

                    end if

                    if (param_anim2D_1%export_q3) then
                        attribute_name='q3_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q3'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"

                        nline_file=nline_file+3
                    end if

                    if (param_anim2D_1%export_pr) then
                        attribute_name='pr_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/pr'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"

                        nline_file=nline_file+3
                    end if

                    if ((param_anim2D_1%export_sca).and.(SCA_state/=0)) then
                        attribute_name='sca_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/sca'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"

                        nline_file=nline_file+3
                    end if


                end do


                if(nrank==0)  write(xdmf_id,*)
                if(nrank==0)  write(xdmf_id,*)
                if(nrank==0)  write(xdmf_id,*)"</Grid>"
                if(nrank==0)  write(xdmf_id,*)

                nline_file=nline_file+4
!            end do


            if(nrank==0)  write(xdmf_id,*)"</Grid>"
            if(nrank==0)  write(xdmf_id,*)


            ! WRITTING MASKS SECTION


            if (IBM_activated) then

                frame_label='MASKS'
                if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Grid Name="', trim(frame_label), '" GridType="Uniform">'
                if(nrank==0)  write(xdmf_id,*)

                if(nrank==0)  write(xdmf_id,'(a,I4,I4,a)')'<Topology TopologyType="2DRectMesh" NumberOfElements="', n3, n2,'"/>'
                if(nrank==0)  write(xdmf_id,*)'<Geometry GeometryType="VXVY">'

                if(nrank==0)  write(xdmf_id,'(a,I4,a)')'<DataItem Dimensions="',n3,'" NumberType="Float" Precision="8" Format="HDF">'
                if(nrank==0)  write(xdmf_id,*) trim(file_name)//".h5:/x3"
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"

                if(nrank==0)  write(xdmf_id,'(a,I4,a)')'<DataItem Dimensions="',n2,'" NumberType="Float" Precision="8" Format="HDF">'
                if(nrank==0)  write(xdmf_id,*)trim(file_name)//".h5:/x2"
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"


                if(nrank==0)  write(xdmf_id,*)"</Geometry>"

                if(nrank==0)  write(xdmf_id,*)

                do s = 1, param_anim2D_1%nb_slices

                    i=param_anim2D_1%slices(s)
                    write(tmp_str, "(i10)")i
                    slice_name=trim(frame_label)//'/slice_'//trim(adjustl(tmp_str))

                    if(nrank==0)  write(xdmf_id,*)

                    if (param_anim2D_1%export_q1) then
                        attribute_name='q1_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q1'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
                    end if

                    if (param_anim2D_1%export_q2) then
                        attribute_name='q2_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q2'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"

                    end if

                    if (param_anim2D_1%export_q3) then
                        attribute_name='q3_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q3'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n3,n2, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
                    end if

                end do


                if(nrank==0)  write(xdmf_id,*)
                if(nrank==0)  write(xdmf_id,*)
                if(nrank==0)  write(xdmf_id,*)"</Grid>"
                if(nrank==0)  write(xdmf_id,*)

            endif


            if(nrank==0)  write(xdmf_id,*)"</Domain>"
            if(nrank==0)  write(xdmf_id,*)"</Xdmf>"
            if(nrank==0)  close(xdmf_id)

        end subroutine write_xdmf2_multi_HDF5

    end subroutine anim2D_addframe1


    subroutine anim2D_addframe3(anim_dir)
        use animation
        use anim2D
        use HDF5_IO
        use IBM_data
        use physical_fields
        use mesh
        use mpi
        use IBM_settings
        use SCALAR_data
        use mesh_generator
        use DNS_settings

        implicit none
        character(*)            :: anim_dir
        integer, save           :: file_nb=1
        integer, save           :: nbsteps=1
        character(200), save    :: file_path
        character(20), save     :: file_name
        character(20)           :: frame_name, slice_name
        character(10)           :: tmp_str
        integer                 :: i,j

        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))  :: IBM_mask1_z, IBM_mask2_z, IBM_mask3_z
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))  :: IBM_mask1_y, IBM_mask2_y, IBM_mask3_y

        integer         :: k,s, xdmf_id, ierr, anim_id


        character(20), dimension(3) :: global_fields
        character(20), dimension(4) :: frame_fields
        integer, save               :: nbexec=0

        xdmf_id=513

        nbexec=nbexec+1


        global_fields(1)="IBM_mask1"
        global_fields(2)="IBM_mask2"
        global_fields(3)="IBM_mask3"

        frame_fields(1)="q1"
        frame_fields(2)="q2"
        frame_fields(3)="q3"
        frame_fields(4)="pr"

        if (nbexec==1) then
            call open(trim(anim_dir), trim("anim1"), param_anim2D_3%nb_steps, (/xsize(1), xsize(2)/), anim_id)
            call declare_globalfields(anim_id, global_fields, 3)
            call declare_framefields(anim_id, frame_fields, 4)

            call see(anim_id)
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! ACTIVATE this part if you want to export animation to various .h5 files !!!!
        !!!! AVOIDS problems while computing on various NODES !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! USE the appropriate write_xdmf2 subroutine !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (anim2D_3==2) then
            write(tmp_str, "(i10)")nbsteps
            file_name='Frame'//trim(adjustl(tmp_str))
            frame_name=''
            file_path=trim(anim_dir)//'X3/'//trim(file_name)
            if(nrank==0)  call hdf_create_file_with_2Dmesh(trim(file_path), Z, Y, "x1", "x2", n1, n2)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if ((file_nb==1).and.(nbsteps==1)) then

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! ACTIVATE this part if you want to export animation to a single .h5 file !!!!
        !!!! You also need to activate the creation of the Frame group !!!!!!!!!!!!!!!!!!
        !!!! MIGHT cause problems while computing on various NODES !!!!!!!!!!!!!!!!!!!!!!
        !!!! USE the appropriate write_xdmf2 subroutine !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(anim2D_3==1) then
                file_name="anim1"
                file_path=trim(anim_dir)//'X3/'//trim(file_name)
                if(nrank==0)  call hdf_create_file_with_2Dmesh(trim(file_path), Z, Y, "x1", "x2", n1, n2)
            endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if(nrank==0)  call hdf_addgroup(file_path, "MASKS")

            if (IBM_activated) then

                call transpose_x_to_y(IBM_mask1, IBM_mask1_y)
                call transpose_y_to_z(IBM_mask1_y, IBM_mask1_z)

                call transpose_x_to_y(IBM_mask2, IBM_mask2_y)
                call transpose_y_to_z(IBM_mask2_y, IBM_mask2_z)

                call transpose_x_to_y(IBM_mask3, IBM_mask3_y)
                call transpose_y_to_z(IBM_mask3_y, IBM_mask3_z)

                do s = 1, param_anim2D_3%nb_slices

                    k=param_anim2D_3%slices(s)
                    write(tmp_str, "(i10)")k
                    slice_name='MASKS/'//'slice_'//trim(adjustl(tmp_str))
                    if (nrank==0) write(*,*)'ADD GROUPE 1', slice_name
                    if(nrank==0)  call hdf_addgroup(file_path, slice_name)

                    if (param_anim2D_3%export_q1) then
                        call hdf_add_2Dfield(trim(file_path), IBM_mask1_z(:, :,k), trim(slice_name)//"/q1", n1, n2, zstart(1),zend(1), zstart(2),zend(2))
                    end if

                    if (param_anim2D_3%export_q2) then
                        call hdf_add_2Dfield(trim(file_path), IBM_mask2_z(:, :,k), trim(slice_name)//"/q2", n1, n2, zstart(1),zend(1), zstart(2),zend(2))
                    end if

                    if (param_anim2D_3%export_q3) then
                        call hdf_add_2Dfield(trim(file_path), IBM_mask3_z(:, :,k), trim(slice_name)//"/q3", 1, n2, 1,1, zstart(2),zend(2))
                    end if

                end do

            endif

        end if


        ! ADDING THE CURRENT STEP
        if (nrank==0) write(*,*) 'ADDING FRAME n°', nbsteps, 'in file : ', trim(file_path)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! ACTIVATE this part if you want to export animation to a single .h5 file !!!!
        !!!! MIGHT cause problems while computing on various NODES !!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (anim2D_3==1) then
            write(tmp_str, "(i10)")nbsteps
            frame_name='Frame'//trim(adjustl(tmp_str))//'/'
            if (nrank==0) write(*,*)'ADD_GROUPE', frame_name
            if(nrank==0)  call hdf_addgroup(file_path, frame_name)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do s = 1, param_anim2D_3%nb_slices

            k=param_anim2D_3%slices(s)
            write(tmp_str, "(i10)")k
            slice_name=trim(frame_name)//'slice_'//trim(adjustl(tmp_str))
            !if (nrank==0) write(*,*)'ADD_GROUPE', slice_name
            if(nrank==0)  call hdf_addgroup(file_path, slice_name)


            if (param_anim2D_3%export_q1) then
                call hdf_add_2Dfield(trim(file_path), q1_z(:, :,k), trim(slice_name)//"/q1", n1, n2, zstart(1),zend(1), zstart(2),zend(2))
            end if

            if (param_anim2D_3%export_q2) then
                call hdf_add_2Dfield(trim(file_path), q2_z(:, :,k), trim(slice_name)//"/q2", n1, n2, zstart(1),zend(1), zstart(2),zend(2))
            end if

            if (param_anim2D_3%export_q3) then
                call hdf_add_2Dfield(trim(file_path), q3_z(:, :,k), trim(slice_name)//"/q3", n1, n2, zstart(1),zend(1), zstart(2),zend(2))
            end if

            if (param_anim2D_3%export_pr) then
                call hdf_add_2Dfield(trim(file_path), dp_z(:, :,k), trim(slice_name)//"/pr", n1, n2, zstart(1),zend(1), zstart(2),zend(2))
            end if

            if ((param_anim2D_3%export_sca).and.(SCA_state/=0)) then
                call hdf_add_2Dfield(trim(file_path), sca_z(:, :,k), trim(slice_name)//"/sca", n1, n2, zstart(1),zend(1), zstart(2),zend(2))
            end if

        end do

        if(anim2D_3==1) call write_xdmf2_single_HDF5
        if(anim2D_3==2) call write_xdmf2_multi_HDF5

        nbsteps=nbsteps+1

!        if (nbsteps>param_anim2D_3%nb_steps) then
!
!
!            nbsteps=1
!            file_nb=file_nb+1
!
!
!            write(tmp_str, "(i10)")file_nb
!            file_name="anim"//trim(adjustl(tmp_str))
!            file_path=trim(anim_dir)//'X3/'//trim(file_name)
!            if(nrank==0)  call hdf_create_file_with_2Dmesh(trim(file_path), Z, Y, "x1", "x2", n1, n2)
!
!
!
!
!            if (IBM_activated) then
!
!                if(nrank==0)  call hdf_addgroup(file_path, "MASKS")
!
!                call transpose_x_to_y(IBM_mask1, IBM_mask1_y)
!                call transpose_y_to_z(IBM_mask1_y, IBM_mask1_z)
!
!                call transpose_x_to_y(IBM_mask2, IBM_mask2_y)
!                call transpose_y_to_z(IBM_mask2_y, IBM_mask2_z)
!
!                call transpose_x_to_y(IBM_mask3, IBM_mask3_y)
!                call transpose_y_to_z(IBM_mask3_y, IBM_mask3_z)
!
!                do s = 1, param_anim2D_3%nb_slices
!
!                    k=param_anim2D_3%slices(s)
!                    write(tmp_str, "(i10)")k
!                    slice_name='MASKS/'//'slice_'//trim(adjustl(tmp_str))
!                    if(nrank==0)  call hdf_addgroup(file_path, slice_name)
!
!                    if (param_anim2D_3%export_q1) then
!                        call hdf_add_2Dfield(trim(file_path), IBM_mask1_z(:, :,k), trim(slice_name)//"/q1", n1, n2, zstart(1),zend(1), zstart(2),zend(2))
!                    end if
!
!                    if (param_anim2D_3%export_q2) then
!                        call hdf_add_2Dfield(trim(file_path), IBM_mask2_z(:, :,k), trim(slice_name)//"/q2", n1, n2, zstart(1),zend(1), zstart(2),zend(2))
!                    end if
!
!                    if (param_anim2D_3%export_q3) then
!                        call hdf_add_2Dfield(trim(file_path), IBM_mask3_z(:, :,k), trim(slice_name)//"/q3", n1, n2, zstart(1),zend(1), zstart(2),zend(2))
!                    end if
!
!                end do
!
!            end if
!
!        end if

    contains


        subroutine write_xdmf2_single_HDF5()
            use file_copy
            use DNS_settings, only:dt
            use run_ctxt_data, only:tinit
            implicit none
            integer     :: f
            character(200)  :: xdmf_path, attribute_path
            character(20)   :: frame_label, attribute_name
            character(80)   :: tmp_str, tmp_str2, tmp_str3, tmp_str4

            write(tmp_str, "(i10)")file_nb
            xdmf_path=trim(anim_dir)//'X3/anim'//trim(adjustl(tmp_str))//".xdmf"

            if(nrank==0)  call copy_ascii_file(file_dest=trim(xdmf_path), file_src=trim(anim_dir)//"X3/modele2.xdmf")

            if(nrank==0)  open(xdmf_id, file=trim(xdmf_path), position="append")
            if(nrank==0)  write(xdmf_id,'(F10.5,F10.5,I4)')tinit, dt*param_anim2D_1%step_size, param_anim2D_3%nb_steps
            if(nrank==0)  write(xdmf_id,*)"</DataItem>"
            if(nrank==0)  write(xdmf_id,*)"</Time>"
            if(nrank==0)  write(xdmf_id,*)


            do f = 1, nbsteps

                write(tmp_str, "(i10)")f
                frame_label='Frame'//trim(adjustl(tmp_str))
                if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Grid Name="', trim(frame_label), '" GridType="Uniform">'
                if(nrank==0)  write(xdmf_id,*)

                if(nrank==0)  write(xdmf_id,'(a,I4,I4,a)')'<Topology TopologyType="2DRectMesh" NumberOfElements="', n2, n1,'"/>'
                if(nrank==0)  write(xdmf_id,*)'<Geometry GeometryType="VXVY">'

                if(nrank==0)  write(xdmf_id,'(a,I4,a)')'<DataItem Dimensions="',n2,'" NumberType="Float" Precision="8" Format="HDF">'
                if(nrank==0)  write(xdmf_id,*) trim(file_name)//".h5:/x2"
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"

                if(nrank==0)  write(xdmf_id,'(a,I4,a)')'<DataItem Dimensions="',n1,'" NumberType="Float" Precision="8" Format="HDF">'
                if(nrank==0)  write(xdmf_id,*)trim(file_name)//".h5:/x1"
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"


                if(nrank==0)  write(xdmf_id,*)"</Geometry>"

                if(nrank==0)  write(xdmf_id,*)

                do s = 1, param_anim2D_3%nb_slices

                    k=param_anim2D_3%slices(s)
                    write(tmp_str, "(i10)")k
                    slice_name=trim(frame_label)//'/slice_'//trim(adjustl(tmp_str))

                    if(nrank==0)  write(xdmf_id,*)

                    if (param_anim2D_3%export_q1) then
                        attribute_name='q1_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q1'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n2,n1, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
                    end if

                    if (param_anim2D_3%export_q2) then
                        attribute_name='q2_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q2'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n2,n1, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"

                    end if

                    if (param_anim2D_3%export_q3) then
                        attribute_name='q3_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q3'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2),n2,n1, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
                    end if

                    if (param_anim2D_3%export_pr) then
                        attribute_name='pr_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/pr'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2),n2,n1, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
                    end if

                    if ((param_anim2D_3%export_sca).and.(SCA_state/=0)) then
                        attribute_name='sca_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/sca'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2),n2,n1, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
                    end if


                end do


                if(nrank==0)  write(xdmf_id,*)
                if(nrank==0)  write(xdmf_id,*)
                if(nrank==0)  write(xdmf_id,*)"</Grid>"
                if(nrank==0)  write(xdmf_id,*)
            end do


            if(nrank==0)  write(xdmf_id,*)"</Grid>"

            ! WRITTING MASK

            if (IBM_activated) then

                frame_label='MASKS'
                if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Grid Name="', trim(frame_label), '" GridType="Uniform">'
                if(nrank==0)  write(xdmf_id,*)

                if(nrank==0)  write(xdmf_id,'(a,I4,I4,a)')'<Topology TopologyType="2DRectMesh" NumberOfElements="', n2, n1,'"/>'
                if(nrank==0)  write(xdmf_id,*)'<Geometry GeometryType="VXVY">'

                if(nrank==0)  write(xdmf_id,'(a,I4,a)')'<DataItem Dimensions="',n2,'" NumberType="Float" Precision="8" Format="HDF">'
                if(nrank==0)  write(xdmf_id,*) trim(file_name)//".h5:/x2"
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"

                if(nrank==0)  write(xdmf_id,'(a,I4,a)')'<DataItem Dimensions="',n1,'" NumberType="Float" Precision="8" Format="HDF">'
                if(nrank==0)  write(xdmf_id,*)trim(file_name)//".h5:/x1"
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"


                if(nrank==0)  write(xdmf_id,*)"</Geometry>"

                if(nrank==0)  write(xdmf_id,*)

                do s = 1, param_anim2D_3%nb_slices

                    k=param_anim2D_3%slices(s)
                    write(tmp_str, "(i10)")k
                    slice_name=trim(frame_label)//'/slice_'//trim(adjustl(tmp_str))

                    if(nrank==0)  write(xdmf_id,*)

                    if (param_anim2D_3%export_q1) then
                        attribute_name='q1_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q1'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n2,n1, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
                    end if

                    if (param_anim2D_3%export_q2) then
                        attribute_name='q2_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q2'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n2,n1, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"

                    end if

                    if (param_anim2D_3%export_q3) then
                        attribute_name='q3_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q3'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2),n2,n1, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
                    end if

                end do


                if(nrank==0)  write(xdmf_id,*)
                if(nrank==0)  write(xdmf_id,*)
                if(nrank==0)  write(xdmf_id,*)"</Grid>"
                if(nrank==0)  write(xdmf_id,*)

            endif


            if(nrank==0)  write(xdmf_id,*)
            if(nrank==0)  write(xdmf_id,*)"</Domain>"
            if(nrank==0)  write(xdmf_id,*)"</Xdmf>"
            if(nrank==0)  close(xdmf_id)

        end subroutine write_xdmf2_single_HDF5

        subroutine write_xdmf2_multi_HDF5()
            use file_copy
            use DNS_settings, only:dt
            use run_ctxt_data, only:tinit
            implicit none
!            integer     :: f
            character(200)  :: xdmf_path, attribute_path
            character(20)   :: frame_label, attribute_name
            character(80)   :: tmp_str, tmp_str2, tmp_str3, tmp_str4
            integer,save    :: nline_file

            write(tmp_str, "(i10)")file_nb
            xdmf_path=trim(anim_dir)//'X3/anim'//trim(adjustl(tmp_str))//".xdmf"

            if ((nbsteps==1).and.(nrank==0))  call copy_ascii_file(file_dest=trim(xdmf_path), &
                                                                   file_src=trim(anim_dir)//"X3/modele2.xdmf",nline=nline_file)
            if(nrank==0)  open(xdmf_id, file=trim(xdmf_path), position="rewind")

            if (nbsteps==1) then
                !!!! SKIPPING the already existing text !!!!
                do i=1,nline_file
                    if(nrank==0) read(xdmf_id,*)
                enddo

                if(nrank==0)  write(xdmf_id,'(F10.5,F10.5,I4)')tinit, dt*param_anim2D_1%step_size, param_anim2D_3%nb_steps
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"
                if(nrank==0)  write(xdmf_id,*)"</Time>"
                if(nrank==0)  write(xdmf_id,*)

                nline_file=nline_file+4
            else
                !!!! SKIPPING the already existing text, excepting the last lines (Grid, Domain and Xdmf) !!!!
                do i=1,nline_file
                    if(nrank==0) read(xdmf_id,*)
                enddo
            endif

!            do f = 1, nbsteps

                write(tmp_str, "(i10)")nbsteps
                frame_label='Frame'//trim(adjustl(tmp_str))
                if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Grid Name="', trim(frame_label), '" GridType="Uniform">'
                if(nrank==0)  write(xdmf_id,*)

                if(nrank==0)  write(xdmf_id,'(a,I4,I4,a)')'<Topology TopologyType="2DRectMesh" NumberOfElements="', n2, n1,'"/>'
                if(nrank==0)  write(xdmf_id,*)'<Geometry GeometryType="VXVY">'

                if(nrank==0)  write(xdmf_id,'(a,I4,a)')'<DataItem Dimensions="',n2,'" NumberType="Float" Precision="8" Format="HDF">'
                if(nrank==0)  write(xdmf_id,*) trim(file_name)//".h5:/x2"
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"

                if(nrank==0)  write(xdmf_id,'(a,I4,a)')'<DataItem Dimensions="',n1,'" NumberType="Float" Precision="8" Format="HDF">'
                if(nrank==0)  write(xdmf_id,*)trim(file_name)//".h5:/x1"
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"


                if(nrank==0)  write(xdmf_id,*)"</Geometry>"

                if(nrank==0)  write(xdmf_id,*)

                nline_file=nline_file+12

                do s = 1, param_anim2D_3%nb_slices

                    k=param_anim2D_3%slices(s)
                    write(tmp_str, "(i10)")k
                    slice_name='slice_'//trim(adjustl(tmp_str))

                    if(nrank==0)  write(xdmf_id,*)

                    nline_file=nline_file+1

                    if (param_anim2D_3%export_q1) then
                        attribute_name='q1_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q1'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n2,n1, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"

                        nline_file=nline_file+3
                    end if

                    if (param_anim2D_3%export_q2) then
                        attribute_name='q2_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q2'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n2,n1, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"

                        nline_file=nline_file+3
                    end if

                    if (param_anim2D_3%export_q3) then
                        attribute_name='q3_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q3'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2),n2,n1, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"

                        nline_file=nline_file+3
                    end if

                    if (param_anim2D_3%export_pr) then
                        attribute_name='pr_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/pr'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2),n2,n1, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"

                        nline_file=nline_file+3
                    end if

                    if ((param_anim2D_3%export_sca).and.(SCA_state/=0)) then
                        attribute_name='sca_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/sca'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2),n2,n1, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"

                        nline_file=nline_file+3
                    end if


                end do


                if(nrank==0)  write(xdmf_id,*)
                if(nrank==0)  write(xdmf_id,*)
                if(nrank==0)  write(xdmf_id,*)"</Grid>"
                if(nrank==0)  write(xdmf_id,*)

                nline_file=nline_file+4
!            end do


            if(nrank==0)  write(xdmf_id,*)"</Grid>"

            ! WRITTING MASK

            if (IBM_activated) then

                frame_label='MASKS'
                if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Grid Name="', trim(frame_label), '" GridType="Uniform">'
                if(nrank==0)  write(xdmf_id,*)

                if(nrank==0)  write(xdmf_id,'(a,I4,I4,a)')'<Topology TopologyType="2DRectMesh" NumberOfElements="', n2, n1,'"/>'
                if(nrank==0)  write(xdmf_id,*)'<Geometry GeometryType="VXVY">'

                if(nrank==0)  write(xdmf_id,'(a,I4,a)')'<DataItem Dimensions="',n2,'" NumberType="Float" Precision="8" Format="HDF">'
                if(nrank==0)  write(xdmf_id,*) trim(file_name)//".h5:/x2"
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"

                if(nrank==0)  write(xdmf_id,'(a,I4,a)')'<DataItem Dimensions="',n1,'" NumberType="Float" Precision="8" Format="HDF">'
                if(nrank==0)  write(xdmf_id,*)trim(file_name)//".h5:/x1"
                if(nrank==0)  write(xdmf_id,*)"</DataItem>"


                if(nrank==0)  write(xdmf_id,*)"</Geometry>"

                if(nrank==0)  write(xdmf_id,*)

                do s = 1, param_anim2D_3%nb_slices

                    k=param_anim2D_3%slices(s)
                    write(tmp_str, "(i10)")k
                    slice_name=trim(frame_label)//'/slice_'//trim(adjustl(tmp_str))

                    if(nrank==0)  write(xdmf_id,*)

                    if (param_anim2D_3%export_q1) then
                        attribute_name='q1_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q1'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n2,n1, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
                    end if

                    if (param_anim2D_3%export_q2) then
                        attribute_name='q2_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q2'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2), n2,n1, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"

                    end if

                    if (param_anim2D_3%export_q3) then
                        attribute_name='q3_'//trim(adjustl(tmp_str))
                        attribute_path=trim(file_name)//".h5:"//trim(slice_name)//'/q3'

                        tmp_str2='<DataItem Dimensions="'
                        tmp_str3='  " NumberType="Float" Precision="8" Format="HDF">'
                        tmp_str4='</DataItem>'
                        if(nrank==0)  write(xdmf_id,'(a,a,a)')'<Attribute Name="', trim(attribute_name), '" AttributeType="Scalar" Center="Node">'
                        if(nrank==0)  write(xdmf_id,'(a,I4,I4,a,a,a)')trim(tmp_str2),n2,n1, trim(tmp_str3), trim(attribute_path), trim(tmp_str4)
                        if(nrank==0)  write(xdmf_id,*)"</Attribute>"
                    end if

                end do


                if(nrank==0)  write(xdmf_id,*)
                if(nrank==0)  write(xdmf_id,*)
                if(nrank==0)  write(xdmf_id,*)"</Grid>"
                if(nrank==0)  write(xdmf_id,*)

            endif


            if(nrank==0)  write(xdmf_id,*)
            if(nrank==0)  write(xdmf_id,*)"</Domain>"
            if(nrank==0)  write(xdmf_id,*)"</Xdmf>"
            if(nrank==0)  close(xdmf_id)

        end subroutine write_xdmf2_multi_HDF5

    end subroutine anim2D_addframe3

end module anim2D_writer


