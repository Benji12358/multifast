module SCALAR_solver
    use decomp_2d
    use mesh
    use boundaries
    use schemes3D_interface
    use irregular_derivative_coefficients
    use DNS_settings
    use SCALAR_data
    use physical_fields

    implicit none

    real*8, dimension(:,:,:), allocatable       :: previousRHS1, previousRHS2, previoussource_term
    real*8, dimension(:,:), allocatable, save :: diffo, conv4o

    real*8, dimension(:,:,:), allocatable     :: Fextsca, previousFextsca

    real*8  :: gam,rom


    logical, private     :: previousRHS_are_available=.false.

contains

    subroutine init()

        use run_ctxt_data
        use SCALAR_workspace_view, only: recovery_RHS_dir

        implicit none


        allocate(previousRHS1(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
        previousRHS1=0.d0
        allocate(previousRHS2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
        previousRHS2=0.d0
        allocate(previoussource_term(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
        previoussource_term=0.d0
        
        allocate(Fextsca(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
        Fextsca=0.d0
        allocate(previousFextsca(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
        previousFextsca=0.d0

        if ((run_ctxt==CONTINUE_FROM_PREVIOUS_RUN).or.(run_ctxt==RECOVERY_A_RUN)) then

            inquire( file=trim(recovery_RHS_dir)//"/RHS1.h5", exist=previousRHS_are_available)
            if (previousRHS_are_available) then
                call read_previousRHS
            end if
        endif

        contains

            subroutine read_previousRHS()

                use HDF5_IO

                use COMMON_fieldstools

                implicit none
                real*8  :: cs1, cs2
                character(200)    :: file_path


                call perform_checksum_y(previousRHS1, cs1)
                call perform_checksum_y(previousRHS2, cs2)

                if (nrank==0) write(*,*)'SCALAR readpreviousRHS A:', cs1, cs2


                file_path=trim(recovery_RHS_dir)//"/RHS1"
                call hdf_read_3Dfield(file_path, previousRHS1, "RHS1", nx_global,ny_global,nz_global, ystart(1),yend(1),ystart(2),yend(2),ystart(3),yend(3))

                file_path=trim(recovery_RHS_dir)//"/RHS2"
                call hdf_read_3Dfield(file_path, previousRHS2, "RHS2", nx_global,ny_global,nz_global, ystart(1),yend(1),ystart(2),yend(2),ystart(3),yend(3))

                call perform_checksum_y(previousRHS1, cs1)
                call perform_checksum_y(previousRHS2, cs2)

                if (nrank==0) write(*,*)'SCALAR readpreviousRHS B:', cs1, cs2


            end subroutine read_previousRHS

    end subroutine init

    subroutine solve_scalar(ns)

        use time_schemes
        use run_ctxt_data
        use inflow_settings
        use numerical_methods_settings

        use mesh, only: Yc
        use DNS_settings, only: ren

        use SCALAR_data

        use SCALAR_inout_flow_old, only:    &
        OPEN_OLD_add_outflow=>add_outflow,     &
        OPEN_OLD_get_inflow=>get_inflow

        use SCALAR_inout_flow, only:    &
        OPEN_add_outflow=>add_outflow,     &
        OPEN_get_inflow=>get_inflow, &
        OPEN_newversion=>inout_newversion

        use IBM_settings, only: IBM_activated

        implicit none

        integer, intent(in) :: ns

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: sca_y0

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: q1S   ! scalar * Q1
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: q2S   ! scalar * Q2
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: q3S   ! scalar * Q3
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: RHS
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: RHS1
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: RHS2, RHS2_1
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: source_term

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: u_x
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: u_y

        real*8                                                  :: pr

        if (.not. inout_newversion) then

            if (outflow_buff>0) call OPEN_OLD_add_outflow(ntime, ns)
            if (BC1==OPEN) then
                call OPEN_OLD_get_inflow(ntime, ns)
            endif

        else

            if (BC1==OPEN) then
                call OPEN_get_inflow()
            endif

        endif

        call set_time_coeffs
        call perform_velocityscalar_products(q1_x, q2_y, q3_z, q1S, q2S, q3S)
        call perform_RHS1(q1S, q2S, q3S, RHS1)
        call perform_RHS2(RHS2)

        if (init_type==KAWAMURA_INIT) then
            ! If Kawamura equations are used, we need to compute a source term for the transport of the passive scalar
            ! Check Kawamura 1998 (doi: 10.1016/S0142-727X(98)10026-7) and Kasagi 1992 (doi: 10.1115/1.2911323)

            ! Caution ! 
            ! Works only for streamwise = 1
            call D0s_3Dx(q1_x, u_x, n1,xsize(2),xsize(2),xsize(3),xsize(3), NS_DEF_BC1)
            call transpose_x_to_y(u_x,u_y)
            call perform_source_term(u_x, u_y, source_term, dx3)

        else

            source_term = 0.d0

        endif

        ! RHS=RHS1+RHS2

        if (BC1==OPEN) call outflowing
        call next_scalar

        if (init_type==CONSTANT_HEAT_FLUX) call next_wall_scalar

        previousRHS1=RHS1
        previousRHS2=RHS2
        previoussource_term=source_term
        previousFextsca=Fextsca 

        call transpose_y_to_x(sca_y, sca_x)
        call transpose_y_to_z(sca_y, sca_z)

        if ((inout_newversion).and.(outflow_buff>0)) call OPEN_add_outflow

    contains

        subroutine set_time_coeffs()

            logical, save :: first_time_advancement=.true.

!            previousRHS_are_available=.false.

            if (first_time_advancement) then

                if (previousRHS_are_available) then

                    gam=ga(ns)
                    rom=ro(ns)

                else
                    gam=1.d0
                    rom=0.d0

                end if

            else
                gam=ga(ns)
                rom=ro(ns)

            end if

            first_time_advancement=.false.

        end subroutine set_time_coeffs

        subroutine next_scalar
            use IBM_data
            use IBM_settings
            use IBM

            implicit none
            integer :: i,j,k
            integer :: n1s, n1e, n2s,n2e, n3s,n3e
            real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: sca_term_x
            real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: sca_term
            real*8  :: RHSi

            if (BC1==OPEN) then
                n1s=max(2, ystart(1))
                n1e=min(n1-2, yend(1))
            else
                n1s=ystart(1)
                n1e=min(n1-1, yend(1))
            endif
            n2s=1
            n2e=n2-1

            n3s=ystart(3)
            n3e=min(n3-1, yend(3))

            if (IBM_activated) then

                call compute_antisymmetric_temperature(sca_x, sca_y, sca_z, xsize, sca_term_x, delta_T, init_type)
                call transpose_x_to_y(sca_term_x, sca_term)

                if (init_type.eq.CLASSIC_INIT) then
                ! For classic init, Lyons problem to solve with fixed temperatures at walls
                ! We imposed constant temperature in objects and do a "fadlun-like" extrapolation in fluid domain
                ! Thus the mask contains all points in the objects + the first nodes in the fluid domain

                    do k = n3s, n3e
                        do j = n2s, n2e
                            do i = n1s, n1e

                                RHSi = dt*gam*(RHS1(i,j,k)+RHS2(i,j,k)+source_term(i,j,k)+Fextsca(i,j,k))+dt*rom*(previousRHS1(i,j,k)+previousRHS2(i,j,k)+previoussource_term(i,j,k)+previousFextsca(i,j,k))
                                sca_y(i,j,k)=sca_y(i,j,k) + RHSi + (mask_objects_sca_y(i,j,k)+mask_fadlun_sca_y(i,j,k))*(-RHSi + sca_term(i,j,k))

                            end do
                        end do
                    end do

                else if (init_type.eq.KAWAMURA_INIT) then
                ! For Kawamura init, we need to impose theta = 0 on surfaces.
                ! Thus, we use only the first nodes inside the objects as for the velocity.

                ! For now, only the use of fadlun implementation works for the scalar. 
                ! Using Kim implementaion is leading to accumulation problem. Maybe an issue with the mass source/sink term of Kim ...

                    do k = n3s, n3e
                        do j = n2s, n2e
                            do i = n1s, n1e

                                RHSi = dt*gam*(RHS1(i,j,k)+RHS2(i,j,k)+source_term(i,j,k)+Fextsca(i,j,k))+dt*rom*(previousRHS1(i,j,k)+previousRHS2(i,j,k)+previoussource_term(i,j,k)+previousFextsca(i,j,k))
                                ! sca_y(i,j,k)=sca_y(i,j,k) + RHSi + mask_kim_sca_y(i,j,k)*(-RHSi + sca_term(i,j,k))
                                sca_y(i,j,k)=sca_y(i,j,k) + RHSi + (mask_objects_sca_y(i,j,k)+mask_fadlun_sca_y(i,j,k))*(-RHSi + sca_term(i,j,k))

                            end do
                        end do
                    end do

                else
                ! For type of init, i.e. problems, the implementation has still to be done

                if (nrank.eq.0) then
                    write(*,*) "This type of problem (for the temperature) is not implemented yet with the IBM"
                    call exit
                endif

        endif

            else

                do k = n3s, n3e
                    do j = n2s, n2e
                        do i = n1s, n1e
                            sca_y(i,j,k)=sca_y(i,j,k)+dt*gam*(RHS1(i,j,k)+RHS2(i,j,k)+source_term(i,j,k)+Fextsca(i,j,k))+dt*rom*(previousRHS1(i,j,k)+previousRHS2(i,j,k)+previoussource_term(i,j,k)+previousFextsca(i,j,k))

                            ! sca_y(i,j,k)=sca_y(i,j,k)+dt*(gam*RHS1(i,j,k)+rom*previousRHS1(i,j,k))+dt*(gam*RHS2(i,j,k)+rom*previousRHS2(i,j,k))
                            !if (sca_y(i,j,k)<0) sca_y(i,j,k)=0.d0
                            !if (sca_y(i,j,k)>1) sca_y(i,j,k)=1.d0
                        end do
                    end do
                end do

            endif

        end subroutine next_scalar

        subroutine next_wall_scalar
            implicit none
            integer :: i,j,k
            integer :: n1s, n1e, n3s,n3e

            if (BC1==OPEN) then
                n1s=max(2, ystart(1))
                n1e=min(n1-2, yend(1))
            else
                n1s=ystart(1)
                n1e=min(n1-1, yend(1))
            endif

            n3s=ystart(3)
            n3e=min(n3-1, yend(3))

            do k = n3s, n3e
                do i = n1s, n1e

                        sca_wall20(i,k) = sca_y(i,1,k) + heat_flux * Yc(1)
                        sca_wall21(i,k) = sca_y(i,n2-1,k) + heat_flux * Yc(1)

                end do
            end do

        end subroutine next_wall_scalar

    end subroutine solve_scalar

    subroutine outflowing()

        use MPI
        use DNS_settings

        use SCALAR_data

        implicit none
        integer, save                                           :: nb=1
        ! integer, parameter                                      :: outflow_type=2
        integer                                                 :: n2s,n2e, n3s, n3e, j, k
        real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: cx4
        real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: conv4

        real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: diff

        nb=nb+1

        diff=0.d0

        if (outflow_type==1) then
            cx4=dt/dx1
        endif

        if (outflow_type==2) then
            call perform_cx
        endif

        if (outflow_type==3) then
            call perform_cx3
        endif

        if (outflow_type==4) then
            call perform_cx3
            call perform_diff
        endif


        call perform_conv

        n2s=xstart(2)
        n3s=xstart(3)
        n2e=xend(2)
        n3e=xend(3)

        do k=n3s,n3e
            do j=n2s,n2e
                sca_x(n1-1,j,k)=sca_x(n1-1,j,k)+conv4(j,k)+diff(j,k)*dt
                !if(sca_x(n1-1,j,k)<0) sca_x(n1-1,j,k)=0.d0
                !if(sca_x(n1-1,j,k)>1) sca_x(n1-1,j,k)=1.d0
            enddo
        enddo

        call transpose_x_to_y(sca_x(:,:,:), sca_y(:,:,:))

    contains

        subroutine perform_cx()
            implicit none
            real*8  :: num, den

            do k=xstart(3),min(xend(3), n3-1)
                do j=xstart(2),min(xend(2), n2-1)

                    cx4(j,k)=sum(q1_x(3:n1-4, j,k))
                    cx4(j,k)=cx4(j,k)/(n1-6)
                    cx4(j,k)=cx4(j,k)*dt/dx1

                enddo
            enddo

        end subroutine

        subroutine perform_cx2()
            implicit none
            real*8  :: num, den

            do k=xstart(3),min(xend(3), n3-1)
                do j=xstart(2),min(xend(2), n2-1)

                    cx4(j,k)=q1_x(n1-1, j,k)
                    cx4(j,k)=cx4(j,k)*dt/dx1

                enddo
            enddo


        end subroutine perform_cx2

        subroutine perform_cx3()
            implicit none
            real*8  :: num, den
            real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: q1_cc

            call D0ssh_3Dx(q1_x, q1_cc, n1,xsize(2),xsize(2),xsize(3),xsize(3), NS_Q3_BC1)

            do k=xstart(3),min(xend(3), n3-1)
                do j=xstart(2),min(xend(2), n2-1)

                    cx4(j,k)=q1_cc(n1-1, j,k)
                    cx4(j,k)=cx4(j,k)*dt/dx1

                enddo
            enddo


        end subroutine perform_cx3

        subroutine perform_diff()
            use boundary_scheme
            use snapshot_writer
            implicit none
            integer :: i,j,k

            real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: diffc
            real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: diff3
            real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: diff2
            real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: diff3y
            real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: diff3x, diff2x

            if (.not. allocated(diffo)) then
                allocate(diffo(xstart(2):xend(2), xstart(3):xend(3)))

                diffo=0.d0
            endif

            call D2c_3Dz(sca_z, diff3, zsize(1),zsize(1),zsize(2),zsize(2),n3, dx3*dsqrt(ren*prandtl), .true., NS_Q1_BC3)

            call transpose_z_to_y(diff3, diff3y)
            call transpose_y_to_x(diff3y, diff3x)

            call D1c_MULT_3Dy(sca_y, diff2, ysize(1),ysize(1),n2,ysize(3),ysize(3), dx2*ren*prandtl, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,1))
            call D2c_MULTACC_3Dy(sca_y, diff2, ysize(1),ysize(1),n2,ysize(3),ysize(3), dx2*dsqrt(ren*prandtl), .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,2))


            do k=ystart(3), min(n3m, yend(3))       !do k=1,n3m
                do i=ystart(1), min(n1m, yend(1))       !do i=1,n1m
                    diff2(i,1,k)=( sca_y(i,2,k)*a3_d + sca_y(i,1,k)*a2_d + sca_wall20(i,k)*a1_d)            /(dx2*dsqrt(ren*prandtl))**2
                    diff2(i,n2m,k)=( sca_y(i,n2m-1,k)*a1_u + sca_y(i,n2m,k)*a2_u + sca_wall21(i,k)*a3_u ) /(dx2*dsqrt(ren*prandtl))**2
                enddo
            enddo


            call transpose_y_to_x(diff2, diff2x)

            diffc=diff2x(n1-1,:,:)+diff3x(n1-1,:,:)

            diff=gam*diffc+rom*diffo

            diffo=diffc

        end subroutine perform_diff

        subroutine perform_conv()
            use snapshot_writer
            implicit none
            real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: conv4c

            if (.not. allocated(conv4o)) then

                allocate(conv4o(xstart(2):xend(2), xstart(3):xend(3)))
                conv4o=0.d0

            endif

            conv4c=0.d0
            do k=xstart(3),min(xend(3), n3-1)
                do j=xstart(2),min(xend(2), n2-1)
                    conv4c(j,k)=-cx4(j,k)*(sca_x(n1-1,j,k)-sca_x(n1-2,j,k))
                enddo
            enddo

            conv4=gam*conv4c+rom*conv4o

            conv4o=conv4c


        end subroutine perform_conv
    end subroutine outflowing

    subroutine perform_RHS1(q1S, q2S, q3S, RHS1)

        use SCALAR_data

        implicit none

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent (in) :: q1S   ! scalar * Q1
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent (in) :: q2S   ! scalar * Q2
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent (in) :: q3S   ! scalar * Q3

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent (out)    :: RHS1

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: RHS1x   ! scalar * Q1
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))  :: RHS1y   ! scalar * Q2
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))  :: RHS1z   ! scalar * Q3

        integer :: n1e, n2e, n3e

        RHS1x=0.d0
        RHS1y=0.d0
        RHS1z=0.d0
        RHS1=0.d0

        ! ! Q3*scalar ************************************************
        ! ! velocity interpolation: j (1:n2-1), k(1,n3-1))
        ! n1e=(min(n1-1, zend(1))-zstart(1))+1
        ! n2e=(min(n2-1, zend(2))-zstart(2))+1
        ! call D1c_3Dz(q3S, RHS1z, zsize(1),n1e,zsize(2),n2e,n3, -dx3, .true., TRANSPORT_Q3S_BC3)


        ! ! Q2*scalar ************************************************
        ! ! velocity interpolation: i (1:n1-1), k(1,n3-1)
        ! n1e=(min(n1-1, yend(1))-ystart(1))+1
        ! n3e=(min(n3-1, yend(3))-ystart(3))+1


        ! call D1c_MULTACC_3Dy(q2S, RHS1y, ysize(1),n1e,n2,ysize(3),n3e, -dx2,.true., TRANSPORT_Q2S_BC2, Yc_to_YcTr_for_D1)

        ! !if (SCA_BC2==FIXED_VALUE) call d1c_wall2
        ! call d1c_wall2

        ! ! Q1*scalar ************************************************
        ! ! velocity interpolation: i (1:n1-1), j(1,n2-1)
        ! n2e=(min(n2m, xend(2))-xstart(2))+1
        ! n3e=(min(n3m, xend(3))-xstart(3))+1
        ! call D1c_3Dx(q1S, RHS1x, n1,xsize(2),n2e,xsize(3),n3e, -dx1, .true., TRANSPORT_Q1S_BC1)

        ! ! RHS1 et RHS2 Table compilation
        ! RHS1=RHS1y
        ! call transpose_x_to_y(RHS1x, RHS1y)
        ! RHS1=RHS1+RHS1y
        ! call transpose_z_to_y(RHS1z, RHS1y)
        ! RHS1=RHS1+RHS1y


        ! NON CONSERVATIVE

        ! Q3*scalar ************************************************
        ! velocity interpolation: j (1:n2-1), k(1,n3-1))
        n1e=(min(n1-1, zend(1))-zstart(1))+1
        n2e=(min(n2-1, zend(2))-zstart(2))+1
        call D1c_3Dz(sca_z, RHS1z, zsize(1),n1e,zsize(2),n2e,n3, -dx3, .true., TRANSPORT_Q3S_BC3)


        ! Q2*scalar ************************************************
        ! velocity interpolation: i (1:n1-1), k(1,n3-1)
        n1e=(min(n1-1, yend(1))-ystart(1))+1
        n3e=(min(n3-1, yend(3))-ystart(3))+1


        call D1c_MULTACC_3Dy(sca_y, RHS1y, ysize(1),n1e,n2,ysize(3),n3e, -dx2,.true., TRANSPORT_Q2S_BC2, Yc_to_YcTr_for_D1)

        !if (SCA_BC2==FIXED_VALUE) call d1c_wall2
        call d1c_wall2

        ! Q1*scalar ************************************************
        ! velocity interpolation: i (1:n1-1), j(1,n2-1)
        n2e=(min(n2m, xend(2))-xstart(2))+1
        n3e=(min(n3m, xend(3))-xstart(3))+1
        call D1c_3Dx(sca_x, RHS1x, n1,xsize(2),n2e,xsize(3),n3e, -dx1, .true., TRANSPORT_Q1S_BC1)

        ! RHS1 et RHS2 Table compilation
        RHS1x = RHS1x*q1S
        RHS1y = RHS1y*q2S
        RHS1z = RHS1z*q3S

        RHS1=RHS1y
        call transpose_x_to_y(RHS1x, RHS1y)
        RHS1=RHS1+RHS1y
        call transpose_z_to_y(RHS1z, RHS1y)
        RHS1=RHS1+RHS1y

    contains

        subroutine d1c_wall2()
            use boundary_scheme
            implicit none
            integer ::i,j,k

            if (SCA_BC2==FIXED_VALUE) then

                ! do k=ystart(3), min(n3m, yend(3))       !do k=1,n3m
                !     do i=ystart(1), min(n1m, yend(1))       !do i=1,n1m
                !         ! RHS1y(i,1,k)=-(q2S(i,1,k)-q2_wall20(i,k)*sca_wall20(i,k))/(Yc(1)-Y(1))
                !         ! RHS1y(i,n2m,k)=-(q2_wall21(i,k)*sca_wall21(i,k)-q2S(i,n2m,k))/(Y(n2)-Yc(n2-1))

                !         RHS1y(i,1,k)= - ( ( 0.5d0 * ( q2S(i,1,k) + q2S(i,2,k) ) - q2_wall20(i,k)*sca_wall20(i,k) )/dx2 ) * Yc_to_YcTr_for_D1(1)
                !         RHS1y(i,n2m,k)= - ( ( q2_wall21(i,k)*sca_wall21(i,k) - 0.5d0 * ( q2S(i,n2-1,k) + q2S(i,n2-2,k) ) )/dx2 ) * Yc_to_YcTr_for_D1(n2m)

                !     enddo
                ! enddo

                ! NON CONSERVATIVE

                do k=ystart(3), min(n3m, yend(3))       !do k=1,n3m
                    do i=ystart(1), min(n1m, yend(1))       !do i=1,n1m
                        ! RHS1y(i,1,k)=-(q2S(i,1,k)-q2_wall20(i,k)*sca_wall20(i,k))/(Yc(1)-Y(1))
                        ! RHS1y(i,n2m,k)=-(q2_wall21(i,k)*sca_wall21(i,k)-q2S(i,n2m,k))/(Y(n2)-Yc(n2-1))

                        RHS1y(i,1,k)= - ( ( 0.5d0 * ( sca_y(i,1,k) + sca_y(i,2,k) ) - sca_wall20(i,k) )/dx2 ) * Yc_to_YcTr_for_D1(1)
                        RHS1y(i,n2m,k)= - ( ( sca_wall21(i,k) - 0.5d0 * ( sca_y(i,n2-1,k) + sca_y(i,n2-2,k) ) )/dx2 ) * Yc_to_YcTr_for_D1(n2m)

                    enddo
                enddo


            end if
            
            if (SCA_BC2==FIXED_FLUX) then
                do k=ystart(3), min(n3m, yend(3))       !do k=1,n3m
                    do i=ystart(1), min(n1m, yend(1))       !do i=1,n1m
                        RHS1y(i,1,k)=-(q2S(i,1,k)-q2_wall20(i,k)*(heat_flux*(Yc(1)-Y(1))+sca_y(i,1,k)))/(Yc(1)-Y(1))
                        !RHS1y(i,n2m,k)=-(q2_wall21(i,k)*(heat_flux*(Yc(1)-Y(1))+sca_y(i,1,k))-q2S(i,n2m,k))/(Y(n2)-Yc(n2-1))
                        RHS1y(i,n2m,k)=-q2_wall21(i,k)*heat_flux
                    enddo
                enddo
            end if

        end subroutine d1c_wall2

    end subroutine perform_RHS1

    subroutine perform_RHS2(RHS2)
        implicit none

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent (out)    :: RHS2

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: RHS2x   ! scalar * Q1
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))  :: RHS2y   ! scalar * Q2
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))  :: RHS2z   ! scalar * Q3

        integer :: n1e, n2e, n3e

        RHS2x=0.d0
        RHS2y=0.d0
        RHS2z=0.d0
        RHS2=0.d0

        ! Diffusive scalar - Z direction ************************************************
        n1e=(min(n1-1, zend(1))-zstart(1))+1
        n2e=(min(n2-1, zend(2))-zstart(2))+1
        call D2c_3Dz(sca_z, RHS2z, zsize(1),n1e,zsize(2),n2e,n3, dx3*dsqrt(ren*prandtl), .true., TRANSPORT_SCA_BC3)


        ! Diffusive scalar - Y direction ************************************************
        n1e=(min(n1-1, yend(1))-ystart(1))+1
        n3e=(min(n3-1, yend(3))-ystart(3))+1


        call D1c_MULTACC_3Dy(sca_y, RHS2y, ysize(1),n1e,n2,ysize(3),n3e, dx2*ren*prandtl,.true., TRANSPORT_SCA_BC2, Yc_to_YcTr_for_D2(:,1))
        call D2c_MULTACC_3Dy(sca_y, RHS2y, ysize(1),n1e,n2,ysize(3),n3e, dx2*dsqrt(ren*prandtl),.true., TRANSPORT_SCA_BC2, Yc_to_YcTr_for_D2(:,2))

        !if (SCA_BC2==FIXED_VALUE) call d2c_wall2
        call d2c_wall2

        ! Diffusive scalar - X direction ************************************************
        n2e=(min(n2m, xend(2))-xstart(2))+1
        n3e=(min(n3m, xend(3))-xstart(3))+1
        call D2c_3Dx(sca_x, RHS2x, n1,xsize(2),n2e,xsize(3),n3e, dx1*dsqrt(ren*prandtl), .true., TRANSPORT_SCA_BC1)

        ! RHS2 Table compilation

        RHS2=RHS2y
        call transpose_x_to_y(RHS2x, RHS2y)
        RHS2=RHS2+RHS2y
        call transpose_z_to_y(RHS2z, RHS2y)
        RHS2=RHS2+RHS2y

    contains

        subroutine d2c_wall2()
            use boundary_scheme
            implicit none
            integer ::i,j,k

            if (SCA_BC2==FIXED_VALUE) then
                do k=ystart(3), min(n3m, yend(3))       !do k=1,n3m
                    do i=ystart(1), min(n1m, yend(1))       !do i=1,n1m
                        RHS2y(i,1,k)=RHS2y(i,1,k)+( sca_y(i,2,k)*a3_d + sca_y(i,1,k)*a2_d + sca_wall20(i,k)*a1_d)            /(dx2*dsqrt(ren*prandtl))**2
                        RHS2y(i,n2m,k)=RHS2y(i,n2m,k)+( sca_y(i,n2m-1,k)*a1_u + sca_y(i,n2m,k)*a2_u + sca_wall21(i,k)*a3_u ) /(dx2*dsqrt(ren*prandtl))**2
                    enddo
                enddo
            end if

            if (SCA_BC2==FIXED_FLUX) then
                do k=ystart(3), min(n3m, yend(3))       !do k=1,n3m
                    do i=ystart(1), min(n1m, yend(1))       !do i=1,n1m
                        RHS2y(i,1,k)=RHS2y(i,1,k)+( sca_y(i,2,k)*a3_d + sca_y(i,1,k)*(a2_d+a1_d) +heat_flux*(Yc(1)-Y(1))*a1_d)            /(dx2*dsqrt(ren*prandtl))**2
                        RHS2y(i,n2m,k)=RHS2y(i,n2m,k)+( sca_y(i,n2m-1,k)*a1_u + sca_y(i,n2m,k)*(a2_u+a3_u) + heat_flux*(Y(n2)-Yc(n2-1))*a3_u ) /(dx2*dsqrt(ren*prandtl))**2
                    enddo
                enddo
            end if

!            open(1516, file="rhs2", position="append")
!            write(1516,*)sum(RHS2y(:,n2m,:))
!            close(1516)

        end subroutine d2c_wall2

    end subroutine perform_RHS2

    subroutine perform_source_term(q_x, q_y, source_term, dx_spanwise)

        use mpi
        use COMMON_workspace_view, only: COMMON_log_path
        use run_ctxt_data, only: ntime

        use IBM_settings, only: IBM_activated, ren_tau_ibm, height_object, number_of_objects, i_start_obj_sca, i_end_obj_sca
        use IBM_data, only: mask_objects_sca_x

        implicit none

        real*8, intent (in) :: dx_spanwise
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent (in) :: q_y
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent (in) :: q_x
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent (out) :: source_term
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: source_term_x
        real*8, dimension(n1) :: ut1, ut1_global, IBM_presence
        real*8 :: umean1, umean_glob1, ren_tau
        integer:: i,j,k,n
        integer     :: mpi_err
        character*200       :: temperature_history_file

        ! If we use IBM, the flowrate is computed at each x-position in the fluid domain
        if (IBM_activated) then

            ! First, we define the positions where an object is present
            IBM_presence=0.d0
            do n=1,number_of_objects

                do i=i_start_obj_sca(n),i_end_obj_sca(n)
                    IBM_presence(i)=1.d0
                enddo

            enddo

            ! Then, we perform the integral over y-direction
            ut1=0.d0
            do i = xstart(1),min(xend(1),n1-1)
                do j=xstart(2),min(xend(2),n2-1)
                    do k=xstart(3),min(xend(3),n3-1)

                        ! Caution, only values of velocity in fluid domain are taken
                        ! And, we divide by the real section of the fluid domain.
                        ! BIG CAUTION !
                        ! This naive implementation only works for staggered or face-to-face roughness elements 
                        ut1(i)=ut1(i)+(1-mask_objects_sca_x(i,j,k))*q_x(i,j,k)*dx_spanwise*cell_size_Y(j)/(L2-IBM_presence(i)*height_object)

                    enddo
                enddo
            enddo

            ut1_global=0.d0
            call MPI_ALLREDUCE(ut1, ut1_global, n1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
            ut1_global=ut1_global/L3

            do i = xstart(1),min(xend(1),n1-1)
                do j=xstart(2),min(xend(2),n2-1)
                    do k=xstart(3),min(xend(3),n3-1)

                        source_term_x(i,j,k) = q_x(i,j,k) / ut1_global(i)

                    enddo
                enddo
            enddo

            call transpose_x_to_y(source_term_x,source_term)

            source_term = source_term * ren_tau_ibm/ren

        else
            ! We compute as usual, bulk velocity at different x-positions and u_tau from viscous drag

            ! Perform the integral over y-direction
            ut1=0.d0
            do i = xstart(1),min(xend(1),n1-1)
                do j=xstart(2),min(xend(2),n2-1)
                    do k=xstart(3),min(xend(3),n3-1)
                        ut1(i)=ut1(i)+q_x(i,j,k)*cell_size_Y(j)*dx_spanwise
                    enddo
                enddo
            enddo

            call MPI_ALLREDUCE(ut1, ut1_global, n1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
            ut1_global=ut1_global/(L3*L2)

            do i = xstart(1),min(xend(1),n1-1)
                do j=xstart(2),min(xend(2),n2-1)
                    do k=xstart(3),min(xend(3),n3-1)

                        source_term_x(i,j,k) = q_x(i,j,k) / ut1_global(i)

                    enddo
                enddo
            enddo

            call transpose_x_to_y(source_term_x,source_term)

            ! compute u_tau
            umean1 = 0.d0

            do k=ystart(3), min(n3m, yend(3))
                do i=ystart(1), min(n1m, yend(1))
                    umean1=umean1+q_y(i,1,k)
                    umean1=umean1+q_y(i,n2-1,k)
                enddo
            enddo

            call MPI_ALLREDUCE (umean1, umean_glob1, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

            umean_glob1=umean_glob1/dfloat(2*n1m*n3m)
            ren_tau=dsqrt(ren*(umean_glob1/Yc(1)))

            source_term = source_term * ren_tau/ren

        endif

        ! In case we want to monitor this source term over time

        ! if (nrank==0) then

        !     write(*,*) 'ut1', sum(ut1)/n1m
        !     write(*,*) 'ut1_global', sum(ut1_global)/n1m
        !     write(*,*) 'ut1_global', ut1_global
        !     write(*,*) 'source_term', sum(source_term)/(n1m*n2m*n3m)

        ! endif

        ! write ren_tau and u_bulk to file

        ! temperature_history_file = trim(COMMON_log_path)//'Velocity/TemperatureKawamura.dat'

        ! open(700, file=trim(temperature_history_file), position="append")

        ! write(700,*) ntime, ut1_global, ren_tau
        ! close(700)

    end subroutine perform_source_term

    subroutine perform_velocityscalar_products(q1_x, q2_y, q3_z, q1S, q2S, q3S)

        use SCALAR_data

        implicit none
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent (in) :: q1_x
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent (in) :: q2_y
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent (in) :: q3_z

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent (out)    :: q1S   ! scalar * Q1
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent (out)    :: q2S   ! scalar * Q2
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent (out)    :: q3S   ! scalar * Q3

        integer :: n1e, n2e, n3e
        integer :: i,j,k

        ! ! Q1*scalar ************************************************
        ! ! velocity interpolation: j (1:n2-1), k(1,n3-1)
        ! n2e=(min(n2-1, xend(2))-xstart(2))+1
        ! n3e=(min(n3-1, xend(3))-xstart(3))+1
        ! q1S=0.d0
        ! call D0s_3Dx(q1_x, q1S, n1,xsize(2),n2e,xsize(3),n3e, NS_Q1_BC1)
        ! q1S=q1S*sca_x


        ! ! Q2*scalar ************************************************
        ! ! velocity interpolation: i (1:n1-1), k(1,n3-1)
        ! n1e=(min(n1-1, yend(1))-ystart(1))+1
        ! n3e=(min(n3-1, yend(3))-ystart(3))+1
        ! q2S=0.d0
        ! call D0s_3Dy(q2_y, q2S, ysize(1),n1e,n2,ysize(3),n3e, NS_Q2_BC2)
        ! q2S=q2S*sca_y


        ! ! Q3*scalar ************************************************
        ! ! velocity interpolation: i (1:n1-1), j(1,n2-1)
        ! n1e=(min(n1-1, zend(1))-zstart(1))+1
        ! n2e=(min(n2-1, zend(2))-zstart(2))+1
        ! q3S=0.d0
        ! call D0s_3Dz(q3_z, q3S, zsize(1),n1e,zsize(2),n2e,n3, NS_Q2_BC3)
        ! q3S=q3S*sca_z

        ! if (IBM_activated) then

        !     ! do i=xstart(1),xend(1)
        !     !     do j=xstart(2),xend(2)
        !     !         do k=xstart(3),xend(3)

        !     !             q1S(i,j,k) = q1S(i,j,k) + IBM_mask_boundscc_x(i,j,k)*(-q1S(i,j,k))

        !     !         enddo
        !     !     enddo
        !     ! enddo

        !     ! do i=ystart(1),yend(1)
        !     !     do j=ystart(2),yend(2)
        !     !         do k=ystart(3),yend(3)

        !     !             q2S(i,j,k) = q2S(i,j,k) + IBM_mask_boundscc_y(i,j,k)*(-q2S(i,j,k))

        !     !         enddo
        !     !     enddo
        !     ! enddo

        !     ! do i=zstart(1),zend(1)
        !     !     do j=zstart(2),zend(2)
        !     !         do k=zstart(3),zend(3)

        !     !             q3S(i,j,k) = q3S(i,j,k) + IBM_mask_boundscc_z(i,j,k)*(-q3S(i,j,k))

        !     !         enddo
        !     !     enddo
        !     ! enddo

        !     q1S=0.d0
        !     q2S=0.d0
        !     q3S=0.d0

        ! endif



        ! NON CONSERVATIVE


        ! Q1*scalar ************************************************
        ! velocity interpolation: j (1:n2-1), k(1,n3-1)
        n2e=(min(n2-1, xend(2))-xstart(2))+1
        n3e=(min(n3-1, xend(3))-xstart(3))+1
        q1S=0.d0
        call D0s_3Dx(q1_x, q1S, n1,xsize(2),n2e,xsize(3),n3e, NS_Q1_BC1)
        ! q1S=q1S*sca_x


        ! Q2*scalar ************************************************
        ! velocity interpolation: i (1:n1-1), k(1,n3-1)
        n1e=(min(n1-1, yend(1))-ystart(1))+1
        n3e=(min(n3-1, yend(3))-ystart(3))+1
        q2S=0.d0
        call D0s_3Dy(q2_y, q2S, ysize(1),n1e,n2,ysize(3),n3e, NS_Q2_BC2)
        ! q2S=q2S*sca_y


        ! Q3*scalar ************************************************
        ! velocity interpolation: i (1:n1-1), j(1,n2-1)
        n1e=(min(n1-1, zend(1))-zstart(1))+1
        n2e=(min(n2-1, zend(2))-zstart(2))+1
        q3S=0.d0
        call D0s_3Dz(q3_z, q3S, zsize(1),n1e,zsize(2),n2e,n3, NS_Q2_BC3)
        ! q3S=q3S*sca_z

        ! if (IBM_activated) then

        !     ! do i=xstart(1),xend(1)
        !     !     do j=xstart(2),xend(2)
        !     !         do k=xstart(3),xend(3)

        !     !             q1S(i,j,k) = q1S(i,j,k) + IBM_mask_boundscc_x(i,j,k)*(-q1S(i,j,k))

        !     !         enddo
        !     !     enddo
        !     ! enddo

        !     ! do i=ystart(1),yend(1)
        !     !     do j=ystart(2),yend(2)
        !     !         do k=ystart(3),yend(3)

        !     !             q2S(i,j,k) = q2S(i,j,k) + IBM_mask_boundscc_y(i,j,k)*(-q2S(i,j,k))

        !     !         enddo
        !     !     enddo
        !     ! enddo

        !     ! do i=zstart(1),zend(1)
        !     !     do j=zstart(2),zend(2)
        !     !         do k=zstart(3),zend(3)

        !     !             q3S(i,j,k) = q3S(i,j,k) + IBM_mask_boundscc_z(i,j,k)*(-q3S(i,j,k))

        !     !         enddo
        !     !     enddo
        !     ! enddo

        !     q1S=0.d0
        !     q2S=0.d0
        !     q3S=0.d0

        ! endif

    end subroutine perform_velocityscalar_products



    subroutine check_solver()
        use mesh
        use snapshot_writer
        use COMMON_workspace_view
        use mathematical_constants
        implicit none
        real*8, parameter   :: T1=0.5d0, T0=-0.5d0
        real*8, dimension(1:n3) :: x3
        real*8, dimension(1:n2) :: x2
        real*8, dimension(1:n1) :: x1

        real*8, dimension(1:n3-1) :: x3c
        real*8, dimension(1:n2-1) :: x2c
        real*8, dimension(1:n1-1) :: x1c

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: q1, q1S, q1S_expected, delta_q1S
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: q2, q2S, q2S_expected, delta_q2S
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: q3, q3S, q3S_expected, delta_q3S

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: Tx
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: Ty
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: Tz

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: RHS1, RHS2, RHS1_expected, delta_RHS1, RHS2_expected, delta_RHS2
        real*8  :: dT
        integer :: i, j, k, n1e, n3e

        ! x & xc **************************
        do i = 1, n3-1
            x3(i)=L3*(i-1.d0)/(n3-1.d0)
            x3c(i)=L3*(i-0.5d0)/(n3-1.d0)
        enddo
        x3(n3)=L3*(i-1.d0)/(n3-1.d0)

        do i = 1, n2-1
            x2(i)=L2*(i-1.d0)/(n2-1.d0)
            x2c(i)=L2*(i-0.5d0)/(n2-1.d0)
        enddo
        x2(n2)=L2*(i-1.d0)/(n2-1.d0)

        do i = 1, n1-1
            x1(i)=L1*(i-1.d0)/(n1-1.d0)
            x1c(i)=L1*(i-0.5d0)/(n1-1.d0)
        enddo
        x1(n1)=L1*(i-1.d0)/(n1-1.d0)

        ! q ******************************
        do k = xstart(3), min(xend(3), n3-1)
            do j = xstart(2), min(xend(2), n2-1)
                do i = xstart(1), xend(1)
                    q1(i,j,k)=dsin(x1(i))
                end do
            end do
        end do

        do k = ystart(3), min(yend(3), n3-1)
            do j = ystart(2), yend(2)
                do i = ystart(1), min(yend(1), n1-1)
                    q2(i,j,k)=dsin(pi*Y(j))
                end do
            end do
        end do

        do k = zstart(3), zend(3)
            do j = zstart(2), min(zend(2), n2-1)
                do i = zstart(1), min(zend(1), n1-1)
                    q3(i,j,k)=dsin(x3(k))
                end do
            end do
        end do

        ! T ********************************
        do k = ystart(3), min(yend(3), n3-1)
            do j = ystart(2), min(yend(2), n2-1)
                do i = ystart(1), min(yend(1), n1-1)
                    Ty(i,j,k)=dcos(x1c(i))*dsin(pi*Yc(j))*dcos(x3c(k))+(T0+(T1-T0)*Yc(j)/L2)
                end do
            end do
        end do

        call transpose_y_to_x(Ty, Tx)
        call transpose_y_to_z(Ty, Tz)

        ! qS_expected***********************
        do k = xstart(3), min(xend(3), n3-1)
            do j = xstart(2), min(xend(2), n2-1)
                do i = xstart(1), xend(1)
                    q1S_expected(i,j,k)=dsin(x1c(i))*Tx(i,j,k)
                end do
            end do
        end do

        do k = ystart(3), min(yend(3), n3-1)
            do j = ystart(2), yend(2)
                do i = ystart(1), min(yend(1), n1-1)
                    q2S_expected(i,j,k)=dsin(pi*Yc(j))*Ty(i,j,k)
                end do
            end do
        end do

        do k = zstart(3), zend(3)
            do j = zstart(2), min(zend(2), n2-1)
                do i = zstart(1), min(zend(1), n1-1)
                    q3S_expected(i,j,k)=dsin(x3c(k))*Tz(i,j,k)
                end do
            end do
        end do

        ! Snapshots of q, Ty and deltas
        call create_snapshot(COMMON_snapshot_path, "DEBUG", q1, "q1", 1)
        call create_snapshot(COMMON_snapshot_path, "DEBUG", q2, "q2", 2)
        call create_snapshot(COMMON_snapshot_path, "DEBUG", q3, "q3", 3)
        call create_snapshot(COMMON_snapshot_path, "DEBUG", Ty, "Ty", 2)


        !call perform_velocityscalar_products(q1, q2, q3, Tx, Ty, Tz, q1S, q2S, q3S)

        ! Checking of qS
        delta_q1S=q1S_expected-q1S
        delta_q2S=q2S
        delta_q3S=q3S_expected-q3S
        call create_snapshot(COMMON_snapshot_path, "DEBUG", delta_q1S, "delta_q1S", 1)
        call create_snapshot(COMMON_snapshot_path, "DEBUG", delta_q2S, "delta_q2S", 2)
        call create_snapshot(COMMON_snapshot_path, "DEBUG", delta_q3S, "delta_q3S", 3)

        !call perform_RHS1(q1S, q2S, q3S, RHS1)
        call perform_RHS2(RHS2)

        RHS1_expected=0.d0
        RHS2_expected=0.d0
        ! Expected RHS1 X-Direction
        do k = ystart(3), min(yend(3), n3-1)
            do j = ystart(2), min(yend(2), n2-1)
                do i = ystart(1), min(yend(1), n1-1)
                    RHS1_expected(i,j,k)=-(dcos(x1c(i))*Ty(i,j,k)-dsin(x1c(i))*dsin(x1c(i))*dsin(pi*Yc(j))*dcos(x3c(k)))
                    RHS2_expected(i,j,k)=-dcos(x1c(i))*dsin(pi*Yc(j))*dcos(x3c(k))/ren*prandtl
                end do
            end do
        end do

        ! Expected RHS1 Y-Direction
        do k = ystart(3), min(yend(3), n3-1)
            do j = ystart(2), min(yend(2), n2-1)
                do i = ystart(1), min(yend(1), n1-1)
                    dT=pi*dcos(pi*Yc(j))*dcos(x1c(i))*dcos(x3c(k))+1.d0/L2
                    RHS1_expected(i,j,k)=RHS1_expected(i,j,k)-(pi*dcos(pi*Yc(j))*Ty(i,j,k)+dT*dsin(pi*Yc(j)))
                    RHS2_expected(i,j,k)=RHS2_expected(i,j,k)-((1.d0/ren*prandtl)*(pi**2)*dsin(pi*Yc(j))*dcos(x1c(i))*dcos(x3c(k)))
                end do
            end do
        end do

        ! Expected RHS1 Z-Direction
        do k = ystart(3), min(yend(3), n3-1)
            do j = ystart(2), min(yend(2), n2-1)
                do i = ystart(1), min(yend(1), n1-1)
                    RHS1_expected(i,j,k)=RHS1_expected(i,j,k)-(dcos(x3c(k))*Ty(i,j,k)-dsin(x3c(k))*dsin(x3c(k))*dcos(x1c(i))*dsin(pi*Yc(j)))
                    RHS2_expected(i,j,k)=RHS2_expected(i,j,k)-((1.d0/ren*prandtl)*dcos(x3c(k))*dcos(x1c(i))*dsin(pi*Yc(j)))
                end do
            end do
        end do

        ! Checking of RHS1
        delta_RHS1=RHS1-RHS1_expected
        delta_RHS2=RHS2-RHS2_expected
        call create_snapshot(COMMON_snapshot_path, "DEBUG", delta_RHS1, "delta_RHS1", 2)
        call create_snapshot(COMMON_snapshot_path, "DEBUG", delta_RHS2, "delta_RHS2", 2)
        call create_snapshot(COMMON_snapshot_path, "DEBUG", RHS2_expected, "RHS2_exp", 2)
        call create_snapshot(COMMON_snapshot_path, "DEBUG", RHS2, "RHS2", 2)

    end subroutine check_solver


    subroutine save_state()

        use SCALAR_workspace_view, only: recovery_RHS_dir
        use HDF5_IO

        implicit none

        character(200)    :: file_path

        file_path=trim(recovery_RHS_dir)//"/RHS1"
        call hdf_write_3Dfield(file_path, previousRHS1(:,:,:), "RHS1", nx_global, ny_global, nz_global, ystart(1),yend(1),ystart(2),yend(2),ystart(3),yend(3))

        file_path=trim(recovery_RHS_dir)//"/RHS2"
        call hdf_write_3Dfield(file_path, previousRHS2(:,:,:), "RHS2", nx_global, ny_global, nz_global, ystart(1),yend(1),ystart(2),yend(2),ystart(3),yend(3))

    end subroutine save_state
!
!
end module SCALAR_solver
