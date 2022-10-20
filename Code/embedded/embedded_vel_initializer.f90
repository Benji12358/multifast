module embedded_turbulence_generator

    use embedded_physical_fields
    use embedded_snapshot_writer
    use DNS_settings
    use embedded_mesh

    implicit none

contains

    ! *****************************************************************************************************************************************
    ! *****************************************************************************************************************************************

    ! "Generic" generator of a turbulent field. The generated turbulent field depends on a various set of parameters
    ! (boundary condition, inflow type etc.):

    ! This subroutine can generate the following initial flows :
    !   - Poiseuille profile in y direction
    !   - 2D Poiseuille profile in x-y or y-z directions
    !   - Building of an inner flow from an input flow in x1 direction for open case (when no 3D starting file is given by user)
    !   - A basic jet support in x1 direction (due to the test kind of this feature, the jet parameters are defined directly in subroutine)

    !   Your are free to implement your own subroutine to study other cases (Coutte flow for example)

    ! *****************************************************************************************************************************************
    ! *****************************************************************************************************************************************

    subroutine generate_meanprofile(disturbance_intensity) !! y dummy argument

        use embedded_irregular_derivative_coefficients
        use mathematical_constants
        use DNS_settings
        use boundaries


        use mpi
        use decomp_2d

        implicit none

        real*8, intent(in)  :: disturbance_intensity
        real*8              :: modulation_j

        integer j,k,i,l, mpi_err
        integer :: n1s, n1e, n2s, n2e, n3s, n3e

        real*8, dimension(ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)) :: ts1
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2)) :: ts3
        real*8  :: glob_value

        real*8  :: yr, ycenter
        integer :: k1, k2, j1, j2

        if (nrank==0) write(*,*)"Use of the generic field generator"


        ! **************************************************************************************
        ! Velocity field initialisation --------------------------------------------------------
        ! **************************************************************************************

        ts1=0.d0
        ts3=0.d0

        ! Generation of an uniform flow in xy or yz planes
        if (flow_type==CONSTANT_FLOW) then
            if (streamwise==1)  ts1=u_bulk
        end if

        if(nrank==0) write(*,*)'streamwise =',streamwise

        ! Generation of 1D or 2D Poiseuille profile
        if ((flow_type==CHANNEL_FLOW).or.(flow_type==VORTICES)) then
            if (streamwise==1) call perform_stream1(ts1, BC2, BC3)
        end if

        ! Generation of inner flow from the first input flow
        ! Useful for open flows
        ! if (flow_type==FLOW_FROM_INFLOW) then
        !     call get_inflow(ts1)
        ! end if

        if (flow_type==BOUNDARY_LAYER) then
            if (streamwise==1) call perform_boundary_layer_1(ts1, delta_BL)
        end if

        ! At this point, the 2D array ts1 and ts3 contain the correct q1 and q3 profiles
        ! and can be used to define the inner velocity field.
        do j=1,n2
            do k=ystart_e(3), yend_e(3)
                do i=ystart_e(1), yend_e(1)

                    q3_y(i,j,k)= ts3(i,j)
                    q2_y(i,j,k)= 0.d0
                    q1_y(i,j,k)= ts1(j,k)        ! Poiseuille

                enddo
            enddo
        enddo

        if (yend_e(3)==n3) then
            do j = 1, n2
                do i = ystart_e(1), yend_e(1)
                    q1_y(i,j,n3)=0.d0
                    q2_y(i,j,n3)=0.d0
                    q3_y(i,j,n3)=0.d0
                end do
            end do
        end if

        do k = ystart_e(3), yend_e(3)
            do i = ystart_e(1), yend_e(1)
                q1_y(i,n2,k)=0.d0
                q2_y(i,n2,k)=0.d0
                q3_y(i,n2,k)=0.d0
            end do
        end do

        call MPI_BARRIER(MPI_COMM_WORLD , mpi_err)

        !Impermeability conditions ********************************************************

        ! X1 direction
        if ((BC1==FREESLIP).or.(BC1==NOSLIP)) then

            q3_wall10(xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)) =0.d0
            q2_wall10(xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)) =0.d0
            q1_wall10(xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)) =0.d0
            q1_x(1, xstart_e(2):xend_e(2), xstart_e(3):xend_e(3))   =0.d0

            q3_wall11(xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)) =0.d0
            q2_wall11(xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)) =0.d0
            q1_wall11(xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)) =0.d0
            q1_x(n1, xstart_e(2):xend_e(2), xstart_e(3):xend_e(3))  =0.d0

        end if

        ! X2 direction
        if ((BC2==FREESLIP).or.(BC2==NOSLIP)) then

            q3_wall20(ystart_e(1):yend_e(1), ystart_e(3):yend_e(3)) =0.d0
            q2_wall20(ystart_e(1):yend_e(1), ystart_e(3):yend_e(3)) =0.d0
            q1_wall20(ystart_e(1):yend_e(1), ystart_e(3):yend_e(3)) =0.d0
            q2_y(ystart_e(1):yend_e(1), 1, ystart_e(3):yend_e(3))   =0.d0

            q3_wall21(ystart_e(1):yend_e(1), ystart_e(3):yend_e(3)) =0.d0
            q2_wall21(ystart_e(1):yend_e(1), ystart_e(3):yend_e(3)) =0.d0
            q1_wall21(ystart_e(1):yend_e(1), ystart_e(3):yend_e(3)) =0.d0
            q2_y(ystart_e(1):yend_e(1), n2, ystart_e(3):yend_e(3))  =0.d0

        end if

        ! Warning : not checked
        ! X3 direction
        if ((BC3==FREESLIP).or.(BC3==NOSLIP)) then

            q3_wall30(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2))=0.d0
            q2_wall30(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2))=0.d0
            q1_wall30(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2))=0.d0
            q3_z(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), 1) =0.d0

            q3_wall31(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2))=0.d0
            q2_wall31(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2))=0.d0
            q1_wall31(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2))=0.d0
            q3_z(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), n3) =0.d0

        end if

        return


    contains
        subroutine perform_stream1(stream1, BC2, BC3)
            implicit none

            real*8, dimension(ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)) :: stream1
            integer                                                 :: BC2, BC3

            real*8                                                  :: f2(n2), f3(n3), c2, c3
            integer                                                 :: j, k

            f3=1.d0

            if (BC2==NOSLIP) then


                if (mod(n2, 2)==0) then
                    c2=Yc(n2/2)
                else
                    c2=Y((n2-1)/2+1)
                end if

                do j=1,n2-1
                    f2(j)=(1.d0-((Yc(j)-c2)/(0.5d0*L2))**2)
                enddo
                f2(n2)=f2(n2-1)

                do k = ystart_e(3),yend_e(3)
                    do j= ystart_e(2),yend_e(2)
                        stream1(j,k)=f2(j)*f3(k)
                    end do
                end do

            else
                stream1=0.d0
            end if

        end subroutine perform_stream1


        subroutine perform_boundary_layer_1(stream1, delta_BL)
            implicit none

            real*8, dimension(ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)) :: stream1

            real*8                                                  :: delta_BL
            integer                                                 :: j
            logical                                                 :: pair_n2

            pair_n2 = (mod(n2, 2)==0)

            do j=ystart_e(2),n2/2
                if (Yc(j)<delta_BL) then
                    !stream1(j,:) = 1.5d0 * (Yc(j)/delta_BL) - 0.5d0 * (Yc(j)/delta_BL)**2
                    !stream1(n2-j,:) = 1.5d0 * (Yc(j)/delta_BL) - 0.5d0 * (Yc(j)/delta_BL)**2
                    stream1(j,:) = 2.d0 * (Yc(j)/delta_BL) - 2.d0 * (Yc(j)/delta_BL)**3 + 1.d0 * (Yc(j)/delta_BL)**4
                    stream1(n2-j,:) = 2.d0 * (Yc(j)/delta_BL) - 2.d0 * (Yc(j)/delta_BL)**3 + 1.d0 * (Yc(j)/delta_BL)**4
                else
                    stream1(j,:) = 1.d0
                    stream1(n2-j,:) = 1.d0
                endif
            enddo

            if (pair_n2) then
                stream1(n2/2+1,:) = 1.d0
            endif

        end subroutine perform_boundary_layer_1

    end subroutine

    subroutine add_noise_rand_sin_new(disturbance_intensity)
        use embedded_mesh
        use embedded_irregular_derivative_coefficients
        use mathematical_constants
        use DNS_settings
        use boundaries
        use embedded_fringe_data


        use mpi
        use decomp_2d

        implicit none

        real*8, intent(in)  :: disturbance_intensity
        real*8              :: disturbance_intensity_at_j


        real*8, dimension (n2)    :: cac

        real*8, dimension(n2) :: v30
        real*8, dimension(n2)   :: v1m, v2m, v3m, q1m, q2m
        integer j,k,i,l, mpi_err
        real*8 vol,velocity_poisseuille,velocity_Z_disturbance,velocity_Y_disturbance,velocity_X_disturbance
        integer, parameter  :: ndv=3
        real*8 vmax(ndv)

        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)) :: rand_u, rand_v, rand_w
        real*8  :: glob_value, ph1, ph2, ph3

        integer :: writer_proc=0
        integer :: log_file_id=20
        character(200)  :: field_generator_path
        integer :: meanXZ_perturbation_file_id=21, zone_id, var_id
        integer, dimension(n2m) :: j_array
        integer             :: n1s, n2s, n3s, n1e, n2e, n3e

        real*8  :: ycenter

        call random_number(rand_u)
        call random_number(rand_v)
        call random_number(rand_w)
        ! **************************************************************************************
        ! Pressure field initialisation --------------------------------------------------------
        ! **************************************************************************************

!        pr_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))=0.d0



        ! **************************************************************************************
        ! Velocity field initialisation --------------------------------------------------------
        ! **************************************************************************************
        vol=dfloat(n1m*n2m*n3m)*2.d0


        cac(1)  =(y(2)-y(1))*dx2
        do j=2,n2m
            cac(j)=(y(j+1)-y(j-1))*dx2*0.5d0
        enddo
        cac(n2) =(y(n2)-y(n2m))*dx2


        !  v30(j): Poiseuille profile at the center of the cell
        !  v1m(j): Mean following X, Z of the Z component of the random perturbation
        !  v2m(j): Mean following X, Z of the Y component of the random perturbation
        !  v3m(j): Mean following X, Z of the X component of the random perturbation

        !        TO NOTE : The random number generator depends on the computer

        n1s = max(1, ystart_e(1))
        n3s = max(1, ystart_e(3))

        if (streamwise==3) n3s=max(2, ystart_e(3))
        if (streamwise==1) n1s=max(1, ystart_e(1))
        n2s = 1

        n1e = min(n1m, yend_e(1))
        n2e = n2m
        n3e = min(n3m, yend_e(3))

        do j=n2s,n2e
            v1m(j)=0.d0
            v2m(j)=0.d0
            v3m(j)=0.d0

            if (mod(n2, 2)==0) then
                ycenter=Yc(n2/2)
            else
                ycenter=Y((n2-1)/2+1)
            end if

            ! POISEUILLE PROFILE v30 !
            v30(j)=(1.d0-((Yc(j)-ycenter)/(0.5d0*L2))**2)

            disturbance_intensity_at_j=disturbance_intensity
            if(dabs(1.d0-Yc(j)).lt.0.025) disturbance_intensity_at_j=disturbance_intensity/5.d0

            do k=n3s,n3e
                do i=n1s,n1e

                    ph1=4.d0*2.d0*pi*(i-1)/(n1-1)
                    ph2=6.d0*2.d0*pi*(j-1)/(n2-1)
                    ph3=4.d0*2.d0*pi*(k-1)/(n3-1)

                    !   random disturbance  on u, v and w
                    q1_y(i,j,k)=disturbance_intensity_at_j*(-1.d0+2.d0*rand_w(i,j,k))
                    q2_y(i,j,k)=disturbance_intensity_at_j*(-1.d0+2.d0*rand_v(i,j,k))
                    q3_y(i,j,k)=disturbance_intensity_at_j*(-1.d0+2.d0*rand_u(i,j,k))

                    !  adding sinusoidal disturbance on u, v, w
                    q1_y(i,j,k)=q1_y(i,j,k) + disturbance_intensity*dsin(ph1)*dsin(ph2)*dsin(ph3)
                    q2_y(i,j,k)=q2_y(i,j,k) + disturbance_intensity*(dsin(ph1)*dsin(ph2)*dsin(ph3))**2.d0
                    q3_y(i,j,k)=q3_y(i,j,k) + disturbance_intensity*dsin(ph1)*dsin(ph2)*dsin(ph3)

                    v1m(j)=v1m(j)+q1_y(i,j,k)
                    v2m(j)=v2m(j)+q2_y(i,j,k)
                    v3m(j)=v3m(j)+q3_y(i,j,k)


                enddo
            enddo

            call MPI_ALLREDUCE (v1m(j), glob_value, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
            v1m(j)=glob_value/dfloat(n1m*n3m)

            call MPI_ALLREDUCE (v2m(j), glob_value, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
            v2m(j)=glob_value/dfloat(n1m*n3m)

            call MPI_ALLREDUCE (v3m(j), glob_value, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
            v3m(j)=glob_value/dfloat(n1m*n3m)

        enddo

        velocity_poisseuille    =0.d0
        velocity_Z_disturbance  =0.d0
        velocity_Y_disturbance  =0.d0
        velocity_X_disturbance  =0.d0

        do j=n2s, n2e
            do k=n3s, n3e
                do i=n1s, n1e
!
                    !! Setting the mean along (X, Z) directions to 0
                    q3_y(i,j,k)=q3_y(i,j,k) - v3m(j)
                    q2_y(i,j,k)=q2_y(i,j,k) - v2m(j)
                    q1_y(i,j,k)=q1_y(i,j,k) - v1m(j)

                    ! adding the poiseuille profile for streamwise velocity
                    if (streamwise==3)  q3_y(i,j,k)=q3_y(i,j,k) + v30(j)
                    if (streamwise==1)  q1_y(i,j,k)=q1_y(i,j,k) + v30(j)


                    velocity_poisseuille    = velocity_poisseuille      +   v30(j)*cell_size_Y(j)/dx2
                    velocity_Z_disturbance  = velocity_Z_disturbance    +   q1_y(i,j,k)*cell_size_Y(j)/dx2
                    velocity_Y_disturbance  = velocity_Y_disturbance    +   (q2_y(i,j,k)*cac(j) + q2_y(i,j+1,k)*cac(j+1) )*0.5d0
                    velocity_X_disturbance  = velocity_X_disturbance    +   q3_y(i,j,k)*cell_size_Y(j)/dx2

                enddo
            enddo
        enddo


        call MPI_ALLREDUCE (velocity_poisseuille, glob_value, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        velocity_poisseuille=glob_value/vol

        call MPI_ALLREDUCE (velocity_Z_disturbance, glob_value, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        velocity_Z_disturbance=glob_value/vol

        call MPI_ALLREDUCE (velocity_Y_disturbance, glob_value, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        velocity_Y_disturbance=glob_value/vol

        call MPI_ALLREDUCE (velocity_X_disturbance, glob_value, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        velocity_X_disturbance=glob_value/vol


        !Impermeability conditions
        do k=ystart_e(3), min(n3m, yend_e(3))
            do i=ystart_e(1), min(n1m, yend_e(1))
                q2_y(i,1,k)=0.d0
                q2_y(i,n2,k)=0.d0
            enddo
        enddo

        vmax(1)=MAXVAL(abs(q1_y))
        vmax(2)=MAXVAL(abs(q2_y))
        vmax(3)=MAXVAL(abs(q3_y))


        call MPI_ALLREDUCE (vmax(1), glob_value, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
        vmax(1)=glob_value


        call MPI_ALLREDUCE (vmax(2), glob_value, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
        vmax(2)=glob_value


        call MPI_ALLREDUCE (vmax(3), glob_value, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
        vmax(3)=glob_value

        call transpose_y_to_z(q3_y, q3_z, decomp_embedded)
        call transpose_y_to_x(q1_y, q1_x, decomp_embedded)


!         ! **************************************************************************************
!         ! Pressure field initialisation --------------------------------------------------------
!         ! **************************************************************************************
! 115     format(' Bulk velocity of the Poiseuille profile velocity_poisseuille=',e10.3)
! 116     format(' velocity_X_disturbance=',e10.3,2x,'velocity_Y_disturbance=',e10.3,2x,'velocity_Z_disturbance=',e10.3)
! 117     format(' Maximum Velocity w max=',e11.4,2x, 'v max=',e11.4,2x,'u max=',e11.4)

!         if (nrank==writer_proc) then

!             field_generator_path=trim(log_path)//"Fields_generator/"

!             open(log_file_id,file=trim(field_generator_path)//'overview.out')

!             write(log_file_id,115) velocity_poisseuille
!             write(log_file_id,116) velocity_X_disturbance,velocity_Y_disturbance,velocity_Z_disturbance
!             write(log_file_id,117) (vmax(l),l=1,3)
!             write(log_file_id,*)'Mean along x1 & x3 of the initial perturbation'

!             write(log_file_id,*)'j, q1m(j), q2m(j), v3m(j)'
!             do j=1,n2m
!                 write(log_file_id,615)j,q1m(j),q2m(j),v3m(j)
!             enddo

! 615         format(2x,i3,2x,3e12.4)
!             close(log_file_id)

!         end if

        call MPI_BARRIER(MPI_COMM_WORLD , mpi_err)

        return

!
    end subroutine add_noise_rand_sin_new


    subroutine init_turbulent_field(intensity)
        use COMMON_workspace_view, only: COMMON_snapshot_path
        implicit none
        real*8, intent(in)      :: intensity

        ! Build the mean profile
        call generate_meanprofile(intensity)

        ! Add noise
        if(intensity/=0.d0) call add_noise_rand_sin_new(intensity)   ! to do/check

        ! Transpose data in 2D-decomposition stencil

        call transpose_y_to_x(q3_y, q3_x, decomp_embedded)
        call transpose_y_to_z(q3_y, q3_z, decomp_embedded)

        call transpose_y_to_x(q1_y, q1_x, decomp_embedded)
        call transpose_y_to_z(q1_y, q1_z, decomp_embedded)

        call transpose_y_to_x(q2_y, q2_x, decomp_embedded)
        call transpose_y_to_z(q2_y, q2_z, decomp_embedded)

        ! The final 3D field are exported for checking purposes
        call create_snapshot(COMMON_snapshot_path, "EMBEDDED_INIT", q1_y, "W", 2)
        call create_snapshot(COMMON_snapshot_path, "EMBEDDED_INIT", q2_y, "V", 2)
        call create_snapshot(COMMON_snapshot_path, "EMBEDDED_INIT", q3_y, "U", 2)

        !write(*,*) 'q3_x', q3_x(32,:,n3m)


    end subroutine init_turbulent_field


end module embedded_turbulence_generator
