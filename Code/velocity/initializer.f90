module Turbulence_generator

    use physical_fields
    use snapshot_writer
    use DNS_settings
    use mesh

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

        use irregular_derivative_coefficients
        use mathematical_constants
        use DNS_settings
        use inflow_settings
        use boundaries


        use mpi
        use decomp_2d

        implicit none

        real*8, intent(in)  :: disturbance_intensity
        real*8              :: modulation_j

        integer j,k,i,l, mpi_err
        integer :: n1s, n1e, n2s, n2e, n3s, n3e

        real*8, dimension(ystart(2):yend(2), ystart(3):yend(3)) :: ts1
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2)) :: ts3
        real*8  :: glob_value

        real*8  :: yr, ycenter
        integer :: k1, k2, j1, j2

        if (nrank==0) write(*,*)"Use of the generic field generator"

        ! TEST of a jet flow: jet parameter are hard coded here
        if (.true.) then
            k1=22
            k2=42
            j1=22
            j2=42
        endif
        if(.false.) then
            k1=1!22
            k2=n3!42
            j1=1!22
            j2=n2!42
        endif

!        pr_x=0.d0


        ! **************************************************************************************
        ! Velocity field initialisation --------------------------------------------------------
        ! **************************************************************************************


        ts1=0.d0
        ts3=0.d0

        ! Generation of an uniform flow in xy or yz planes
        if (flow_type==CONSTANT_FLOW) then
            if (streamwise==1)  ts1=1.d0
            if (streamwise==3)  ts3=1.d0
        end if

        if(nrank==0) write(*,*)'streamwise =',streamwise

        ! Generation of 1D or 2D Poiseuille profile
        if ((flow_type==CHANNEL_FLOW).or.(flow_type==VORTICES)) then
            if (streamwise==1) call perform_stream1(ts1, BC2, BC3)
            if (streamwise==3) call perform_stream3(ts3, BC1, BC2)
        end if

        ! Generation of inner flow from the first input flow
        ! Useful for open flows
        if (flow_type==FLOW_FROM_INFLOW) then
            call get_inflow(ts1)
        end if



        ! At this point, the 2D array ts1 and ts3 contain the correct q1 and q3 profiles
        ! and can be used to define the inner velocity field.
        do j=1,n2
            do k=ystart(3), yend(3)
                do i=ystart(1), yend(1)

                    q3_y(i,j,k)= ts3(i,j)
                    q2_y(i,j,k)= 0.d0
                    q1_y(i,j,k)= ts1(j,k)        ! Poiseuille

                enddo
            enddo
        enddo

        ! Generation of vortices
        if (flow_type==VORTICES) then
            call perform_vortices
        end if

        ! CLEAN FIELDS ***************************************************************
        if (yend(1)==n1) then
            do k = ystart(3), yend(3)
                do j = 1, n2
                    q1_y(n1,j,k)=0.d0
                    q2_y(n1,j,k)=0.d0
                    q3_y(n1,j,k)=0.d0
                    !if ((BC1==OPEN).or.(BC1==UNBOUNDED).or.(BC1==UNBOUNDED)
                end do
            end do
        end if

        if (yend(3)==n3) then
            do j = 1, n2
                do i = ystart(1), yend(1)
                    q1_y(i,j,n3)=0.d0
                    q2_y(i,j,n3)=0.d0
                    q3_y(i,j,n3)=0.d0
                end do
            end do
        end if

        do k = ystart(3), yend(3)
            do i = ystart(1), yend(1)
                q1_y(i,n2,k)=0.d0
                q2_y(i,n2,k)=0.d0
                q3_y(i,n2,k)=0.d0
            end do
        end do

        call MPI_BARRIER(MPI_COMM_WORLD , mpi_err)

        !Impermeability conditions ********************************************************

        ! X1 direction
        if ((BC1==FREESLIP).or.(BC1==NOSLIP)) then

            q3_wall10(xstart(2):xend(2), xstart(3):xend(3)) =0.d0
            q2_wall10(xstart(2):xend(2), xstart(3):xend(3)) =0.d0
            q1_wall10(xstart(2):xend(2), xstart(3):xend(3)) =0.d0
            q1_x(1, xstart(2):xend(2), xstart(3):xend(3))   =0.d0

            q3_wall11(xstart(2):xend(2), xstart(3):xend(3)) =0.d0
            q2_wall11(xstart(2):xend(2), xstart(3):xend(3)) =0.d0
            q1_wall11(xstart(2):xend(2), xstart(3):xend(3)) =0.d0
            q1_x(n1, xstart(2):xend(2), xstart(3):xend(3))  =0.d0

        end if

        ! X2 direction
        if ((BC2==FREESLIP).or.(BC2==NOSLIP)) then

            q3_wall20(ystart(1):yend(1), ystart(3):yend(3)) =0.d0
            q2_wall20(ystart(1):yend(1), ystart(3):yend(3)) =0.d0
            q1_wall20(ystart(1):yend(1), ystart(3):yend(3)) =0.d0
            q2_y(ystart(1):yend(1), 1, ystart(3):yend(3))   =0.d0

            q3_wall21(ystart(1):yend(1), ystart(3):yend(3)) =0.d0
            q2_wall21(ystart(1):yend(1), ystart(3):yend(3)) =0.d0
            q1_wall21(ystart(1):yend(1), ystart(3):yend(3)) =0.d0
            q2_y(ystart(1):yend(1), n2, ystart(3):yend(3))  =0.d0

        end if

        ! Warning : not checked
        ! X3 direction
        if ((BC3==FREESLIP).or.(BC3==NOSLIP)) then

            q3_wall30(zstart(1):zend(1), zstart(2):zend(2))=0.d0
            q2_wall30(zstart(1):zend(1), zstart(2):zend(2))=0.d0
            q1_wall30(zstart(1):zend(1), zstart(2):zend(2))=0.d0
            q3_z(zstart(1):zend(1), zstart(2):zend(2), 1) =0.d0

            q3_wall31(zstart(1):zend(1), zstart(2):zend(2))=0.d0
            q2_wall31(zstart(1):zend(1), zstart(2):zend(2))=0.d0
            q1_wall31(zstart(1):zend(1), zstart(2):zend(2))=0.d0
            q3_z(zstart(1):zend(1), zstart(2):zend(2), n3) =0.d0

        end if

        ! ATTENTION
        if (BC1==OPEN) then

            call random_number(q2_wall10(xstart(2):xend(2), xstart(3):xend(3)))
            call random_number(q3_wall10(xstart(2):xend(2), xstart(3):xend(3)))

            do k=xstart(3), xend(3)
                do j=xstart(2), min(xend(2),n2-1)

                    yr=abs(Yc(j)-L2/2.d0)
                    modulation_j=disturbance_intensity*dexp(-0.2d0*yr**2)

                    if (inflow_mode==INFLOW_SQUARE) then

                        if ((k.ge.inflow_sq_k1).and.(k.lt.inflow_sq_k2).and.(j.ge.inflow_sq_j1).and.(j.lt.inflow_sq_j2)) then
                            q1_wall10(j, k)=inflow_int
                            q1_x(1,j,k)=q1_wall10(j, k)
                            q2_wall10(j, k)=0.d0
                            q3_wall10(j, k)=0.d0
                        else
                            q1_wall10(j, k)=q1_x(1,j,k)
                            q2_wall10(j, k)=q2_wall10(j, k)*modulation_j
                            q3_wall10(j, k)=q3_wall10(j, k)*modulation_j
                        endif

                    endif

                    if (inflow_mode==INFLOW_NONE) then

                        q1_wall10(j, k)=q1_x(1,j,k)
                        q2_wall10(j, k)=q2_wall10(j, k)*modulation_j
                        q3_wall10(j, k)=q3_wall10(j, k)*modulation_j

                    endif

                    q3_wall11(j, k) = q3_x(n1-1, j, k)
                    q2_wall11(j, k) = q2_x(n1-1, j, k)
                    q1_wall11(j, k) = q1_x(n1-1, j, k)
                    q1_x(n1, j, k)  = q1_wall11(j, k)

                enddo
            enddo

        end if

        return


    contains
        subroutine perform_stream1(stream1, BC2, BC3)
            implicit none

            real*8, dimension(ystart(2):yend(2), ystart(3):yend(3)) :: stream1
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

                do k = ystart(3),yend(3)
                    do j= ystart(2),yend(2)
                        stream1(j,k)=f2(j)*f3(k)
                    end do
                end do

            else
                stream1=0.d0
            end if

        end subroutine perform_stream1


        subroutine perform_stream3(stream3, BC1, BC2)
            implicit none

            real*8, dimension(ystart(1):yend(1), ystart(2):yend(2)) :: stream3
            integer                                                 :: BC1, BC2

            real*8                                                  :: f1(n1), f2(n2), c1, c2
            integer                                                 :: i, j,temp_test

            f1=1.d0


            if (BC2==NOSLIP) then


                do j=1,n2

                    if (mod(n2, 2)==0) then
                       c2=Yc(n2/2)
                    else
                        c2=Y((n2-1)/2+1)
                    end if

                    f2(j)=(1.d0-((Y(j)-c2)/(0.5d0*L2))**2)

                enddo

               if (BC1==NOSLIP) then
               !     do j = ystart(2),yend(2)
               !      do i= ystart(1),yend(1)
               !         stream3(i,j)=f2(j)
               !      end do
               !    end do

                    if (mod(n1, 2)==0) then
                       c1=(n1/2-0.5d0)*dx1
                    else
                      c1=((n1-1)/2+1-0.5d0)*dx1
                    end if


                    call rectangle_channel_flow_init(stream3,c1,c2)
                    
                 else !

                    do j = ystart(2),yend(2)
                     do i= ystart(1),yend(1)
                        stream3(i,j)=f2(j)
                     end do
                   end do

               endif !BC1

            else !
                stream3=0.d0
            end if !BC2

        end subroutine perform_stream3


        subroutine rectangle_channel_flow_init(stream3_temp,c1,c2)
        use mathematical_constants
        implicit none
        real*8, intent(in)      :: c1,c2
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2)) :: stream3_temp
        real*8  :: fact1,fact2,fact3,fact4,fact4p,term1,term2,term1p,term2p,Umax_lim
        real*8  :: k
        integer :: i,j,n,n_iter

        n_iter = 10

        do j = ystart(2),yend(2)
           do i= ystart(1),yend(1)
              stream3_temp(i,j) = (1.d0-((Y(j)-c2)/(0.5d0*L2))**2)
           enddo
        enddo


        do n = 1, n_iter
           k = real(2*n-1) ! odd numbers
           do j = ystart(2),yend(2)
             do i= ystart(1),yend(1)

             fact1 = k*PI*(Z(i)-c1)/(2.d0*(0.5d0*L2))
             fact2 = k*PI*(0.5d0*L1)/(2.d0*(0.5d0*L2))
             fact3 = k*PI*(Y(j)-c2)/(2.d0*(0.5d0*L2))
             fact4 = (-1.d0)**((k-1.d0)/2.d0)
             term1 = dcos(fact3)/(PI*k/2.d0)**3.d0
             term2 = (( dcosh(fact1)/dcosh(fact2) ))

             stream3_temp(i,j) = stream3_temp(i,j) - 4.d0*fact4*term2*term1
             enddo
           enddo
        enddo
!
!         write(*,*)'stream3_temp(32,32)=', stream3_temp(32,32)
!         write(*,*)'stream3_temp(1,1)=', stream3_temp(1,1)
!
        Umax_lim = 1.d0


        do n = 1,100

           k = real(2*n-1) !odd numbers only

           fact4p = (-1.d0)**((k-1.d0)/2.d0)

           term2p = (-( 1.d0/dcosh(k*PI*(L1*0.5d0)/(2.d0*(0.5d0*L2))) ))

           term1p = 1.d0/((PI*k/2.d0)**3.d0)

           Umax_lim = Umax_lim + 4.d0*fact4p*term1p*term2p

        enddo

       write(6,*)'Umax_lim =', Umax_lim


       do j = ystart(2),yend(2)
         do i= ystart(1),yend(1)
            if (j==1 .or. j==n2m) then
                stream3_temp(i,j) = 0.d0
            elseif (i==1 .or. i==n1m) then
                stream3_temp(i,j) = 0.d0
            else
                stream3_temp(i,j) = stream3_temp(i,j)*(1.d0/Umax_lim)
            endif
         enddo
       enddo

    end subroutine


        subroutine get_inflow(ts11, ntime)
            use HDF5_IO

            use start_settings, only:start_it
            use COMMON_workspace_view, only: COMMON_inflow_path

            implicit none
            integer, optional   :: ntime
            integer, save   :: inflow_nb=1
            character(200)       :: current_inflow_path

            real*8, dimension(1, ystart(2):yend(2), ystart(3):yend(3))      :: ts11_tmp, ts12_tmp, ts13_tmp
            real*8, dimension(ystart(2):yend(2), ystart(3):yend(3))         :: ts11, ts12, ts13

            character*10 tmp_str

            write(tmp_str, "(i10)")inflow_nb+start_it-1
            current_inflow_path=trim(COMMON_inflow_path)//'outflow_'//trim(adjustl(tmp_str))

            if (nrank==0) write(*,*)
            if (nrank==0) write(*,*) "Reading inflow from file:", trim(current_inflow_path)//".h5"

            call hdf_read_3Dfield(current_inflow_path, ts11_tmp, "q1_out", 1, ny_global, nz_global, 1,1, ystart(2),yend(2), ystart(3),yend(3))
            !call hdf_read_3Dfield(current_inflow_path, ts12_tmp, "q2_out", 1, ny_global, nz_global, 1,1, xstart(2),xend(2), xstart(3),xend(3))
            !call hdf_read_3Dfield(current_inflow_path, ts13_tmp, "q3_out", 1, ny_global, nz_global, 1,1, xstart(2),xend(2), xstart(3),xend(3))

            ts11(:,:)=ts11_tmp(1, :,:)
            !ts12(:,:)=ts12_tmp(1, :,:)
            !ts13(:,:)=ts13_tmp(1, :,:)

        end subroutine get_inflow

        subroutine perform_vortices()
            use HDF5_IO
            use DNS_settings
            use VELOCITY_operations, only:perform_velocity_at_center
            implicit none

            real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))      :: wc_x
            real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))      :: uc_y, vc_y, wc_y
            real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))      :: uc_z

            real*8, dimension(2, ystart(2):yend(2), ystart(3):yend(3)) :: uv1
            real*8, dimension(2, ystart(1):yend(1), ystart(2):yend(2)) :: wv3
            integer     :: i,j,k
            integer     :: n1s, n1e, n2s, n2e, n3s, n3e
            real*8      :: x1trans, x2trans, x3trans, x1ctrans, x2ctrans, x3ctrans

            uv1=0.d0
            wv3=0.d0

            n1s=ystart(1)
            n1e=min(yend(1),n1m)
            n2s=1
            n2e=n2m
            n3s=ystart(3)
            n3e=min(yend(3),n3m)

            if (vort_dir==3)then

                do i= n1s,n1e
                    do j= n2s,n2e
                        x1ctrans= (Zc(i) - x1vs)*2.d0/dvort
                        x2ctrans= (Yc(j) - hvort)*2.d0/dvort
                        x1trans= (Z(i) - x1vs)*2.d0/dvort
                        x2trans= (Y(j) - hvort)*2.d0/dvort

                        if (sqrt(x1trans**2+x2ctrans**2)<=1.d0) then
                            wv3(1,i,j)= -Urot*x2ctrans
                        endif

                        if (sqrt(x1ctrans**2+x2trans**2)<=1.d0) then
                            wv3(2,i,j)=  Urot*x1ctrans
                        endif
                    enddo
                enddo

                do k = n3s,n3e
                    x3ctrans= (Xc(k) - x3vs)/Lvort

                    if ((x3ctrans>=0.d0).and.(x3ctrans<=1.d0)) then

                        do i= n1s,n1e
                            do j= n2s,n2e
                                q1_y(i,j,k) = q1_y(i,j,k) + wv3(1,i,j)
                                q2_y(i,j,k) = q2_y(i,j,k) + wv3(2,i,j)
                            enddo
                        enddo

                    endif

                enddo
            endif

            call transpose_y_to_x(q1_y,q1_x)

            call perform_velocity_at_center(q3_z, q2_y, q1_x, uc_z, vc_y, wc_x)
            call transpose_x_to_y(wc_x,wc_y)
            call transpose_z_to_y(uc_z,uc_y)


            if(nrank==0)  call hdf_create_file_with_3Dmesh('Test_vort', Zc, Yc,Xc, "x1", "x2", "x3", n1m, n2m,n3m)
            if(nrank==0)  call hdf_addgroup('Test_vort', "U")
            call hdf_add_3Dfield('Test_vort', wc_y(n1s:n1e,n2s:n2e,n3s:n3e), "U/U1", n1m, n2m, n3m, n1s, n1e, n2s, n2e, n3s, n3e)
            call hdf_add_3Dfield('Test_vort', vc_y(n1s:n1e,n2s:n2e,n3s:n3e), "U/U2", n1m, n2m, n3m, n1s, n1e, n2s, n2e, n3s, n3e)
            call hdf_add_3Dfield('Test_vort', uc_y(n1s:n1e,n2s:n2e,n3s:n3e), "U/U3", n1m, n2m, n3m, n1s, n1e, n2s, n2e, n3s, n3e)


        end subroutine perform_vortices

    end subroutine
!
!
! sin pertubation
    subroutine add_noise(intensity)
        use mathematical_constants
        implicit none
        real*8, intent(in)      :: intensity
        real*8  :: ph1, ph2, ph3
        integer :: i,j,k
!

        do j=1,n2
            do k=ystart(3), yend(3)
                do i=ystart(1), yend(1)
                    ph1=4.d0*2.d0*pi*(i-1)/(n1-1)
                    ph2=6.d0*2.d0*pi*(j-1)/(n2-1)
                    ph3=4.d0*2.d0*pi*(k-1)/(n3-1)
                    q1_y(i,j,k)=q1_y(i,j,k)+intensity*dsin(ph1)*dsin(ph2)*dsin(ph3)
                    q2_y(i,j,k)=intensity*(dsin(ph1)*dsin(ph2)*dsin(ph3))**2.d0
                    q3_y(i,j,k)=q3_y(i,j,k)+intensity*dsin(ph1)*dsin(ph2)*dsin(ph3)
                enddo
            enddo
        enddo
!
    end subroutine add_noise
!
!
! Adding a random disturbance to the mean flow
    subroutine add_noise_rand(disturbance_intensity)
        use mesh
        use irregular_derivative_coefficients
        use mathematical_constants
        use DNS_settings
        use workspace_view


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

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: rand_u, rand_v
        real*8  :: glob_value

        integer :: writer_proc=0
        integer :: log_file_id=20
        character(200)  :: field_generator_path
        integer :: meanXZ_perturbation_file_id=21, zone_id, var_id
        integer, dimension(n2m) :: j_array

        real*8  :: ycenter

        call random_number(rand_u)
        call random_number(rand_v)

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


        do j=1,n2m
            v1m(j)=0.d0
            v2m(j)=0.d0
            v3m(j)=0.d0

            disturbance_intensity_at_j=disturbance_intensity
            if((1.d0-dabs(Yc(j))).lt.0.025) disturbance_intensity_at_j=disturbance_intensity/5.d0

            do k=ystart(3), min(n3m, yend(3))
                do i=ystart(1), min(n1m, yend(1))

                    !   random disturbance  on u, v and w
                    q1_y(i,j,k)=disturbance_intensity_at_j
                    q2_y(i,j,k)=disturbance_intensity_at_j*(-1.d0+2.d0*rand_v(i,j,k))
                    q3_y(i,j,k)=disturbance_intensity_at_j*(-1.d0+2.d0*rand_u(i,j,k))

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


     do j=1,n2m

            if (mod(n2, 2)==0) then
                ycenter=Yc(n2/2)
            else
                ycenter=Y((n2-1)/2+1)
            end if

            v30(j)=(1.d0-((Yc(j)-ycenter)/(0.5d0*L2))**2)

            do k=ystart(3), min(n3m, yend(3))
                do i=ystart(1), min(n1m, yend(1))

                    q3_y(i,j,k)=q3_y(i,j,k)-v3m(j)
                    q3_y(i,j,k)=q3_y(i,j,k) + v30(j)        ! Poiseuille + disturbance

                    q2_y(i,j,k)=q2_y(i,j,k)-v2m(j)

                    q1_y(i,j,k)=q1_y(i,j,k)-v1m(j)


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
        do k=ystart(3), min(n3m, yend(3))
            do i=ystart(1), min(n1m, yend(1))
                q2_y(i,1,k)=0.d0
                q2_y(i,n2,k)=0.d0
            enddo
        enddo

        !   q1m(j): mean along (X, Z) of transverse disturbance    (=0)
        !   q2m(j): mean along (X, Z) of normal disturbance        (=0)
        do j=1,n2m
            q1m(j)=0.d0
            q2m(j)=0.d0

            do k=ystart(3), min(n3m, yend(3))
                do i=ystart(1), min(n1m, yend(1))
                    q1m(j)=q1m(j)+q1_y(i,j,k)
                    q2m(j)=q2m(j)+q2_y(i,j,k)
                enddo
            enddo

            call MPI_ALLREDUCE (q1m(j), glob_value, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
            q1m(j)=glob_value/dfloat(n1m*n3m)

            call MPI_ALLREDUCE (q2m(j), glob_value, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
            q2m(j)=glob_value/dfloat(n1m*n3m)
        enddo


        do j=2,n2m

            do k=ystart(3), min(n3m, yend(3))
                do i=ystart(1), min(n1m, yend(1))
                    q2_y(i,j,k)=q2_y(i,j,k)-q2m(j)
                enddo
            enddo

        enddo


        do j=1,n2m
            do k=ystart(3), min(n3m, yend(3))
                do i=ystart(1), min(n1m, yend(1))

                    q1_y(i,j,k)=q1_y(i,j,k)-q1m(j)

                enddo
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



        call transpose_y_to_z(q3_y, q3_z)
        call transpose_y_to_x(q1_y, q1_x)


        ! **************************************************************************************
        ! Pressure field initialisation --------------------------------------------------------
        ! **************************************************************************************
115     format(' Bulk velocity of the Poiseuille profile velocity_poisseuille=',e10.3)
116     format(' velocity_X_disturbance=',e10.3,2x,'velocity_Y_disturbance=',e10.3,2x,'velocity_Z_disturbance=',e10.3)
117     format(' Maximum Velocity w max=',e11.4,2x, 'v max=',e11.4,2x,'u max=',e11.4)

        if (nrank==writer_proc) then

            field_generator_path=trim(log_path)//"Fields_generator/"

            open(log_file_id,file=trim(field_generator_path)//'overview.out')

            write(log_file_id,115) velocity_poisseuille
            write(log_file_id,116) velocity_X_disturbance,velocity_Y_disturbance,velocity_Z_disturbance
            write(log_file_id,117) (vmax(l),l=1,3)
            write(log_file_id,*)'Mean along x1 & x3 of the initial perturbation'

            write(log_file_id,*)'j, q1m(j), q2m(j), v3m(j)'
            do j=1,n2m
                write(log_file_id,615)j,q1m(j),q2m(j),v3m(j)
            enddo

615         format(2x,i3,2x,3e12.4)
            close(log_file_id)

        end if

        call MPI_BARRIER(MPI_COMM_WORLD , mpi_err)

        return

!
    end subroutine add_noise_rand

! Adding a random disturbance to the mean flow and the sin pertubations (mixed strategy)
! Created OD on 13/02/2017
!
    subroutine add_noise_rand_sin(disturbance_intensity)
        use mesh
        use irregular_derivative_coefficients
        use mathematical_constants
        use DNS_settings
        use workspace_view


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

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: rand_u, rand_v
        real*8  :: glob_value, ph1, ph2, ph3

        integer :: writer_proc=0
        integer :: log_file_id=20
        character(200)  :: field_generator_path
        integer :: meanXZ_perturbation_file_id=21, zone_id, var_id
        integer, dimension(n2m) :: j_array

        real*8  :: ycenter

        call random_number(rand_u)
        call random_number(rand_v)

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


        do j=1,n2m
            v1m(j)=0.d0
            v2m(j)=0.d0
            v3m(j)=0.d0

            disturbance_intensity_at_j=disturbance_intensity
            if((1.d0-dabs(Yc(j))).lt.0.025) disturbance_intensity_at_j=disturbance_intensity/5.d0

            do k=ystart(3), min(n3m, yend(3))
                do i=ystart(1), min(n1m, yend(1))

                    !   random disturbance  on u, v and w
                    q1_y(i,j,k)=disturbance_intensity_at_j
                    q2_y(i,j,k)=disturbance_intensity_at_j*(-1.d0+2.d0*rand_v(i,j,k))
                    q3_y(i,j,k)=disturbance_intensity_at_j*(-1.d0+2.d0*rand_u(i,j,k))

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


     do j=1,n2m

            if (mod(n2, 2)==0) then
                ycenter=Yc(n2/2)
            else
                ycenter=Y((n2-1)/2+1)
            end if

            v30(j)=(1.d0-((Yc(j)-ycenter)/(0.5d0*L2))**2)

!           if (j.le.(n2m/2)) then
!           epsi = 0.01
!           else
!           epsi = 0.01
!           endif

            do k=ystart(3), min(n3m, yend(3))
                do i=ystart(1), min(n1m, yend(1))

                    ph1=4.d0*2.d0*pi*(i-1)/(n1-1)
                    ph2=6.d0*2.d0*pi*(j-1)/(n2-1)
                    ph3=4.d0*2.d0*pi*(k-1)/(n3-1)
!
! WARNING : OD 13/02/2017 : Streamwise = 3 ONLY for starting perturbation
!
                    q3_y(i,j,k)=q3_y(i,j,k) - v3m(j)
                    q3_y(i,j,k)=q3_y(i,j,k) + v30(j)
                    q3_y(i,j,k)=q3_y(i,j,k) + disturbance_intensity*dsin(ph1)*dsin(ph2)*dsin(ph3)

                    q2_y(i,j,k)=q2_y(i,j,k) - v2m(j) + disturbance_intensity*(dsin(ph1)*dsin(ph2)*dsin(ph3))**2.d0

                    q1_y(i,j,k)=q1_y(i,j,k) - v1m(j) + disturbance_intensity*dsin(ph1)*dsin(ph2)*dsin(ph3)


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
        do k=ystart(3), min(n3m, yend(3))
            do i=ystart(1), min(n1m, yend(1))
                q2_y(i,1,k)=0.d0
                q2_y(i,n2,k)=0.d0
            enddo
        enddo

        !   q1m(j): mean along (X, Z) of transverse disturbance    (=0)
        !   q2m(j): mean along (X, Z) of normal disturbance        (=0)
        do j=1,n2m
            q1m(j)=0.d0
            q2m(j)=0.d0

            do k=ystart(3), min(n3m, yend(3))
                do i=ystart(1), min(n1m, yend(1))
                    q1m(j)=q1m(j)+q1_y(i,j,k)
                    q2m(j)=q2m(j)+q2_y(i,j,k)
                enddo
            enddo

            call MPI_ALLREDUCE (q1m(j), glob_value, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
            q1m(j)=glob_value/dfloat(n1m*n3m)

            call MPI_ALLREDUCE (q2m(j), glob_value, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
            q2m(j)=glob_value/dfloat(n1m*n3m)
        enddo


        do j=2,n2m

            do k=ystart(3), min(n3m, yend(3))
                do i=ystart(1), min(n1m, yend(1))
                    q2_y(i,j,k)=q2_y(i,j,k)-q2m(j)
                enddo
            enddo

        enddo


        do j=1,n2m
            do k=ystart(3), min(n3m, yend(3))
                do i=ystart(1), min(n1m, yend(1))

                    q1_y(i,j,k)=q1_y(i,j,k)-q1m(j)

                enddo
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



        call transpose_y_to_z(q3_y, q3_z)
        call transpose_y_to_x(q1_y, q1_x)


        ! **************************************************************************************
        ! Pressure field initialisation --------------------------------------------------------
        ! **************************************************************************************
115     format(' Bulk velocity of the Poiseuille profile velocity_poisseuille=',e10.3)
116     format(' velocity_X_disturbance=',e10.3,2x,'velocity_Y_disturbance=',e10.3,2x,'velocity_Z_disturbance=',e10.3)
117     format(' Maximum Velocity w max=',e11.4,2x, 'v max=',e11.4,2x,'u max=',e11.4)

        if (nrank==writer_proc) then

            field_generator_path=trim(log_path)//"Fields_generator/"

            open(log_file_id,file=trim(field_generator_path)//'overview.out')

            write(log_file_id,115) velocity_poisseuille
            write(log_file_id,116) velocity_X_disturbance,velocity_Y_disturbance,velocity_Z_disturbance
            write(log_file_id,117) (vmax(l),l=1,3)
            write(log_file_id,*)'Mean along x1 & x3 of the initial perturbation'

            write(log_file_id,*)'j, q1m(j), q2m(j), v3m(j)'
            do j=1,n2m
                write(log_file_id,615)j,q1m(j),q2m(j),v3m(j)
            enddo

615         format(2x,i3,2x,3e12.4)
            close(log_file_id)

        end if

        call MPI_BARRIER(MPI_COMM_WORLD , mpi_err)

        return

!
    end subroutine add_noise_rand_sin

    subroutine add_noise_rand_sin_new(disturbance_intensity)
        use mesh
        use irregular_derivative_coefficients
        use mathematical_constants
        use DNS_settings
        use workspace_view


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

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: rand_u, rand_v, rand_w
        real*8  :: glob_value, ph1, ph2, ph3

        integer :: writer_proc=0
        integer :: log_file_id=20
        character(200)  :: field_generator_path
        integer :: meanXZ_perturbation_file_id=21, zone_id, var_id
        integer, dimension(n2m) :: j_array

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


        do j=1,n2m
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

            do k=ystart(3), min(n3m, yend(3))
                do i=ystart(1), min(n1m, yend(1))

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


     do j=1,n2m
          do k=ystart(3), min(n3m, yend(3))
                do i=ystart(1), min(n1m, yend(1))
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
        do k=ystart(3), min(n3m, yend(3))
            do i=ystart(1), min(n1m, yend(1))
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

        call transpose_y_to_z(q3_y, q3_z)
        call transpose_y_to_x(q1_y, q1_x)


        ! **************************************************************************************
        ! Pressure field initialisation --------------------------------------------------------
        ! **************************************************************************************
115     format(' Bulk velocity of the Poiseuille profile velocity_poisseuille=',e10.3)
116     format(' velocity_X_disturbance=',e10.3,2x,'velocity_Y_disturbance=',e10.3,2x,'velocity_Z_disturbance=',e10.3)
117     format(' Maximum Velocity w max=',e11.4,2x, 'v max=',e11.4,2x,'u max=',e11.4)

        if (nrank==writer_proc) then

            field_generator_path=trim(log_path)//"Fields_generator/"

            open(log_file_id,file=trim(field_generator_path)//'overview.out')

            write(log_file_id,115) velocity_poisseuille
            write(log_file_id,116) velocity_X_disturbance,velocity_Y_disturbance,velocity_Z_disturbance
            write(log_file_id,117) (vmax(l),l=1,3)
            write(log_file_id,*)'Mean along x1 & x3 of the initial perturbation'

            write(log_file_id,*)'j, q1m(j), q2m(j), v3m(j)'
            do j=1,n2m
                write(log_file_id,615)j,q1m(j),q2m(j),v3m(j)
            enddo

615         format(2x,i3,2x,3e12.4)
            close(log_file_id)

        end if

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

        !call add_noise(intensity)
        if(intensity/=0.d0) call add_noise_rand_sin_new(intensity)   ! to do/check


        ! Transpose data in 2D-decomposition stencil
!        call transpose_y_to_x(pr_y, pr_x)
!        call transpose_y_to_z(pr_y, pr_z)

        call transpose_y_to_x(q3_y, q3_x)
        call transpose_y_to_z(q3_y, q3_z)

        call transpose_y_to_x(q1_y, q1_x)
        call transpose_y_to_z(q1_y, q1_z)

        call transpose_y_to_x(q2_y, q2_x)
        call transpose_y_to_z(q2_y, q2_z)

        ! The final 3D field are exported for checking purposes
        call create_snapshot(COMMON_snapshot_path, "INITIALIZATION", q1_y, "W", 2)
        call create_snapshot(COMMON_snapshot_path, "INITIALIZATION", q2_y, "V", 2)
        call create_snapshot(COMMON_snapshot_path, "INITIALIZATION", q3_y, "U", 2)


    end subroutine init_turbulent_field


end module Turbulence_generator
