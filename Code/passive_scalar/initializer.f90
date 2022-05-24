module SCALAR_field_generator
    use decomp_2d
    use mesh

    use SCALAR_data

    implicit none

contains

    subroutine generate_fields(sca_x, sca_y, sca_z, sca_down, sca_up)
        use physical_fields, only: q1_y, q3_y 
        use mesh, only: Yc
        use mpi
        use DNS_settings, only: ren, streamwise
        use IBM_settings, only: IBM_activated
        use COMMON_workspace_view, only: COMMON_snapshot_path
        use snapshot_writer, only: create_snapshot
        implicit none
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: sca_x
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))  :: sca_y
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))  :: sca_z
        real*8                                                                      :: sca_down, sca_up
        real*8                                                                      :: rentau_loc, rentau
        logical                                                 :: pair_n2
        real*8                                                  :: pr

        real*8  :: delta

        integer i, j, k
        integer mpi_err

        delta=sca_up-sca_down

        sca_y=0.d0

        if (init_type==KAWAMURA_INIT) then

            if (streamwise==3) then

                do i = ystart(1), min(n1-1, yend(1))
                    do k = ystart(3), min(n3-1, yend(3))

                        rentau_loc = rentau_loc + dsqrt(ren * q3_y(i,1,k)/Yc(1))
                        rentau_loc = rentau_loc + dsqrt(ren * q3_y(i,n2-1,k)/Yc(1))

                    enddo
                enddo

                call MPI_ALLREDUCE (rentau_loc, rentau, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

                rentau=rentau/dfloat(2*n1m*n3m)

                do k = ystart(3), min(n3-1, yend(3))
                    do i = ystart(1), min(n1-1, yend(1))

                        sca_y(i,:,k)=( renprandtl / rentau ) * q3_y(i,:,k)

                    end do
                end do
                
            elseif (streamwise==1) then

                ! do i = ystart(1), min(n1-1, yend(1))
                !     do k = ystart(3), min(n3-1, yend(3))

                !         rentau_loc = rentau_loc + dsqrt(ren * q1_y(i,1,k)/Yc(1))
                !         rentau_loc = rentau_loc + dsqrt(ren * q1_y(i,n2-1,k)/Yc(1))

                !     enddo
                ! enddo

                ! call MPI_ALLREDUCE (rentau_loc, rentau, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

                ! rentau=rentau/dfloat(2*n1m*n3m)

                ! do k = ystart(3), min(n3-1, yend(3))
                !     do i = ystart(1), min(n1-1, yend(1))

                !         sca_y(i,:,k)=( renprandtl / rentau ) * q1_y(i,:,k)

                !     end do
                ! end do

                pair_n2 = (mod(n2, 2)==0)
                pr = renprandtl/ren

                do j=ystart(2),n2/2
                    sca_y(:,j,:) = pr*dsqrt(2*ren) * (1/8.d0 * Yc(j) **4 - 1/2.d0 * Yc(j) **3 + Yc(j))
                    sca_y(:,n2-j,:) = pr*dsqrt(2*ren) * (1/8.d0 * Yc(j) **4 - 1/2.d0 * Yc(j) **3 + Yc(j))
                enddo

                if (.not.pair_n2) then
                    j=n2/2+1
                    sca_y(:,j,:) = pr*dsqrt(2*ren) * (1/8.d0 * Yc(j) **4 - 1/2.d0 * Yc(j) **3 + Yc(j))
                endif

                ! write(*,*) 'theta^+', sca_y(1,:,1)

            endif

        else if (init_type==CONSTANT_HEAT_FLUX) then

            pair_n2 = (mod(n2, 2)==0)

            do j=ystart(2),n2/2
                sca_y(:,j,:) = sca_down+(Yc(j)/L2)*delta
                sca_y(:,n2-j,:) = sca_down+(Yc(j)/L2)*delta
            enddo

            if (.not.pair_n2) then
                j=n2/2+1
                sca_y(:,j,:) = sca_down+(Yc(j)/L2)*delta
            endif

        else

            do k = ystart(3), min(n3-1, yend(3))
                do i = ystart(1), min(n1-1, yend(1))
                    do j = 1, n2-1
                        sca_y(i,j,k)=sca_down+(Yc(j)/L2)*delta
                        !sca_y(i,j,k)=0.d0
                    end do
                end do
            end do

        endif

		if(nrank==0) write(*,*)'SCALAR_generate_fields: Tinf Tsup', sca_down, sca_up, delta

        sca_wall10=0.d0
        sca_wall11=0.d0

        if (init_type==KAWAMURA_INIT) then
            sca_wall20=0.d0
            sca_wall21=0.d0
        else if (init_type==CONSTANT_HEAT_FLUX) then
            sca_wall20=sca_down
            sca_wall21=sca_down
        else
            sca_wall20=sca_down
            sca_wall21=sca_up
        endif

        sca_wall30=0.d0
        sca_wall31=0.d0

        ! if (IBM_activated) call clear_temperature_IBM(sca_down, sca_up)

        call transpose_y_to_x(sca_y, sca_x)
        call transpose_y_to_z(sca_y, sca_z)

        ! The final 3D field are exported for checking purposes
        call create_snapshot(COMMON_snapshot_path, "INITIALIZATION", sca_y, "sca", 2)

    end subroutine generate_fields

    subroutine generate_from_file(sca_x, sca_y, sca_z)
        use HDF5_IO
        use start_settings, only:start_it, half_length_inflow
        use COMMON_workspace_view, only: COMMON_inflow_path
        use DNS_settings, only: nx_start

        implicit none
        real*8, dimension(:,:,:), allocatable   :: sca_tmp
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: sca_x
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))  :: sca_y
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))  :: sca_z

        character(200)       :: current_inflow_path

        integer i

        current_inflow_path = trim(COMMON_inflow_path)//'sca1'

        if (nrank==0) write(*,*)
        if (nrank==0) write(*,*) "Reading inflow from file:", trim(current_inflow_path)//".h5"
        if (nrank==0) write(*,*) "The inflow is get from nx =", nx_start

        if (half_length_inflow.eq.1) then

            allocate(sca_tmp(xstart(1):xend(1)/2 + 1, xstart(2):xend(2), xstart(3):xend(3)))
            sca_tmp=0.d0

            call hdf_read_3Dfield(current_inflow_path, sca_tmp, "sca1", nx_global/2 + 1, ny_global, nz_global, xstart(1), xend(1)/2 + 1, xstart(2),xend(2), xstart(3),xend(3))

        else 

            allocate(sca_tmp(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            sca_tmp=0.d0

            call hdf_read_3Dfield(current_inflow_path, sca_tmp, "sca1", nx_global, ny_global, nz_global, xstart(1), xend(1), xstart(2),xend(2), xstart(3),xend(3))

        endif

        do i=xstart(1), xend(1)
            sca_x(i,:,:) = sca_tmp(nx_start,:,:)
        enddo

        call transpose_y_to_x(sca_y, sca_x)
        call transpose_y_to_z(sca_y, sca_z)

        if(nrank==0) write(*,*)'Scalar fields generated from file'

        sca_wall10=0.d0
        sca_wall11=0.d0

        sca_wall20=sca_y(ystart(1),ystart(2),ystart(3))
        sca_wall21=sca_y(ystart(1),yend(2),ystart(3))

        sca_wall30=0.d0
        sca_wall31=0.d0

    end subroutine generate_from_file

    subroutine clear_temperature_IBM(sca_down, sca_up)
        use IBM_settings

        use mpi
        use decomp_2d

        implicit none

        integer             :: i,j,k,n
        real*8              :: sca_down, sca_up


        do n=1,number_of_objects

            if (object_in_current_proc_y_q2(n)) then

                do i=ystart_ibm_q2(n,1), yend_ibm_q2(n,1)
                    do j=j_start_obj_q2(n),j_end_obj_q2(n)
                        do k=ystart_ibm_q2(n,3), yend_ibm_q2(n,3)

                        if (j.le.n2/2) then
                            sca_y(i,j,k) = sca_down

                        else
                            sca_y(i,j,k) = sca_up

                        endif

                        enddo
                    enddo
                enddo

            endif

        enddo

        return

!
    end subroutine clear_temperature_IBM

end module SCALAR_field_generator
