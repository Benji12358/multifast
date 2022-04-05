module VELOCITY_operations

    use mpi
    use decomp_2d

    use mesh
    use boundaries
    use irregular_derivative_coefficients
    use DNS_settings
    use schemes3D_interface

    implicit none 
contains 

    subroutine velocity_transpose_x_to_z(q1x, q2x, q3x)

        use physical_fields
        implicit none

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: q3x
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: q2x
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: q1x


        call transpose_x_to_y(q1x, q1_y)
        call transpose_y_to_z(q1_y, q1_z)
        ! q2
        call transpose_x_to_y(q2x, q2_y)
        call transpose_y_to_z(q2_y, q2_z)
        ! q3
        call transpose_x_to_y(q3x, q3_y)
        call transpose_y_to_z(q3_y, q3_z)

    end subroutine velocity_transpose_x_to_z

    subroutine velocity_transpose_z_to_x(q1z, q2z, q3z)

        use physical_fields
        implicit none

        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in)    :: q3z
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in)    :: q2z
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in)    :: q1z


        call transpose_z_to_y(q1z, q1_y)
        call transpose_y_to_x(q1_y, q1_x)
        ! q2
        call transpose_z_to_y(q2z, q2_y)
        call transpose_y_to_x(q2_y, q2_x)
        ! q3
        call transpose_z_to_y(q3z, q3_y)
        call transpose_y_to_x(q3_y, q3_x)

    end subroutine velocity_transpose_z_to_x

    subroutine pressure_gradient_transpose_x_to_z(dphidx1x, dphidx2x, dphidx3x)

        use physical_fields
        implicit none

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: dphidx3x
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: dphidx2x
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: dphidx1x


        call transpose_x_to_y(dphidx1x, dphidx1_y)
        call transpose_y_to_z(dphidx1_y, dphidx1_z)
        ! q2
        call transpose_x_to_y(dphidx2x, dphidx2_y)
        call transpose_y_to_z(dphidx2_y, dphidx2_z)
        ! q3
        call transpose_x_to_y(dphidx3x, dphidx3_y)
        call transpose_y_to_z(dphidx3_y, dphidx3_z)

    end subroutine pressure_gradient_transpose_x_to_z

    subroutine spread_to_all_pencil(q3z, q2y, q1x, dpz)

        use physical_fields
        implicit none

        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in)    :: q3z
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)    :: q2y
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: q1x
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in)    :: dpz

        call transpose_z_to_y(q3z, q3_y)
        call transpose_y_to_x(q3_y, q3_x)

        call transpose_y_to_z(q2y, q2_z)
        call transpose_y_to_x(q2_y, q2_x)

        call transpose_x_to_y(q1x, q1_y)
        call transpose_y_to_z(q1_y, q1_z)

        call transpose_z_to_y(dpz, dp_y)
        call transpose_y_to_x(dp_y, dp_x)

    end subroutine spread_to_all_pencil


    subroutine perform_divergence(div_z, q3_z, q2_y, q1_x, divy_mean)
        use numerical_methods_settings
        use snapshot_writer
        use mpi
        use schemes_interface
        use IBM

        use COMMON_workspace_view, only: COMMON_snapshot_path
        use snapshot_writer
        use HDF5_IO

        implicit none
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in)    :: q3_z
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)    :: q2_y
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: q1_x

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))              :: div1_x
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(out) :: div_z

        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))              :: div3_z, div2_z, div1_z, div_z2

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))   :: div2_y, div1_y
        real*8      :: divy_sum
        real*8, optional :: divy_mean

        integer k,j,i,n, mpi_err, n1e,n2e,n3e
        !!!$    integer :: nb_taches
        !!!$    common/nb_total_threads/nb_taches


        ! X orientation ---------------------------------------------------------
        div1_x=0.d0
        !do k=1,n3m
        n2e=(min(n2m, xend(2))-xstart(2))+1
        n3e=(min(n3m, xend(3))-xstart(3))+1

        call D1s_3Dx(q1_x, div1_x, n1,xsize(2),n2e,xsize(3),n3e, dx1, .false., POISSON_VEL_BC1)

        ! if (IBM_activated) then

        !     do n=1,number_of_objects

        !         if (object_in_current_proc_x_q1(n)) then

        !             do i=i_start_obj_q1(n)-5, i_start_obj_q1(n)
        !                 do j=xstart_ibm_q1(n,2), xend_ibm_q1(n,2)
        !                     do k=xstart_ibm_q1(n,3), xend_ibm_q1(n,3)

        !                         ! remove upward term and add backward term
        !                         div1_x(i,j,k) = div1_x(i,j,k) - (1.d0/dx1) * ( q1_x(i+1,j,k) - q1_x(i,j,k) ) + (1.d0/dx1) * ( q1_x(i-1,j,k) - q1_x(i,j,k) ) 

        !                     enddo
        !                 enddo
        !             enddo

        !         endif

        !     enddo

        ! endif

        ! Y orientation ---------------------------------------------------------

        div2_y=0.d0
        n1e=(min(n1m, yend(1))-ystart(1))+1
        n3e=(min(n3m, yend(3))-ystart(3))+1

        if (use_generic_poisson) then

            call D1s_MULT_3Dy(q2_y, div2_y, ysize(1),n1e,n2,ysize(3),n3e, dx2, .false., POISSON_VEL_BC2, Yc_to_YcTr_for_D1)

        else    ! ONLY FOR CHANNEL FLOWS

            do k=ystart(3), min(n3m, yend(3))   !do k=1,n3m
                do i=ystart(1), min(n1m, yend(1))   !do i=1,n1m
                    div2_y(i,1,k)=(q2_y(i,2,k) - q2_y(i,1,k))*Yc_to_YcTr_for_D1(1)/dx2
                    div2_y(i,n2m,k)=(q2_y(i,n2,k) - q2_y(i,n2m,k))*Yc_to_YcTr_for_D1(n2m)/dx2
                    call D1s_Tamm_MULT(q2_y(i,2:n2m,k), div2_y(i,2:n2m-1,k), n2m-1, dx2, .false., Dirichlet, Yc_to_YcTr_for_D1(2:n2m-1))
                enddo
            enddo
        
        end if

        ! Z orientation ---------------------------------------------------------

        div3_z=0.d0
        n1e=(min(n1m, zend(1))-zstart(1))+1
        n2e=(min(n2m, zend(2))-zstart(2))+1

        call D1s_ACC_3Dz(q3_z, div3_z, zsize(1),n1e,zsize(2),n2e,n3, dx3, .false., POISSON_VEL_BC3)

        ! if (IBM_activated) then

        !     do n=1,number_of_objects

        !         if (object_in_current_proc_z_q3(n)) then

        !             if (k_start_obj_q3(n).gt.1) then

        !                 do k=max(1,k_start_obj_q3(n)-5), k_start_obj_q3(n)
        !                     do j=zstart_ibm_q3(n,2), zend_ibm_q3(n,2)
        !                         do i=zstart_ibm_q3(n,1), zend_ibm_q3(n,1)

        !                             ! remove upward term and add backward term
        !                             div3_z(i,j,k) = div3_z(i,j,k) - (1.d0/dx3) * ( q3_z(i,j,k+1) - q3_z(i,j,k) ) + (1.d0/dx3) * ( q3_z(i,j,k-1) - q3_z(i,j,k) ) 

        !                         enddo
        !                     enddo
        !                 enddo

        !             endif

        !         endif

        !     enddo

        ! endif

        call transpose_x_to_y(div1_x, div1_y)
        call transpose_y_to_z(div1_y, div1_z)

        call transpose_y_to_z(div2_y, div2_z)

        div_z=div1_z+div2_z+div3_z

        ! write(*,*) 'div_z', div_z(:,10,10)

        if (present(divy_mean)) then

            divy_mean=0.d0
            divy_sum=sum(div2_z)

            call MPI_ALLREDUCE (divy_sum, divy_mean, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
            divy_mean=divy_mean/((n1-1)*(n2-1)*(n3-1))

        endif

        return

    end subroutine perform_divergence


    ! subroutine perform_divergence_ibm(div_z_ibm, q3_z, q2_y, q1_x)
    !     use numerical_methods_settings
    !     use snapshot_writer
    !     use mpi
    !     use schemes_interface
    !     use IBM
    !     use IBM_data
    !     use DRP_IBM
    !     use HDF5_IO

    !     use COMMON_workspace_view, only: COMMON_snapshot_path
    !     use snapshot_writer

    !     implicit none

    !     character(200)                :: snaps_dir, snap_dir
    !     character(200)    :: file_path, snap_path
    !     real*8, dimension(zstart_fine(1):zend_fine(1), zstart_fine(2):zend_fine(2), zstart_fine(3):zend_fine(3)), intent(out) :: div_z_ibm

    !     real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in)    :: q3_z
    !     real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)    :: q2_y
    !     real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: q1_x

    !     real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))              :: tmp1_x
    !     real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))              :: tmp1_y, tmp2_y
    !     real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))              :: tmp1_z, tmp2_z, tmp3_z

    !     real*8, dimension(:,:,:), allocatable                                                   :: tmp1_x_ibm
    !     real*8, dimension(:,:,:), allocatable                                                   :: tmp1_y_ibm, tmp2_y_ibm
    !     real*8, dimension(:,:,:), allocatable                                                   :: tmp1_z_ibm, tmp2_z_ibm, tmp3_z_ibm

    !     real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: div_z

    !     integer k,j,i, mpi_err, n1e,n2e,n3e

    !     real*8                  :: div_min, div_sum, div_mean, div_max, div_diff
    !     real*8      :: divy_sum
    !     real*8      :: divy_mean

    !     real*8  :: div_max_glob=0.d0, div_min_glob
    !     integer,save    :: nb=1

    !     !!!$    integer :: nb_taches
    !     !!!$    common/nb_total_threads/nb_taches


    !     ! X orientation ---------------------------------------------------------
    !     tmp1_x=0.d0
    !     !do k=1,n3m
    !     n2e=(min(n2-1, xend(2))-xstart(2))+1
    !     n3e=(min(n3-1, xend(3))-xstart(3))+1

    !     ! first compute on coarsest grid
    !     call D1s_3Dx(q1_x, tmp1_x, n1,xsize(2),n2e,xsize(3),n3e, dx1, .false., POISSON_VEL_BC1)

    !     call transpose_x_to_y(tmp1_x, tmp1_y)
    !     call transpose_y_to_z(tmp1_y, tmp1_z)

    !     allocate(tmp1_x_ibm(xstart_ibm_fcm(1):xend_ibm_fcm(1), xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)))
    !     allocate(tmp1_y_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(2):yend_ibm_fcm(2), ystart_ibm_fcm(3):yend_ibm_fcm(3)))
    !     allocate(tmp1_z_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2), zstart_ibm_fcm(3):zend_ibm_fcm(3)))
    !     tmp1_x_ibm = 0.d0
    !     tmp1_y_ibm = 0.d0
    !     tmp1_z_ibm = 0.d0

    !     n2e = xsize_ibm_fcm(2)
    !     n3e = xsize_ibm_fcm(3)

    !     call D1s_IBM_O2_ACC_3Dx(q1_x_ibm,tmp1_x_ibm,xsize_ibm_fcm(1), xsize_ibm_fcm(2),n2e,xsize_ibm_fcm(3),n3e,i_start_ibm_2nd_fc,i_end_ibm_2nd_fc,dx1_ibm)

    !     call transpose_x_to_y(tmp1_x_ibm, tmp1_y_ibm, decomp_fine)
    !     call transpose_y_to_z(tmp1_y_ibm, tmp1_z_ibm, decomp_fine)

    !     ! Y orientation ---------------------------------------------------------

    !     tmp2_y=0.d0
    !     n1e=(min(n1-1, yend(1))-ystart(1))+1
    !     n3e=(min(n3-1, yend(3))-ystart(3))+1

    !     ! first compute on coarsest grid
    !     if (use_generic_poisson) then

    !         call D1s_MULT_3Dy(q2_y, tmp2_y, ysize(1),n1e,n2,ysize(3),n3e, dx2, .false., POISSON_VEL_BC2, Yc_to_YcTr_for_D1)

    !     else    ! ONLY FOR CHANNEL FLOWS

    !         do k=ystart(3), min(n3-1, yend(3))   !do k=1,n3m
    !             do i=ystart(1), min(n1-1, yend(1))   !do i=1,n1m
    !                 tmp2_y(i,1,k)=(q2_y(i,2,k) - q2_y(i,1,k))*Yc_to_YcTr_for_D1(1)/dx2
    !                 tmp2_y(i,n2-1,k)=(q2_y(i,n2,k) - q2_y(i,n2-1,k))*Yc_to_YcTr_for_D1(n2-1)/dx2
    !                 call D1s_Tamm_MULT(q2_y(i,2:n2-1,k), tmp2_y(i,2:n2-2,k), n2-2, dx2, .false., Dirichlet, Yc_to_YcTr_for_D1(2:n2-2))
    !             enddo
    !         enddo

    !     end if

    !     call transpose_y_to_z(tmp2_y, tmp2_z)

    !     allocate(tmp2_y_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(2):yend_ibm_fcm(2), ystart_ibm_fcm(3):yend_ibm_fcm(3)))
    !     allocate(tmp2_z_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2), zstart_ibm_fcm(3):zend_ibm_fcm(3)))
    !     tmp2_y_ibm = 0.d0
    !     tmp2_z_ibm = 0.d0

    !     n1e = ysize_ibm_fcm(1)
    !     n3e = ysize_ibm_fcm(3)

    !     call D1s_IBM_O2_MULTACC_3Dy(q2_y_ibm,tmp2_y_ibm,ysize_ibm_fcm(1), n1e, ysize_ibm_fcm(2),ysize_ibm_fcm(3),n3e,j_start_ibm_2nd_fc,j_end_ibm_2nd_fc,dx2_ibm,Yc_to_YcTr_for_D1_ibm)

    !     call transpose_y_to_z(tmp2_y_ibm, tmp2_z_ibm,subdecomp_fine_computational)

    !     ! Z orientation ---------------------------------------------------------

    !     tmp3_z=0.d0
    !     n1e=(min(n1, zend(1))-zstart(1))+1
    !     n2e=(min(n2, zend(2))-zstart(2))+1

    !     ! first compute on coarsest grid
    !     call D1s_ACC_3Dz(q3_z, tmp3_z, zsize(1),n1e,zsize(2),n2e,n3, dx3, .false., POISSON_VEL_BC3)

    !     ! and on the fine grid
    !     allocate(tmp3_z_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2), zstart_ibm_fcm(3):zend_ibm_fcm(3)))
    !     tmp3_z_ibm = 0.d0

    !     n1e = zsize_ibm_fcm(1)
    !     n2e = zsize_ibm_fcm(2)

    !     call D1s_IBM_O2_ACC_3Dz(q3_z_ibm,tmp3_z_ibm,zsize_ibm_fcm(1), n1e, zsize_ibm_fcm(2),n2e,zsize_ibm_fcm(3),k_start_ibm_2nd_fc,k_end_ibm_2nd_fc,dx3_ibm)

    !     ! do i=zstart(1), zend(1)
    !     !     do j=zstart(2), zend(2)
    !     !         do k=zstart(3), zend(3)

    !     !             if ((i.ge.i_start_ibm_2nd).and.(i.le.i_end_ibm_2nd).and.(j.ge.j_start_ibm_2nd).and.(j.le.j_end_ibm_2nd).and.(k.ge.k_start_ibm_2nd).and.(k.le.k_end_ibm_2nd)) then

    !     !                 ! 8 fine cells in 1 coarse cell
    !     !                 div_z_ibm(2*i-1,2*j-1,2*k-1)                               = tmp1_z_ibm(2*i-1,2*j-1,2*k-1)                                  &
    !     !                                                                            + tmp2_z_ibm(2*i-1,2*j-1,2*k-1)                                  &
    !     !                                                                            + tmp3_z_ibm(2*i-1,2*j-1,2*k-1)                                  
                        
    !     !                 div_z_ibm(min(2*i,n1_ibm),2*j-1,2*k-1)                     = tmp1_z_ibm(min(2*i,n1_ibm),2*j-1,2*k-1)                        & 
    !     !                                                                            + tmp2_z_ibm(min(2*i,n1_ibm),2*j-1,2*k-1)                        &
    !     !                                                                            + tmp3_z_ibm(min(2*i,n1_ibm),2*j-1,2*k-1)                        
                        
    !     !                 div_z_ibm(2*i-1,min(2*j,n2_ibm),2*k-1)                     = tmp1_z_ibm(2*i-1,min(2*j,n2_ibm),2*k-1)                        & 
    !     !                                                                            + tmp2_z_ibm(2*i-1,min(2*j,n2_ibm),2*k-1)                        &
    !     !                                                                            + tmp3_z_ibm(2*i-1,min(2*j,n2_ibm),2*k-1)                        
                        
    !     !                 div_z_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),2*k-1)           = tmp1_z_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),2*k-1)              & 
    !     !                                                                            + tmp2_z_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),2*k-1)              &
    !     !                                                                            + tmp3_z_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),2*k-1)              
                        
    !     !                 div_z_ibm(2*i-1,2*j-1,min(2*k,n3_ibm))                     = tmp1_z_ibm(2*i-1,2*j-1,min(2*k,n3_ibm))                        & 
    !     !                                                                            + tmp2_z_ibm(2*i-1,2*j-1,min(2*k,n3_ibm))                        & 
    !     !                                                                            + tmp3_z_ibm(2*i-1,2*j-1,min(2*k,n3_ibm))                        
                        
    !     !                 div_z_ibm(2*i-1,min(2*j,n2_ibm),min(2*k,n3_ibm))           = tmp1_z_ibm(2*i-1,min(2*j,n2_ibm),min(2*k,n3_ibm))              & 
    !     !                                                                            + tmp2_z_ibm(2*i-1,min(2*j,n2_ibm),min(2*k,n3_ibm))              &
    !     !                                                                            + tmp3_z_ibm(2*i-1,min(2*j,n2_ibm),min(2*k,n3_ibm))              
                        
    !     !                 div_z_ibm(min(2*i,n1_ibm),2*j-1,min(2*k,n3_ibm))           = tmp1_z_ibm(min(2*i,n1_ibm),2*j-1,min(2*k,n3_ibm))              & 
    !     !                                                                            + tmp2_z_ibm(min(2*i,n1_ibm),2*j-1,min(2*k,n3_ibm))              & 
    !     !                                                                            + tmp3_z_ibm(min(2*i,n1_ibm),2*j-1,min(2*k,n3_ibm))              
                        
    !     !                 div_z_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),min(2*k,n3_ibm)) = tmp1_z_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),min(2*k,n3_ibm))    & 
    !     !                                                                            + tmp2_z_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),min(2*k,n3_ibm))    &
    !     !                                                                            + tmp3_z_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),min(2*k,n3_ibm))    

    !     !             else
    !     !                 ! 8 fine cells in 1 coarse cell
    !     !                 div_z_ibm(2*i-1,2*j-1,2*k-1)                               = tmp1_z(i,j,k) + tmp2_z(i,j,k) + tmp3_z(i,j,k) 
    !     !                 div_z_ibm(min(2*i,n1_ibm),2*j-1,2*k-1)                     = tmp1_z(i,j,k) + tmp2_z(i,j,k) + tmp3_z(i,j,k) 
    !     !                 div_z_ibm(2*i-1,min(2*j,n2_ibm),2*k-1)                     = tmp1_z(i,j,k) + tmp2_z(i,j,k) + tmp3_z(i,j,k) 
    !     !                 div_z_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),2*k-1)           = tmp1_z(i,j,k) + tmp2_z(i,j,k) + tmp3_z(i,j,k) 
    !     !                 div_z_ibm(2*i-1,2*j-1,min(2*k,n3_ibm))                     = tmp1_z(i,j,k) + tmp2_z(i,j,k) + tmp3_z(i,j,k) 
    !     !                 div_z_ibm(2*i-1,min(2*j,n2_ibm),min(2*k,n3_ibm))           = tmp1_z(i,j,k) + tmp2_z(i,j,k) + tmp3_z(i,j,k) 
    !     !                 div_z_ibm(min(2*i,n1_ibm),2*j-1,min(2*k,n3_ibm))           = tmp1_z(i,j,k) + tmp2_z(i,j,k) + tmp3_z(i,j,k) 
    !     !                 div_z_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),min(2*k,n3_ibm)) = tmp1_z(i,j,k) + tmp2_z(i,j,k) + tmp3_z(i,j,k) 

    !     !             endif

    !     !         enddo
    !     !     enddo
    !     ! enddo

    !     return

    ! end subroutine perform_divergence_ibm



    subroutine perform_velocity_at_center(q3_z, q2_y, q1_x, q3c_z, q2c_y, q1c_x)

        implicit none

        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in)    :: q3_z
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)    :: q2_y
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: q1_x

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(out)     :: q1c_x
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(out)     :: q2c_y
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(out)     :: q3c_z


        ! VALGRIND
        q1c_x=0.d0
        q2c_y=0.d0
        q3c_z=0.d0

        call D0s_3Dx(q1_x, q1c_x, n1,xsize(2),xsize(2),xsize(3),xsize(3), NS_DEF_BC1)
        call D0s_3Dy(q2_y, q2c_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_DEF_BC2)
        call D0s_3Dz(q3_z, q3c_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_DEF_BC3)

    end subroutine perform_velocity_at_center


    subroutine perform_passive_scalar_at_center(sca_y, scaC_y, delta_T)

        implicit none

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)    :: sca_y

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(out)     :: scaC_y

        real*8, intent(in)  :: delta_T
        integer :: i,k
        real(kind=8) A(5)

        A(1) =  1.234102595369109d0     /2.d0
        A(2) =  -0.3184044667712d0     /2.d0
        A(3) =  0.11029870162898d0      /2.d0
        A(4) =  -0.030619166038291d0    /2.d0
        A(5) =  0.0046223358114003d0    /2.d0

        ! be carefull
        call D0s_3Dy(sca_y, scaC_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_DEF_BC2)
        do i=ystart(1),min(yend(1),n1m)
            do k=ystart(3),min(yend(3),n3m)

                scaC_y(i,1,k) = 0.5d0*(sca_y(i,2,k)-delta_T)

                scaC_y(i,n2-1,k) = 0.5d0*(delta_T+sca_y(i,n2-1,k))

                scaC_y(i,n2-2,k)  =   9.0d0/16.d0 * (sca_y(i,n2-1,k)+sca_y(i,n2-2,k))     &
                    -   0.0625d0    * (delta_T+sca_y(i,n2-3,k))

                scaC_y(i,n2-3,k)  =   0.59395104312381d0  *   (sca_y(i,n2-2,k)+sca_y(i,n2-3,k))       &
                    -   0.10967656468571d0  *   (sca_y(i,n2-1,k)+sca_y(i,n2-4,k))             &
                    +   0.015725521561903d0 *   (delta_T + sca_y(i,n2-5,k))


                scaC_y(i,n2-5,k)  =   A(1)  *   (sca_y(i,n2-4,k)+sca_y(i,n2-5,k))       &
                    +   A(2)  *   (sca_y(i,n2-3,k)+sca_y(i,n2-6,k))             &
                    +   A(3) *   (sca_y(i,n2-2,k) + sca_y(i,n2-7,k))           &
                    +   A(4) *   (sca_y(i,n2-1,k) + sca_y(i,n2-8,k))           &
                    +   A(5)*   (delta_T + sca_y(i,n2-9,k))

            enddo
        enddo

    end subroutine perform_passive_scalar_at_center




    subroutine perform_vorticity(u1_y, u2_y, u3_y, vort3_y, vort2_y, vort1_y)
        use physical_fields

        implicit none

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)  :: u1_y, u2_y, u3_y
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(out) :: vort3_y, vort2_y, vort1_y

        real*8, dimension(n3)   :: don3
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))      :: vort3_x, vort2_x
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))      :: vort2_z, vort1_z, tmpz

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))      :: u1_x, u2_x, u3_x
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))      :: u3_z, u2_z, u1_z


        call transpose_y_to_z(u2_y, u2_z)
        call transpose_y_to_z(u1_y, u1_z)

        call transpose_y_to_x(u2_y, u2_x)
        call transpose_y_to_x(u3_y, u3_x)

        call D1c_3Dx(u2_x, vort3_x, n1,xsize(2),xsize(2),xsize(3),xsize(3), dx1, .true., NS_DEF_BC1)
        vort3_x(1,:,:)=(u2_x(2,:,:)-u2_x(1,:,:))/dx1
        vort3_x(n1m,:,:)=(u2_x(n1m,:,:)-u2_x(n1m-1,:,:))/dx1

        call transpose_x_to_y(vort3_x, vort3_y)

        call D1c_MULTACC_3Dy(u1_y, vort3_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), -dx2, .true., NS_DEF_BC2, Yc_to_YcTr_for_D1(:))
        vort3_y(:,1,:)=vort3_y(:,1,:) - (u1_y(:,2,:)-u1_y(:,1,:))  *Yc_to_YcTr_for_D1(1) /dx2
        vort3_y(:,n2m,:)=vort3_y(:,n2m,:) - (u1_y(:,n2m,:)-u1_y(:,n2m-1,:))*Yc_to_YcTr_for_D1(n2m) /dx2

        call D1c_3Dx(u3_x, vort2_x, n1,xsize(2),xsize(2),xsize(3),xsize(3), -dx1, .true., NS_DEF_BC1)
        vort2_x(1,:,:)=-(u3_x(2,:,:)-u3_x(1,:,:))/dx1  ! CL
        vort2_x(n1m,:,:)=- (u3_x(n1m,:,:)-u3_x(n1m-1,:,:))/dx1  ! CL

        call transpose_x_to_y(vort2_x, vort2_y)
        call transpose_y_to_z(vort2_y, vort2_z)

        call D1c_3Dz(u1_z, tmpz, zsize(1),zsize(1),zsize(2),zsize(2),n3, dx3, .true., NS_DEF_BC3)
        tmpz(:,:,1)=(u1_z(:,:,2)-u1_z(:,:,1))/dx3             ! CL
        tmpz(:,:,n3m)=(u1_z(:,:,n3m)-u1_z(:,:,n3m-1))/dx3     ! CL
        vort2_z=vort2_z+tmpz

        call transpose_z_to_y(vort2_z, vort2_y)

        call D1c_3Dz(u2_z, vort1_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, -dx3, .true., NS_DEF_BC3)
        vort1_z(:,:,1)=-(u2_z(:,:,2)-u2_z(:,:,1))/dx3             ! CL
        vort1_z(:,:,n3m)=-(u2_z(:,:,n3m)-u2_z(:,:,n3m-1))/dx3     ! CL

        call transpose_z_to_y(vort1_z, vort1_y)

        call D1c_MULTACC_3Dy(u3_y, vort1_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), dx2, .true., NS_DEF_BC2, Yc_to_YcTr_for_D1(:))
        vort1_y(:,1,:)=vort1_y(:,1,:) + (u3_y(:,2,:)-u3_y(:,1,:))*Yc_to_YcTr_for_D1(1) /dx2
        vort1_y(:,n2m,:)=vort1_y(:,n2m,:) + (u3_y(:,n2m,:)-u3_y(:,n2m-1,:))*Yc_to_YcTr_for_D1(n2m) /dx2

    end subroutine perform_vorticity


end module VELOCITY_operations
