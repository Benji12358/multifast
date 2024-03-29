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

        integer k,j,i, mpi_err, n1e,n2e,n3e
        integer,save    :: nb=1
        !!!$    integer :: nb_taches
        !!!$    common/nb_total_threads/nb_taches


        ! X orientation ---------------------------------------------------------
        div1_x=0.d0
        !do k=1,n3m
        n2e=(min(n2m, xend(2))-xstart(2))+1
        n3e=(min(n3m, xend(3))-xstart(3))+1

        call D1s_3Dx(q1_x, div1_x, n1,xsize(2),n2e,xsize(3),n3e, dx1, .false., POISSON_VEL_BC1)

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

        call transpose_x_to_y(div1_x, div1_y)
        call transpose_y_to_z(div1_y, div1_z)

        call transpose_y_to_z(div2_y, div2_z)

        div_z=div1_z+div2_z+div3_z


        if (present(divy_mean)) then

            divy_mean=0.d0
            divy_sum=sum(div2_z)



            call MPI_ALLREDUCE (divy_sum, divy_mean, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
            divy_mean=divy_mean/((n1-1)*(n2-1)*(n3-1))

        endif

        return

    end subroutine


    subroutine perform_velocity_at_center(q3_z, q2_y, q1_x, q3c_z, q2c_y, q1c_x)

        implicit none

        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in)    :: q3_z
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)    :: q2_y
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: q1_x

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(out)     :: q1c_x
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(out)     :: q2c_y
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(out)     :: q3c_z

        call D0s_3Dx(q1_x, q1c_x, n1,xsize(2),xsize(2),xsize(3),xsize(3), NS_DEF_BC1)
        call D0s_3Dy(q2_y, q2c_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_DEF_BC2)
        call D0s_3Dz(q3_z, q3c_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_DEF_BC3)

    end subroutine perform_velocity_at_center




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
