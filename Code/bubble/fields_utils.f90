


module BUBBLE_fields_utils

    use mesh
    use boundaries
    use irregular_derivative_coefficients
    use schemes3D_interface
    use DNS_settings

    use decomp_2d
    use mpi

    implicit none

contains


        subroutine perform_vorticity_x(u1_x, u2_x, u3_x, vort1_x, vort2_x, vort3_x)
            use physical_fields

            implicit none

            real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)  :: u1_x, u2_x, u3_x
            real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(out) :: vort3_x, vort2_x, vort1_x

            real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))      :: vort3_y, vort2_y, vort1_y
            real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))      :: vort2_z, vort1_z, tmpz

            real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))      :: u1_y, u2_y, u3_y
            real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))      :: u3_z, u2_z, u1_z

            vort1_x=0.d0
            vort2_x=0.d0
            vort3_x=0.d0

            call transpose_x_to_y(u2_x, u2_y)
            call transpose_y_to_z(u2_y, u2_z)

            call transpose_x_to_y(u1_x, u1_y)
            call transpose_y_to_z(u1_y, u1_z)

            call transpose_x_to_y(u3_x, u3_y)

            call D1c_3Dx(u2_x, vort3_x, n1,xsize(2),xsize(2),xsize(3),xsize(3), dx1, .true., NS_DEF_BC1)
            vort3_x(1,:,:)=(u2_x(2,:,:)-u2_x(1,:,:))/dx1
            vort3_x(n1m,:,:)=(u2_x(n1m,:,:)-u2_x(n1m-1,:,:))/dx1

            call transpose_x_to_y(vort3_x, vort3_y)

            call D1c_MULTACC_3Dy(u1_y, vort3_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), -dx2, .true., NS_DEF_BC2, Yc_to_YcTr_for_D1(:))
            vort3_y(:,1,:)=vort3_y(:,1,:) - (u1_y(:,2,:)-u1_y(:,1,:))  *Yc_to_YcTr_for_D1(1) /dx2
            vort3_y(:,n2m,:)=vort3_y(:,n2m,:) - (u1_y(:,n2m,:)-u1_y(:,n2m-1,:))*Yc_to_YcTr_for_D1(n2m) /dx2

            call transpose_y_to_x(vort3_y, vort3_x)

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
            call transpose_y_to_x(vort2_y, vort2_x)

            call D1c_3Dz(u2_z, vort1_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, -dx3, .true., NS_DEF_BC3)
            vort1_z(:,:,1)=-(u2_z(:,:,2)-u2_z(:,:,1))/dx3             ! CL
            vort1_z(:,:,n3m)=-(u2_z(:,:,n3m)-u2_z(:,:,n3m-1))/dx3     ! CL

            call transpose_z_to_y(vort1_z, vort1_y)

            call D1c_MULTACC_3Dy(u3_y, vort1_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), dx2, .true., NS_DEF_BC2, Yc_to_YcTr_for_D1(:))
            vort1_y(:,1,:)=vort1_y(:,1,:) + (u3_y(:,2,:)-u3_y(:,1,:))*Yc_to_YcTr_for_D1(1) /dx2
            vort1_y(:,n2m,:)=vort1_y(:,n2m,:) + (u3_y(:,n2m,:)-u3_y(:,n2m-1,:))*Yc_to_YcTr_for_D1(n2m) /dx2

            call transpose_y_to_x(vort1_y, vort1_x)

        end subroutine perform_vorticity_x

    subroutine perform_velocity_at_center_x(q3_z, q2_y, q1_x, q3c_x, q2c_x, q1c_x)

        implicit none

        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in)    :: q3_z
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)    :: q2_y
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: q1_x

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(out)     :: q1c_x
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(out)     :: q2c_x
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(out)     :: q3c_x

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))                  :: q2c_y, q3c_y
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))                  :: q3c_z

        q1c_x=0.d0
        q2c_x=0.d0
        q3c_x=0.d0

        call D0s_3Dx(q1_x, q1c_x, n1,xsize(2),xsize(2),xsize(3),xsize(3), NS_DEF_BC1)

        call D0s_3Dy(q2_y, q2c_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_DEF_BC2)
        call transpose_y_to_x(q2c_y, q2c_x)

        call D0s_3Dz(q3_z, q3c_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_DEF_BC3)
        call transpose_z_to_y(q3c_z, q3c_y)
        call transpose_y_to_x(q3c_y, q3c_x)

    end subroutine perform_velocity_at_center_x

end module BUBBLE_fields_utils

