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
        
        !-----------------Higher-Order-Derivative-Coeffs
        double precision :: f1d2c0,f1d2c1,f1d2c2  
        double precision :: b1d2c0,b1d2c1,b1d2c2
        
        double precision :: tmpR1,alpha,dz1,dzN 
        double precision :: df21,df31,df41,df32,df42,df43
        double precision :: db12,db13,db14,db23,db24,db34
        
        double precision :: one, two
        double precision :: c1,c2,c3,c4,c5,c6,c7
        
        one=1.0d0; two=2.0d0  
        
        c1=-49.0d0/20.0d0
        c2= 6.0d0
        c3=-15.0d0/2.0d0
        c4= 20.0d0/3.0d0
        c5= 15.0d0/4.0d0
        c6= 6.0d0/5.0d0
        c7= 1.0d0/6.0d0


        vort3_x=0.0d0; vort2_x=0.0d0
        vort1_z=0.0d0; vort2_z=0.0d0
        vort3_y=0.0d0; vort2_y=0.0d0; vort1_y=0.0d0

        call transpose_y_to_z(u2_y, u2_z)
        call transpose_y_to_z(u1_y, u1_z)

        call transpose_y_to_x(u2_y, u2_x)
        call transpose_y_to_x(u3_y, u3_x)

        call D1c_3Dx(u2_x, vort3_x, n1,xsize(2),xsize(2),xsize(3),xsize(3), dx1, .true., NS_DEF_BC1)
        vort3_x(1,:,:)=(u2_x(2,:,:)-u2_x(1,:,:))/dx1
        vort3_x(n1m,:,:)=(u2_x(n1m,:,:)-u2_x(n1m-1,:,:))/dx1

        call transpose_x_to_y(vort3_x, vort3_y)

        call D1c_MULTACC_3Dy(u1_y, vort3_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), -dx2, .true., NS_DEF_BC2, Yc_to_YcTr_for_D1(:))
        !vort3_y(:,1,:)   = vort3_y(:,1,:)   - (u1_y(:,2,:)-u1_y(:,1,:))  *Yc_to_YcTr_for_D1(1) /dx2
        !vort3_y(:,n2m,:) = vort3_y(:,n2m,:) - (u1_y(:,n2m,:)-u1_y(:,n2m-1,:))*Yc_to_YcTr_for_D1(n2m) /dx2

        vort3_y(:,1,:)   = vort3_y(:,1,  :) - ((c1*u1_y(:,1,:)   + c2*u1_y(:,2,:)     + c3*u1_y(:,3,:)     + c4*u1_y(:,4,:)     - c5*u1_y(:,5,:)     + c6*u1_y(:,6,:)     - c7*u1_y(:,7,:))    /dx2)*Yc_to_YcTr_for_D1(1)
        vort3_y(:,n2m,:) = vort3_y(:,n2m,:) - ((c1*u1_y(:,n2m,:) + c2*u1_y(:,n2m-1,:) + c3*u1_y(:,n2m-2,:) + c4*u1_y(:,n2m-3,:) - c5*u1_y(:,n2m-4,:) + c6*u1_y(:,n2m-5,:) - c7*u1_y(:,n2m-6,:))/dx2)*Yc_to_YcTr_for_D1(n2m)


        ! Forward Derivatives    
        !dz1 = yc(2)-yc(1)
        !alpha = (yc(3)-yc(2))/dz1

        ! Forward first derivative second-order accurate
        !tmpR1 = dz1*alpha*(alpha + one)
        !f1d2c0 = (one - (alpha + one)**two)/tmpR1
        !f1d2c1 = ((alpha + one)**two)/tmpR1
        !f1d2c2 = -one / tmpR1

        !vort3_y(:,1,:)   =  vort3_y(:,1,:)    -  (f1d2c0*u1_y(:,1,:)   + f1d2c1*u1_y(:,2,:)     + f1d2c2*u1_y(:,3,:))

        ! Backward Derivative
        !dzN = yc(n2m)-yc(n2m-1)
        !alpha = (yc(n2m-1)-yc(n2m-2))/dzN

        ! Backward first derivative second-order accurate
        !tmpR1 = dzN*alpha*(alpha + one)
        !b1d2c0 = ((alpha + one)**two - one)/tmpR1
        !b1d2c1 = -((alpha + one)**two)/tmpR1  
        !b1d2c2 = one / tmpR1

        !vort3_y(:,n2m,:) =  vort3_y(:,n2m,:)  -  (b1d2c0*u1_y(:,n2m,:) + b1d2c1*u1_y(:,n2m-1,:) + b1d2c2*u1_y(:,n2m-2,:))

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
        !vort1_y(:,1,:)=vort1_y(:,1,:) + (u3_y(:,2,:)-u3_y(:,1,:))*Yc_to_YcTr_for_D1(1) /dx2
        !vort1_y(:,n2m,:)=vort1_y(:,n2m,:) + (u3_y(:,n2m,:)-u3_y(:,n2m-1,:))*Yc_to_YcTr_for_D1(n2m) /dx2
        
        vort1_y(:,1,:)   = vort1_y(:,1,  :) + ((c1*u3_y(:,1,:)   + c2*u3_y(:,2,:)     + c3*u3_y(:,3,:)     + c4*u3_y(:,4,:)     - c5*u3_y(:,5,:)     + c6*u3_y(:,6,:)     - c7*u3_y(:,7,:))    /dx2)*Yc_to_YcTr_for_D1(1)
        vort1_y(:,n2m,:) = vort1_y(:,n2m,:) + ((c1*u3_y(:,n2m,:) + c2*u3_y(:,n2m-1,:) + c3*u3_y(:,n2m-2,:) + c4*u3_y(:,n2m-3,:) - c5*u3_y(:,n2m-4,:) + c6*u3_y(:,n2m-5,:) - c7*u3_y(:,n2m-6,:))/dx2)*Yc_to_YcTr_for_D1(n2m)

        

        ! Forward Derivatives    
        !dz1 = yc(2)-yc(1)
        !alpha = (yc(3)-yc(2))/dz1

        ! Forward first derivative second-order accurate
        !tmpR1 = dz1*alpha*(alpha + one)
        !f1d2c0 = (one - (alpha + one)**two)/tmpR1
        !f1d2c1 = ((alpha + one)**two)/tmpR1
        !f1d2c2 = -one / tmpR1

        !vort1_y(:,1,:) = vort1_y(:,1,:) + (f1d2c0*u3_y(:,1,:) + f1d2c1*u3_y(:,2,:) + f1d2c2*u3_y(:,3,:))

        ! Backward Derivative
        !dzN = yc(n2m)-yc(n2m-1)
        !alpha = (yc(n2m-1)-yc(n2m-2))/dzN

        ! Backward first derivative second-order accurate
        !tmpR1 = dzN*alpha*(alpha + one)
        !b1d2c0 = ((alpha + one)**two - one)/tmpR1
        !b1d2c1 = -((alpha + one)**two)/tmpR1  
        !b1d2c2 = one / tmpR1

        !vort1_y(:,n2m,:) = vort1_y(:,n2m,:) + (b1d2c0*u3_y(:,n2m,:) + b1d2c1*u3_y(:,n2m-1,:) + b1d2c2*u3_y(:,n2m-2,:))




    end subroutine perform_vorticity


    subroutine y_derivative(array_in,array_out,value_at_wall)
        implicit none

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)  :: array_in
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(out) :: array_out
        real*8, dimension(ystart(1):yend(1), ystart(3):yend(3)), intent(in) :: value_at_wall

        !-----------------Higher-Order-Derivative-Coeffs
        double precision :: f1d2c0,f1d2c1,f1d2c2  
        double precision :: b1d2c0,b1d2c1,b1d2c2
        
        double precision :: tmpR1,alpha,dz1,dzN 
        double precision :: df21,df31,df41,df32,df42,df43
        double precision :: db12,db13,db14,db23,db24,db34
        
        double precision :: one, two
        double precision :: c1,c2,c3,c4,c5,c6,c7
        double precision :: al_up,al_down,a1,a2,a3,y10,y20,y21,b1,b2,b3        
        
        integer :: i
        
        one=1.0d0; two=2.0d0
        
        c1=-49.0d0/20.0d0
        c2= 6.0d0
        c3=-15.0d0/2.0d0
        c4= 20.0d0/3.0d0
        c5= 15.0d0/4.0d0
        c6= 6.0d0/5.0d0
        c7= 1.0d0/6.0d0

               
        array_out=0.0d0
        call D1c_MULTACC_3Dy(array_in, array_out, ysize(1),ysize(1),n2,ysize(3),ysize(3), dx2, .true., NS_DEF_BC2, Yc_to_YcTr_for_D1(:))
        
        y10=yc(1)
        y20=yc(2)
        y21=yc(2)-yc(1)

        al_down=y21/y10
        a1=-al_down/(al_down+1)
        a2=(al_down-1)/al_down
        a3=1/(al_down**2+al_down)
        
        a1=a1/y10
        a2=a2/y10
        a3=a3/y10
        
        al_up=(2.0d0-yc(n2m))/(yc(n2m)-yc(n2m-1))
        b1=-al_up/(al_up+1)
        b2=(al_up-1)/al_up
        b3=1/(al_up**2+al_up)

        b1=b1/(y(n2m)-y(n2m-1))
        b2=b2/(y(n2m)-y(n2m-1))
        b3=b3/(y(n2m)-y(n2m-1))

        !array_out(:,1  ,:) = (array_in(:,2,:)   - array_in(:,1,:))     * Yc_to_YcTr_for_D1(1)   /dx2
        !array_out(:,n2m,:) = (array_in(:,n2m,:) - array_in(:,n2m-1,:)) * Yc_to_YcTr_for_D1(n2m) /dx2
        
        !array_out(:,1,:)   = ((c1*array_in(:,1,:)   + c2*array_in(:,2,:)     + c3*array_in(:,3,:)     + c4*array_in(:,4,:)     - c5*array_in(:,5,:)     + c6*array_in(:,6,:)     - c7*array_in(:,7,:))/dx2)    *Yc_to_YcTr_for_D1(1)
        !array_out(:,n2m,:) = ((c1*array_in(:,n2m,:) + c2*array_in(:,n2m-1,:) + c3*array_in(:,n2m-2,:) + c4*array_in(:,n2m-3,:) - c5*array_in(:,n2m-4,:) + c6*array_in(:,n2m-5,:) - c7*array_in(:,n2m-6,:))/dx2)*Yc_to_YcTr_for_D1(n2m)

        do i=ystart(3),yend(3)
          array_out(:,1,i)   =  a1*value_at_wall(:,i) + a2*array_in(:,1,i)   + a3*array_in(:,2,i)
          array_out(:,n2m,i) =  b3*value_at_wall(:,i) + b2*array_in(:,n2m,i) + b1*array_in(:,n2m-1,i)
        end do
                
        ! Forward Derivatives    
        !dz1 = yc(2)-yc(1)
        !alpha = (yc(3)-yc(2))/dz1

        ! Forward first derivative second-order accurate
        !tmpR1 = dz1*alpha*(alpha + one)
        !f1d2c0 = (one - (alpha + one)**two)/tmpR1
        !f1d2c1 = ((alpha + one)**two)/tmpR1
        !f1d2c2 = -one / tmpR1

        !array_out(:,1,:) = f1d2c0*array_in(:,1,:) + f1d2c1*array_in(:,2,:) + f1d2c2*array_in(:,3,:)

        ! Forward Derivatives    
        !dz1 = yc(3)-yc(2)
        !alpha = (yc(4)-yc(3))/dz1

        ! Forward first derivative second-order accurate
        !tmpR1 = dz1*alpha*(alpha + one)
        !f1d2c0 = (one - (alpha + one)**two)/tmpR1
        !f1d2c1 = ((alpha + one)**two)/tmpR1
        !f1d2c2 = -one / tmpR1

        !array_out(:,2,:) = f1d2c0*array_in(:,2,:) + f1d2c1*array_in(:,3,:) + f1d2c2*array_in(:,4,:)


        ! Forward Derivatives    
        !dz1 = yc(4)-yc(3)
        !alpha = (yc(5)-yc(4))/dz1

        ! Forward first derivative second-order accurate
        !tmpR1 = dz1*alpha*(alpha + one)
        !f1d2c0 = (one - (alpha + one)**two)/tmpR1
        !f1d2c1 = ((alpha + one)**two)/tmpR1
        !f1d2c2 = -one / tmpR1

        !array_out(:,3,:) = f1d2c0*array_in(:,3,:) + f1d2c1*array_in(:,4,:) + f1d2c2*array_in(:,5,:)


        ! Backward Derivative
        !dzN = yc(n2m)-yc(n2m-1)
        !alpha = (yc(n2m-1)-yc(n2m-2))/dzN

        ! Backward first derivative second-order accurate
        !tmpR1 = dzN*alpha*(alpha + one)
        !b1d2c0 = ((alpha + one)**two - one)/tmpR1
        !b1d2c1 = -((alpha + one)**two)/tmpR1  
        !b1d2c2 = one / tmpR1

        !array_out(:,n2m,:) = b1d2c0*array_in(:,n2m,:) + b1d2c1*array_in(:,n2m-1,:) + b1d2c2*array_in(:,n2m-2,:)

        ! Backward Derivative
        !dzN = yc(n2m-1)-yc(n2m-2)
        !alpha = (yc(n2m-2)-yc(n2m-3))/dzN

        ! Backward first derivative second-order accurate
        !tmpR1 = dzN*alpha*(alpha + one)
        !b1d2c0 = ((alpha + one)**two - one)/tmpR1
        !b1d2c1 = -((alpha + one)**two)/tmpR1  
        !b1d2c2 = one / tmpR1

        !array_out(:,n2m-1,:) = b1d2c0*array_in(:,n2m-1,:) + b1d2c1*array_in(:,n2m-2,:) + b1d2c2*array_in(:,n2m-3,:)


        ! Backward Derivative
        !dzN = yc(n2m-2)-yc(n2m-3)
        !alpha = (yc(n2m-3)-yc(n2m-4))/dzN

        ! Backward first derivative second-order accurate
        !tmpR1 = dzN*alpha*(alpha + one)
        !b1d2c0 = ((alpha + one)**two - one)/tmpR1
        !b1d2c1 = -((alpha + one)**two)/tmpR1  
        !b1d2c2 = one / tmpR1

        !array_out(:,n2m-2,:) = b1d2c0*array_in(:,n2m-2,:) + b1d2c1*array_in(:,n2m-3,:) + b1d2c2*array_in(:,n2m-4,:)


    end subroutine y_derivative


    subroutine x_derivative(array_in,array_out)
        implicit none

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)  :: array_in
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(out) :: array_out

        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))  :: tmp1
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))  :: tmp2

        array_out=0.0d0
        tmp1=0.0d0
        tmp2=0.0d0

        call transpose_y_to_z(array_in,tmp1)
        call D1c_3Dz(tmp1, tmp2, zsize(1),zsize(1),zsize(2),zsize(2),n3, dx3, .true., NS_DEF_BC3)
        call transpose_z_to_y(tmp2,array_out)

     end subroutine x_derivative


     subroutine z_derivative(array_in,array_out)
        implicit none

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)  :: array_in
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(out) :: array_out

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: tmp1
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: tmp2

        array_out=0.0d0
        tmp1=0.0d0
        tmp2=0.0d0

        call transpose_y_to_x(array_in,tmp1)
        call D1c_3Dx(tmp1, tmp2, n1,xsize(2),xsize(2),xsize(3),xsize(3), dx1, .true., NS_DEF_BC1)
        call transpose_x_to_y(tmp2,array_out)

     end subroutine z_derivative    
     
     subroutine compute_P_wall(v_3D_y,P_3D_y,p_wall)
        implicit none
        
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)  :: v_3D_y, P_3D_y
        real*8, dimension(ystart(1):yend(1), ystart(3):yend(3)), intent(out) :: p_wall
        
        double precision :: f1d2c0,f1d2c1,f1d2c2  
        double precision :: b1d2c0,b1d2c1,b1d2c2
        
        double precision :: f2d2c0,f2d2c1,f2d2c2,f2d2c3
        double precision :: b2d2c0,b2d2c1,b2d2c2,b2d2c3
        
        double precision :: tmpR1,alpha,dz1,dzN 
        double precision :: df21,df31,df41,df32,df42,df43
        double precision :: db12,db13,db14,db23,db24,db34        
        double precision :: one, two
        
        real*8, dimension(ystart(1):yend(1), ystart(3):yend(3)) :: tmp1
     
        one=1.0d0
        two=2.0d0
        
        
        ! Forward Derivatives    
        dz1 = Yc(1)
        alpha = (Yc(2)-Yc(1))/dz1

        ! Forward first derivative second-order accurate
        tmpR1 = dz1*alpha*(alpha + one)
        f1d2c0 = (one - (alpha + one)**two)/tmpR1
        f1d2c1 = ((alpha + one)**two)/tmpR1
        f1d2c2 = -one / tmpR1
        
         
        ! Forward second derivative second-order accurate
        df21 = Yc(1); df32 = Yc(2)-Yc(1)   
        df31 = Yc(2); df42 = Yc(3)-Yc(1)
        df41 = Yc(3); df43 = Yc(3)-Yc(2)

        f2d2c0 = (two*(-df21-df31-df41))/((-df21)*(-df31)*(-df41))         
        f2d2c1 = (two*(-df31-df41)     )/(  df21 *(-df32)*(-df42)) 
        f2d2c2 = (two*(-df21-df41)     )/(  df31 *  df32 *(-df43))
        f2d2c3 = (two*(-df21-df31)     )/(  df41 *  df42 *  df43 )
        

        
        tmp1(:,:) = (f2d2c1*v_3D_y(:,1,:) + f2d2c2*v_3D_y(:,2,:) + f2d2c3*v_3D_y(:,3,:))/ren
        p_wall(:,:) = (tmp1(:,:) - f1d2c2*P_3D_y(:,2,:) - f1d2c1*P_3D_y(:,1,:))/f1d2c0
     
     end subroutine compute_P_wall
     
     
     subroutine compute_OMG_wall(u_3D_y,v_3D_y,w_3D_y,w_wall,ampt,kap,omg,time,OMGX_wall,OMGY_wall,OMGZ_wall)
        implicit none
        
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)  :: u_3D_y,v_3D_y,w_3D_y
        real*8, dimension(ystart(1):yend(1), ystart(3):yend(3)), intent(in)  :: w_wall
        real*8, dimension(ystart(1):yend(1), ystart(3):yend(3)), intent(out) :: OMGX_wall,OMGY_wall,OMGZ_wall
        
        
        double precision :: f1d2c0,f1d2c1,f1d2c2  
        double precision :: b1d2c0,b1d2c1,b1d2c2
        
        double precision :: f2d2c0,f2d2c1,f2d2c2,f2d2c3
        double precision :: b2d2c0,b2d2c1,b2d2c2,b2d2c3
        
        double precision :: tmpR1,alpha,dz1,dzN 
        double precision :: df21,df31,df41,df32,df42,df43
        double precision :: db12,db13,db14,db23,db24,db34        
        double precision :: one, two
        
        double precision, intent(in) :: ampt, kap, omg, time
        integer :: i
        
        one=1.0d0
        two=2.0d0
                
        ! Forward Derivatives    
        dz1 = Yc(1)
        alpha = (Yc(2)-Yc(1))/dz1

        ! Forward first derivative second-order accurate
        tmpR1 = dz1*alpha*(alpha + one)
        f1d2c0 = (one - (alpha + one)**two)/tmpR1
        f1d2c1 = ((alpha + one)**two)/tmpR1
        f1d2c2 = -one / tmpR1
        
        OMGX_wall(:,:) = f1d2c0*w_wall(:,:) + f1d2c1*w_3D_y(:,1,:) + f1d2c2*w_3D_y(:,2,:)
        
        do i=ystart(3),min(yend(3),n3m)
          OMGY_wall(:,i) = -ampt*kap*dcos(kap*Xc(i)-omg*time)
        end do
        
        OMGZ_wall(:,:) = (f1d2c1*u_3D_y(:,1,:) + f1d2c2*u_3D_y(:,2,:)) - OMGY_wall(:,:)
        
     end subroutine compute_OMG_wall
     
     
!     subroutine x_derivative_2D_at_wall(f,d1f)

!        implicit none

!        real*8,  dimension(xstart(1):xend(1),xstart(3),xend(3)), intent(in) :: f
!        real*8,  dimension(xstart(1):xend(1),xstart(3),xend(3)), intent(out):: d1f

!        integer :: i,j,k
!        real(kind=8) A(7)

!        A(1) =  0.91567612963343d0   /h
!        A(2) = -0.34876133244745d0   /h
!        A(3) =  0.14348458980167d0   /h
!        A(4) = -0.050850207601385d0  /h
!        A(5) =  0.013051374066223d0  /h
!        A(6) = -0.0017438790115182d0 /h



!        do k=7,n3-7
!            d1f(:,k)  =   A(6)*(f(:,k+6) - f(:,k-6))     &
!                      +   A(5)*(f(:,k+5) - f(:,k-5))     &
!                      +   A(4)*(f(:,k+4) - f(:,k-4))     &
!                      +   A(3)*(f(:,k+3) - f(:,k-3))     &
!                      +   A(2)*(f(:,k+2) - f(:,k-2))     &
!                      +   A(1)*(f(:,k+1) - f(:,k-1))
!        enddo


!        do k=1, 6
!            d1f(:,k)  =   A(6)*(f(:,k+6) - f(:,mod(n3+k-8,n3-1)+1))  &
!                      +   A(5)*(f(:,k+5) - f(:,mod(n3+k-7,n3-1)+1))  &
!                      +   A(4)*(f(:,k+4) - f(:,mod(n3+k-6,n3-1)+1))  &
!                      +   A(3)*(f(:,k+3) - f(:,mod(n3+k-5,n3-1)+1))  &
!                      +   A(2)*(f(:,k+2) - f(:,mod(n3+k-4,n3-1)+1))  &
!                      +   A(1)*(f(:,k+1) - f(:,mod(n3+k-3,n3-1)+1))
!        enddo

!        do k=n3-6, n3-1
!            d1f(:,k)  =   A(6)*(f(:,mod(k+5,n3-1)+1) - f(:,k-6))     &
!                      +   A(5)*(f(:,mod(k+4,n3-1)+1) - f(:,k-5))     &
!                      +   A(4)*(f(:,mod(k+3,n3-1)+1) - f(:,k-4))     &
!                      +   A(3)*(f(:,mod(k+2,n3-1)+1) - f(:,k-3))     &
!                      +   A(2)*(f(:,mod(k+1,n3-1)+1) - f(:,k-2))     &
!                      +   A(1)*(f(:,mod(k  ,n3-1)+1) - f(:,k-1))
!        enddo



!        return

!    end subroutine x_derivative_2D_at_wall


end module VELOCITY_operations
