module embedded_velocity_operations

    use mpi
    use decomp_2d

    use embedded_mesh
    use boundaries
    use embedded_irregular_derivative_coefficients
    use DNS_settings
    use schemes3D_interface

    implicit none 
contains 

    subroutine spread_to_all_pencil(q3z, q2y, q1x, dpz)

        use embedded_physical_fields
        use embedded_data
        implicit none

        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3)), intent(in)    :: q3z
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)), intent(in)    :: q2y
        real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)), intent(in)    :: q1x
        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3)), intent(in)    :: dpz

        call transpose_z_to_y(q3z, q3_y, decomp_embedded)
        call transpose_y_to_x(q3_y, q3_x, decomp_embedded)

        call transpose_y_to_z(q2y, q2_z, decomp_embedded)
        call transpose_y_to_x(q2_y, q2_x, decomp_embedded)

        call transpose_x_to_y(q1x, q1_y, decomp_embedded)
        call transpose_y_to_z(q1_y, q1_z, decomp_embedded)

        call transpose_z_to_y(dpz, dp_y, decomp_embedded)
        call transpose_y_to_x(dp_y, dp_x, decomp_embedded)

    end subroutine spread_to_all_pencil


    subroutine perform_divergence(div_z, q3_z, q2_y, q1_x, divy_mean)
        use numerical_methods_settings
        use mpi
        use schemes_interface
        use embedded_data

        implicit none
        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3)), intent(in)    :: q3_z
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)), intent(in)    :: q2_y
        real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)), intent(in)    :: q1_x

        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3)), intent(out)    :: div_z

        real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3))              :: div1_x

        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3))              :: div3_z, div2_z, div1_z, div_z2

        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3))              :: div2_y, div1_y
        real*8      :: divy_sum
        real*8, optional :: divy_mean

        integer k,j,i,n, mpi_err, n1e,n2e,n3e
        !!!$    integer :: nb_taches
        !!!$    common/nb_total_threads/nb_taches


        ! X orientation ---------------------------------------------------------
        div1_x=0.d0
        !do k=1,n3m
        n2e=(min(n2m, xend_e(2))-xstart_e(2))+1
        n3e=(min(n3m, xend_e(3))-xstart_e(3))+1

        call D1s_3Dx(q1_x, div1_x, n1,xsize_e(2),n2e,xsize_e(3),n3e, dx1, .false., POISSON_VEL_BC1)

        ! Y orientation ---------------------------------------------------------

        div2_y=0.d0
        n1e=(min(n1m, yend_e(1))-ystart_e(1))+1
        n3e=(min(n3m, yend_e(3))-ystart_e(3))+1

        call D1s_MULT_3Dy(q2_y, div2_y, ysize_e(1),n1e,n2,ysize_e(3),n3e, dx2, .false., POISSON_VEL_BC2, Yc_to_YcTr_for_D1)


        ! Z orientation ---------------------------------------------------------

        div3_z=0.d0
        n1e=(min(n1m, zend_e(1))-zstart_e(1))+1
        n2e=(min(n2m, zend_e(2))-zstart_e(2))+1

        call D1s_ACC_3Dz(q3_z, div3_z, zsize_e(1),n1e,zsize_e(2),n2e,n3, dx3, .false., POISSON_VEL_BC3)

        call transpose_x_to_y(div1_x, div1_y, decomp_embedded)
        call transpose_y_to_z(div1_y, div1_z, decomp_embedded)

        call transpose_y_to_z(div2_y, div2_z, decomp_embedded)

        div_z=div1_z+div2_z+div3_z

        if (present(divy_mean)) then

            divy_mean=0.d0
            divy_sum=sum(div2_z)

            call MPI_ALLREDUCE (divy_sum, divy_mean, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
            divy_mean=divy_mean/((n1-1)*(n2-1)*(n3-1))

        endif 

    end subroutine perform_divergence

    subroutine perform_velocity_at_center(q3_z, q2_y, q1_x, q3c_z, q2c_y, q1c_x)
        use embedded_data
        implicit none

        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3)), intent(in)    :: q3_z
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)), intent(in)    :: q2_y
        real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)), intent(in)    :: q1_x

        real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)), intent(out)     :: q1c_x
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)), intent(out)     :: q2c_y
        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3)), intent(out)     :: q3c_z


        ! VALGRIND
        q1c_x=0.d0
        q2c_y=0.d0
        q3c_z=0.d0

        call D0s_3Dx(q1_x, q1c_x, n1,xsize_e(2),xsize_e(2),xsize_e(3),xsize_e(3), NS_DEF_BC1)
        call D0s_3Dy(q2_y, q2c_y, ysize_e(1),ysize_e(1),n2,ysize_e(3),ysize_e(3), NS_DEF_BC2)
        call D0s_3Dz(q3_z, q3c_z, zsize_e(1),zsize_e(1),zsize_e(2),zsize_e(2),n3, NS_DEF_BC3)

    end subroutine perform_velocity_at_center


end module embedded_velocity_operations