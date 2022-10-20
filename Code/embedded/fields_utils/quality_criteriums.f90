
module embedded_velocity_quality_criteriums

    use mpi
    use decomp_2d

    use embedded_mesh
    use boundaries
    use embedded_irregular_derivative_coefficients
    use DNS_settings
    use embedded_data

    implicit none
contains


    subroutine perform_stability(q3_z, q2_y, q1_x, cflmax, div_max, div_diff, div_mean)

        use embedded_physical_fields, only:q3_y
        use embedded_velocity_operations
        implicit none

        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3)), intent(in)    :: q3_z
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)), intent(in)    :: q2_y
        real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)), intent(in)    :: q1_x
        real*8, intent(out)                     :: cflmax, div_max, div_diff


        real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3))      :: wc_x, vc_x, uc_x
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3))      :: uc_y, vc_y, wc_y
        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3))      :: uc_z, vc_z, wc_z

        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3)) :: div_z
        real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)) :: div_x


        real*8                  :: cfl, div_min, div_sum, div_mean

        integer k, j, i, mpi_err
        real*8  :: div_max_glob=0.d0, div_min_glob, cflmax_glob=0.d0, kin_energy_glob(3)=0.d0, enstrophy_glob=0.d0

        call perform_velocity_at_center(q3_z, q2_y, q1_x, uc_z, vc_y, wc_x)

        call transpose_x_to_y(wc_x, wc_y, decomp_embedded)
        call transpose_z_to_y(uc_z, uc_y, decomp_embedded)

        cflmax=0.d0

        do k=ystart_e(3), min(n3m, yend_e(3))
            !do i=1,n1m
            do j=1, n2m
                do i=ystart_e(1), min(n1m, yend_e(1))
                    cfl=    abs(uc_y(i,j,k))/ dx3   +   &
                    abs(vc_y(i,j,k))/cell_size_Y(j)    +   &
                    abs(wc_y(i,j,k))/ dx1
                    cflmax=max(cfl, cflmax)
                enddo
            enddo
        enddo

        call MPI_ALLREDUCE (cflmax, cflmax_glob, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)

        cflmax=cflmax_glob*dt

        div_max=-100000.d0
        div_min=100000.d0
        div_sum=0.d0

        div_x=0.d0

        call perform_divergence(div_z, q3_z, q2_y, q1_x)

        do k=1, n3m
            do j=zstart_e(2), min(n2m, zend_e(2))
                do i=zstart_e(1), min(n1m, zend_e(1))
                    div_max=max(div_max, div_z(i,j,k))
                    div_min=min(div_min, div_z(i,j,k))
                    div_sum=div_sum+div_z(i,j,k)
                end do
            end do
        end do

        call MPI_ALLREDUCE (div_max, div_max_glob, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (div_min, div_min_glob, 1, MPI_DOUBLE_PRECISION , MPI_MIN , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (div_sum, div_mean, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        if (abs(div_max_glob)>abs(div_min_glob)) then
            div_max=div_max_glob
        else
            div_max=div_min_glob
        endif


        div_mean=div_mean/(n1m*n2m*n3m)
        div_diff=div_max_glob-div_min_glob

        if (nrank==0) write(*,*)'### 2nd channel div_min_glob', div_min_glob
        if (nrank==0) write(*,*)'### 2nd channel div_max_glob', div_max_glob

        return

    end subroutine

end module embedded_velocity_quality_criteriums
