
module VELOCITY_quality_criteriums

    use mpi
    use decomp_2d

    use mesh
    use boundaries
    use irregular_derivative_coefficients
    use DNS_settings

    implicit none
contains


    subroutine perform_global_divergence(div1, div2, q3_z, q2_y, q1_x)
        use numerical_methods_settings
        use mpi
        implicit none
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in)    :: q3_z
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)    :: q2_y
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: q1_x
        integer :: mpi_err

        real*8 :: div1, div1_loc
        real*8 :: div2, div2_loc

        integer k,j,i

        ! X orientation ---------------------------------------------------------
        div1_loc=0.d0
        !do k=1,n3m
        do k=xstart(3), min(n3m, xend(3))
            do j=xstart(2), min(n2m, xend(2))
                div1_loc=div1_loc+q1_x(n1,j,k)-q1_x(1,j,k)
            enddo
        enddo

        call MPI_ALLREDUCE (div1_loc, div1, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        ! Y orientation ---------------------------------------------------------

        !do k=1,n3m
        div2_loc=0.d0
        do k=ystart(3), min(n3m, yend(3))
            do i=ystart(1), min(n1m, yend(1))
                div2_loc=div2_loc+q2_y(i,n2,k)-q2_y(i,1,k)
            enddo
        enddo
        call MPI_ALLREDUCE (div2_loc, div2, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        return

    end subroutine

    subroutine perform_stability(q3_z, q2_y, q1_x, cflmax, div_max, div_diff, div_mean)

        use physical_fields, only:q3_y
        use VELOCITY_operations

        implicit none

        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in)    :: q3_z
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)    :: q2_y
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: q1_x
        real*8, intent(out)                     :: cflmax, div_max, div_diff


        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))      :: wc_x, vc_x, uc_x
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))      :: uc_y, vc_y, wc_y
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))      :: uc_z, vc_z, wc_z

        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: div_z
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: div_x


        real*8                  :: cfl, div_min, div_sum, div_mean

        integer k, j, i, mpi_err
        real*8  :: div_max_glob=0.d0, div_min_glob, cflmax_glob=0.d0, kin_energy_glob(3)=0.d0, enstrophy_glob=0.d0

        call perform_velocity_at_center(q3_z, q2_y, q1_x, uc_z, vc_y, wc_x)

        call transpose_x_to_y(wc_x, wc_y)
        call transpose_z_to_y(uc_z, uc_y)

        cflmax=0.d0

        do k=ystart(3), min(n3m, yend(3))
            !do i=1,n1m
            do j=1, n2m
                do i=ystart(1), min(n1m, yend(1))
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
            do j=zstart(2), min(n2m, zend(2))
                do i=zstart(1), min(n1m, zend(1))
                    div_max=max(div_max, div_z(i,j,k))
                    div_min=min(div_min, div_z(i,j,k))
                    div_sum=div_sum+div_z(i,j,k)
                end do
            end do
        end do

!        write(*,*)'div_min', div_min
 !       write(*,*)'div_max', div_max

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

        if (nrank==0) write(*,*)'div_min_glob', div_min_glob
        if (nrank==0) write(*,*)'div_max_glob', div_max_glob

        return

    end subroutine

    subroutine perform_maxvel(vmax)

        use physical_fields

        implicit none

        real*8,intent(out)      :: vmax(3)
        real*8                  :: glob_value
        integer                 :: mpi_err


        vmax(1)=MAXVAL(abs(q1_x))
        vmax(2)=MAXVAL(abs(q2_y))
        vmax(3)=MAXVAL(abs(q3_x))
        !if(nrank==0) write(*,*)'  vmax(1,2,3 ',(vmax(l),l=1,3)

        call MPI_ALLREDUCE (vmax(1), glob_value, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
        vmax(1)=glob_value

        call MPI_ALLREDUCE (vmax(2), glob_value, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
        vmax(2)=glob_value

        call MPI_ALLREDUCE (vmax(3), glob_value, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
        vmax(3)=glob_value

    end subroutine perform_maxvel

end module VELOCITY_quality_criteriums
