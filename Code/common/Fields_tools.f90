module COMMON_fieldstools

    use mpi
    use decomp_2d

    implicit none

contains

    subroutine perform_checksum_x(tab, cs)
        use mesh

        implicit none
        real*8,dimension (xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in)    :: tab
        real*8, intent(out)                                                                     :: cs
        real*8  :: cs_loc
        integer k,j,i, mpi_err

        cs_loc=0.d0
        cs=0.d0

        do k=xstart(3), min(n3m, xend(3))
            do j=xstart(2), min(n2m, xend(2))
                do i=1, n1m
                    cs_loc=cs_loc+abs(tab(i,j,k))
                enddo
            enddo
        enddo

        call MPI_ALLREDUCE (cs_loc, cs, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        cs=cs/(n1m*n2m*n3m)

    end subroutine perform_checksum_x

    subroutine perform_checksum_y(tab, cs)
        use mesh

        implicit none
        real*8,dimension (ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)), intent(in)    :: tab
        real*8, intent(out)                                                                     :: cs
        real*8  :: cs_loc
        integer k,j,i, mpi_err

        cs_loc=0.d0
        cs=0.d0

        do k=ystart(3), min(n3m, yend(3))
            do j=ystart(2), min(n2m, yend(2))
                do i=ystart(1), min(n1m, yend(1))
                    cs_loc=cs_loc+abs(tab(i,j,k))
                enddo
            enddo
        enddo

        call MPI_ALLREDUCE (cs_loc, cs, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        cs=cs/(n1m*n2m*n3m)

    end subroutine perform_checksum_y



    subroutine perform_checksum2D_x(tab, cs)
        use mesh

        implicit none
        real*8,dimension (xstart(2):xend(2),xstart(3):xend(3)), intent(in)    :: tab
        real*8, intent(out)                                                                     :: cs
        real*8  :: cs_loc
        integer k,j, mpi_err

        cs_loc=0.d0
        cs=0.d0

        do k=xstart(3), min(n3m, xend(3))
            do j=xstart(2), min(n2m, xend(2))
                    cs_loc=cs_loc+abs(tab(j,k))
            enddo
        enddo

        call MPI_ALLREDUCE (cs_loc, cs, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        cs=cs/(n2m*n3m)

    end subroutine perform_checksum2D_x

end module COMMON_fieldstools
