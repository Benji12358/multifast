module embedded_fringe_initializer
    use decomp_2d
    use mpi
    use embedded_mesh
    use embedded_fringe_data
    use DNS_settings
    use embedded_physical_fields
    use embedded_data

    use embedded_scalar_data

!
    implicit none
!
contains
!
!
    subroutine get_inflow()

        implicit none
        integer :: j, k

        if (inflow_type==3) then
            call inflow_from_prev_run
        else
            call default_inflow
        endif

        ! Assume streamwise direction is 1
        ! so streamwise velocity is q1

        if (SCA_state.eq.1) then

            do k = xstart_e(3), min(n3-1, xend_e(3))
                do j = xstart_e(2), min(n2-1, xend_e(2))
                    sca_inflow(j,k)=(-delta_T)+(Yc(j)/L2)*2.d0*delta_T
                end do
            end do

        endif

        return

        contains

            subroutine default_inflow()
                implicit none

                ! The embedded module is intended to work with a square input

                select case (inflow_type)
                    case (SQUARE_INFLOW)

                      if (nrank==0) write(*,*)'use default square inflow for 2nd channel'

                      u_bulk = 1.d0

                endselect

            end subroutine default_inflow

            ! TO CHECK
            subroutine inflow_from_prev_run()

                implicit none
                real*8      :: u_bulk_loc
                integer     :: k,j,mpi_err

                u_bulk_loc=0.d0

                do k=xstart_e(3),min(xend_e(3),n3-1)
                    do j=xstart_e(2),min(xend_e(2),n2-1)

                        u_bulk_loc = u_bulk_loc + q1_x(1,j,k)*(Y(j+1)-Y(j))*dx3

                    enddo
                enddo

                call MPI_ALLREDUCE(u_bulk_loc, u_bulk, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
                u_bulk=u_bulk/(L3*L2)

            end subroutine inflow_from_prev_run

    end subroutine get_inflow

end module embedded_fringe_initializer