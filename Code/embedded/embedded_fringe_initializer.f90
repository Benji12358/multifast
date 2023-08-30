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
        integer :: j, k, mpi_err
        logical :: pair_n2
        real*8 :: delta, ut1
        real*8, dimension(n2) :: tmp_profile

        if (inflow_type==3) then
            call inflow_from_prev_run
        else
            call default_inflow
        endif

        ! Assume streamwise direction is 1
        ! so streamwise velocity is q1

        if (SCA_state.eq.1) then

            if (init_type==KAWAMURA_INIT) then
                
                ! The Kawamura problem is solved (1998)
                pair_n2 = (mod(n2, 2)==0)

                tmp_profile = 0.d0
                do j=1,n2/2
                    tmp_profile(j) = prandtl*dsqrt(2*ren) * (1/8.d0 * Yc(j) **4 - 1/2.d0 * Yc(j) **3 + Yc(j))
                    tmp_profile(n2-j) = prandtl*dsqrt(2*ren) * (1/8.d0 * Yc(j) **4 - 1/2.d0 * Yc(j) **3 + Yc(j))
                enddo

                if (.not.pair_n2) then
                    j=n2/2+1
                    tmp_profile(j) = prandtl*dsqrt(2*ren) * (1/8.d0 * Yc(j) **4 - 1/2.d0 * Yc(j) **3 + Yc(j))
                endif

                do j=xstart_e(2),xend_e(2)
                    sca_inflow(j,:) = tmp_profile(j)
                enddo

            else

                ! The LYONS problem is solved (1991)
                do k = xstart_e(3), min(n3-1, xend_e(3))
                    do j = xstart_e(2), min(n2-1, xend_e(2))
                        sca_inflow(j,k)=(-delta_T)+(Yc(j)/L2)*2.d0*delta_T
                    end do
                end do

            endif

        endif

        ! USING BLASIUS BOUNDARY LAYER
        delta = Y(5) ! profile on 4 nodes

        pair_n2 = (mod(n2, 2)==0)

        tmp_profile = 0.d0
        do j=1,n2/2
            if (Yc(j)<delta) then
                tmp_profile(j) = 2.d0 * (Yc(j)/delta) - 2.d0 * (Yc(j)/delta)**3 + 1.d0 * (Yc(j)/delta)**4
                tmp_profile(n2-j) = 2.d0 * (Yc(j)/delta) - 2.d0 * (Yc(j)/delta)**3 + 1.d0 * (Yc(j)/delta)**4
            else
                tmp_profile(j) = 1.d0
                tmp_profile(n2-j) = 1.d0
            endif
        enddo

        if (pair_n2) then
            tmp_profile(n2/2+1) = 1.d0
        endif

        do j=xstart_e(2),xend_e(2)
            inflow_profile(j) = tmp_profile(j)
        enddo

        ! Perform the flowrate at the inlet...
        ut1=0.d0
        do k=xstart_e(3),min(xend_e(3),n3-1)
            do j=xstart_e(2),min(xend_e(2),n2-1)
                ut1=ut1+inflow_profile(j)*cell_size_Y(j)*dx3
            enddo
        enddo

        u_bulk_theo=0.d0
        call MPI_ALLREDUCE(ut1, u_bulk_theo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
        u_bulk_theo=u_bulk_theo/(L3*L2)

        return

        contains

            subroutine default_inflow()
                implicit none

                ! The embedded module is intended to work with a square input

                select case (inflow_type)
                    case (SQUARE_INFLOW)

                      if (nrank==0) write(*,*)'use default square inflow for embedded channel'

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
