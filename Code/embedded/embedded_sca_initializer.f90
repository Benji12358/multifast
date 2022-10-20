module embedded_scalar_field_generator
    use decomp_2d
    use embedded_mesh

    use embedded_scalar_data
    use embedded_data

    implicit none

contains

    subroutine generate_fields(sca_x, sca_y, sca_z, sca_down, sca_up)
        use mpi
        use DNS_settings, only: ren, streamwise
        use COMMON_workspace_view, only: COMMON_snapshot_path
        use embedded_snapshot_writer, only: create_snapshot
        implicit none
        real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3))  :: sca_x
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3))  :: sca_y
        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3))  :: sca_z
        real*8                                                                      :: sca_down, sca_up
        real*8                                                                      :: rentau_loc, rentau
        logical                                                 :: pair_n2
        real*8                                                  :: pr

        real*8  :: delta

        integer i, j, k
        integer mpi_err

        delta=sca_up-sca_down

        sca_y=0.d0

        ! For now, only the LYONS problem is solved (1991)
        do k = ystart_e(3), min(n3-1, yend_e(3))
            do i = ystart_e(1), min(n1-1, yend_e(1))
                do j = 1, n2-1
                    sca_y(i,j,k)=sca_down+(Yc(j)/L2)*delta
                    !sca_y(i,j,k)=0.d0
                end do
            end do
        end do

		if(nrank==0) write(*,*)'### 2nd channel SCALAR_generate_fields: Tinf Tsup', sca_down, sca_up, delta

        sca_wall10=0.d0
        sca_wall11=0.d0

        sca_wall20=sca_down
        sca_wall21=sca_up

        sca_wall30=0.d0
        sca_wall31=0.d0

        call transpose_y_to_x(sca_y, sca_x, decomp_embedded)
        call transpose_y_to_z(sca_y, sca_z, decomp_embedded)

        ! The final 3D field are exported for checking purposes
        call create_snapshot(COMMON_snapshot_path, "EMBEDDED_INIT", sca_y, "sca", 2)

    end subroutine generate_fields

end module embedded_scalar_field_generator
