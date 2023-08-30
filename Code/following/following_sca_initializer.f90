module following_scalar_field_generator
    use decomp_2d
    use following_mesh

    use following_scalar_data
    use following_data

    implicit none

contains

    subroutine generate_fields(sca_x, sca_y, sca_z, sca_down, sca_up)
        use mpi
        use DNS_settings, only: ren, streamwise
        use COMMON_workspace_view, only: COMMON_snapshot_path
        use following_snapshot_writer, only: create_snapshot
        implicit none
        real*8, dimension(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3))  :: sca_x
        real*8, dimension(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3))  :: sca_y
        real*8, dimension(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3))  :: sca_z
        real*8                                                                      :: sca_down, sca_up
        real*8                                                                      :: rentau_loc, rentau
        logical                                                 :: pair_n2
        real*8                                                  :: pr

        real*8  :: delta

        integer i, j, k
        integer mpi_err

        delta=sca_up-sca_down

        sca_y=0.d0

        if (init_type==KAWAMURA_INIT) then
            
            ! The Kawamura problem is solved (1998)
            pair_n2 = (mod(n2, 2)==0)

            do j=ystart_f(2),n2/2
                sca_y(:,j,:) = prandtl*dsqrt(2*ren) * (1/8.d0 * Yc(j) **4 - 1/2.d0 * Yc(j) **3 + Yc(j))
                sca_y(:,n2-j,:) = prandtl*dsqrt(2*ren) * (1/8.d0 * Yc(j) **4 - 1/2.d0 * Yc(j) **3 + Yc(j))
            enddo

            if (.not.pair_n2) then
                j=n2/2+1
                sca_y(:,j,:) = prandtl*dsqrt(2*ren) * (1/8.d0 * Yc(j) **4 - 1/2.d0 * Yc(j) **3 + Yc(j))
            endif

        else

            ! The LYONS problem is solved (1991)
            do k = ystart_f(3), min(n3-1, yend_f(3))
                do i = ystart_f(1), min(n1-1, yend_f(1))
                    do j = 1, n2-1
                        sca_y(i,j,k)=sca_down+(Yc(j)/L2)*delta
                        !sca_y(i,j,k)=0.d0
                    end do
                end do
            end do

        endif

        if(nrank==0) write(*,*)'SCALAR_generate_fields: Tinf Tsup', sca_down, sca_up, delta

        sca_wall10=0.d0
        sca_wall11=0.d0

        if (init_type==KAWAMURA_INIT) then
            sca_wall20=0.d0
            sca_wall21=0.d0
        else
            sca_wall20=sca_down
            sca_wall21=sca_up
        endif

        sca_wall30=0.d0
        sca_wall31=0.d0

        call transpose_y_to_x(sca_y, sca_x, decomp_following)
        call transpose_y_to_z(sca_y, sca_z, decomp_following)

        ! The final 3D field are exported for checking purposes
        call create_snapshot(COMMON_snapshot_path, "FOLLOWING_INIT", sca_y, "sca", 2)

    end subroutine generate_fields

end module following_scalar_field_generator
