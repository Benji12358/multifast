module SCALAR_field_generator
    use decomp_2d
    use mesh

    use SCALAR_data

    implicit none

contains

    subroutine generate_fields(sca_x, sca_y, sca_z, sca_down, sca_up)
        implicit none
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: sca_x
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))  :: sca_y
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))  :: sca_z
        real*8                                                                      :: sca_down, sca_up

        real*8  :: delta

        integer i, j, k

        delta=sca_up-sca_down

        sca_y=0.d0

        do k = ystart(3), min(n3-1, yend(3))
            do i = ystart(1), min(n1-1, yend(1))
                do j = 1, n2-1
                    sca_y(i,j,k)=sca_down+(Yc(j)/L2)*delta
                    !sca_y(i,j,k)=0.d0
                end do
            end do
        end do

		if(nrank==0) write(*,*)'SCALAR_generate_fields: Tinf Tsup', sca_down, sca_up, delta

        sca_wall10=0.d0
        sca_wall11=0.d0

        sca_wall20=sca_down
        sca_wall21=sca_up

        sca_wall30=0.d0
        sca_wall31=0.d0

        call transpose_y_to_x(sca_y, sca_x)
        call transpose_y_to_z(sca_y, sca_z)

    end subroutine generate_fields

end module SCALAR_field_generator
