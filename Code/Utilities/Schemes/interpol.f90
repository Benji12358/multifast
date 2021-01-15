module interpol

    implicit none

contains

    subroutine interpolate_point(pt, x1, x2, x3, field, value, error_flag)

        implicit none

        real*8, dimension(3), intent(in)        :: pt
        real*8, dimension(:), intent(in)        :: x1, x2, x3
        real*8, dimension(:,:,:), intent(in)    :: field
        integer                    :: n1, n2, n3
        real*8, intent(out)                     :: value
        integer,intent(out), optional           :: error_flag

        integer                                 :: i, x0, y0, z0
        real*8                                  :: c0, c1, c00, c01, c11, c10
        real*8                                  :: xd, yd, zd

        n1=size (field ,1)
        n2=size (field ,2)
        n3=size (field ,3)

        x0=0
        y0=0
        z0=0

        if (present(error_flag)) error_flag=0

        ! Getting the index of origins from all directions, in order to build a 3D cell around the point pt
        do i = 1, n1-1
            if ((pt(1)>=x1(i)).and.(pt(1)<=x1(i+1))) then
                x0=i
            end if
        end do

        do i = 1, n2-1
            if ((pt(2)>=x2(i)).and.(pt(2)<=x2(i+1))) then
                y0=i
            end if
        end do

        do i = 1, n3-1
            if ((pt(3)>=x3(i)).and.(pt(3)<=x3(i+1))) then
                z0=i
            end if
        end do

        if ((x0==0).or.(y0==0).or.(z0==0)) then
            value=0.d0
            if (present(error_flag)) error_flag=1
            return
        endif

        ! Trilinear interpolation
        xd = ( pt(1) - x1(x0) ) / ( x1(x0+1) - x1(x0) )
        yd = ( pt(2) - x2(y0) ) / ( x2(y0+1) - x2(y0) )
        zd = ( pt(3) - x3(z0) ) / ( x3(z0+1) - x3(z0) )

        c00 = field(x0,y0,z0) * (1-xd) + field(x0+1,y0,z0) * xd

        c10 = field(x0,y0+1,z0) * (1-xd) + field(x0+1,y0+1,z0) * xd
        c01 = field(x0,y0,z0+1) * (1-xd) + field(x0+1,y0,z0+1) * xd
        c11 = field(x0,y0+1,z0+1) * (1-xd) + field(x0+1,y0+1,z0+1) * xd

        c0 = c00 * (1-yd) + c10 * yd
        c1 = c01 * (1-yd) + c11 * yd

        value = c0 * (1-zd) + c1 * zd


    end subroutine interpolate_point

    subroutine interpolate_point2(pt, x1, x2, x3, field, value, error_flag, funit)

        implicit none

        real*8, dimension(3), intent(in)        :: pt
        real*8, dimension(:), intent(in)        :: x1, x2, x3
        real*8, dimension(:,:,:), intent(in)    :: field
        integer                    :: n1, n2, n3
        real*8, intent(out)                     :: value
        integer,intent(out), optional           :: error_flag

        integer                                 :: i, x0, y0, z0
        real*8                                  :: c0, c1, c00, c01, c11, c10
        real*8                                  :: xd, yd, zd
        integer                                 :: funit

        n1=size (field ,1)
        n2=size (field ,2)
        n3=size (field ,3)

        x0=0
        y0=0
        z0=0

        if (present(error_flag)) error_flag=0

        ! Getting the index of origins from all directions, in order to build a 3D cell around the point pt
        do i = 1, n1-1
            if ((pt(1)>=x1(i)).and.(pt(1)<=x1(i+1))) then
                x0=i
            end if
        end do

        do i = 1, n2-1
            if ((pt(2)>=x2(i)).and.(pt(2)<=x2(i+1))) then
                y0=i
            end if
        end do

        do i = 1, n3-1
            if ((pt(3)>=x3(i)).and.(pt(3)<=x3(i+1))) then
                z0=i
            end if
        end do

        if ((x0==0).or.(y0==0).or.(z0==0)) then
            value=0.d0
            if (present(error_flag)) error_flag=1
            return
        endif

        write(funit,*)'interpolate_point2: A', x0,y0,z0

        ! Trilinear interpolation
        xd = ( pt(1) - x1(x0) ) / ( x1(x0+1) - x1(x0) )
        yd = ( pt(2) - x2(y0) ) / ( x2(y0+1) - x2(y0) )
        zd = ( pt(3) - x3(z0) ) / ( x3(z0+1) - x3(z0) )

        write(funit,*)'interpolate_point2: B', xd,yd,zd

        write(funit,*)'interpolate_point3: C', x1(x0), x2(y0), x3(z0)

        write(funit,*)'interpolate_point3: D000', field(x0,y0,z0)
        write(funit,*)'interpolate_point3: D100', field(x0+1,y0,z0)
        write(funit,*)'interpolate_point3: D010', field(x0,y0+1,z0)
        write(funit,*)'interpolate_point3: D110', field(x0+1,y0+1,z0)
        write(funit,*)'interpolate_point3: D001', field(x0,y0,z0+1)
        write(funit,*)'interpolate_point3: D101', field(x0+1,y0,z0+1)
        write(funit,*)'interpolate_point3: D011', field(x0,y0+1,z0+1)
        write(funit,*)'interpolate_point3: D111', field(x0+1,y0+1,z0+1)

        c00 = field(x0,y0,z0) * (1-xd) + field(x0+1,y0,z0) * xd

        c10 = field(x0,y0+1,z0) * (1-xd) + field(x0+1,y0+1,z0) * xd
        c01 = field(x0,y0,z0+1) * (1-xd) + field(x0+1,y0,z0+1) * xd
        c11 = field(x0,y0+1,z0+1) * (1-xd) + field(x0+1,y0+1,z0+1) * xd

        c0 = c00 * (1-yd) + c10 * yd
        c1 = c01 * (1-yd) + c11 * yd

        value = c0 * (1-zd) + c1 * zd


    end subroutine interpolate_point2

end module interpol
