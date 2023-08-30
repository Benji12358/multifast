
module following_arctan_transfo
    implicit none
    real*8      :: beta

contains

    subroutine t1_init_transfo(beta1)
        implicit none
        real*8      :: beta1

        beta=beta1

    end subroutine t1_init_transfo


    real*8 function t1_get_computational_coords(i, n)

        use following_mesh
        use boundaries
        implicit none

        real*8, intent(in)      :: i
        integer, intent(in)     :: n

        real*8  :: yeta

        yeta=(i-1)/dfloat(n-1)
        t1_get_computational_coords=yeta

    end function t1_get_computational_coords


    real*8 function t1_get_computational_coords_from_physical(y_physical)

        use following_mesh
        use boundaries
        implicit none

        real*8, intent(in)      :: y_physical
        real*8  :: yeta

        yeta=0.5d0 + (1.d0/beta) * datanh( dtanh(beta*0.5d0)*(y_physical-1.d0) )
        t1_get_computational_coords_from_physical=yeta

    end function t1_get_computational_coords_from_physical


    real*8 function t1_get_physical_coords(yeta)

        implicit none
        real*8  :: yeta

        t1_get_physical_coords=(dtanh(beta*(yeta-0.5d0))/dtanh(beta*0.5)+1.d0)

    end function t1_get_physical_coords


    real*8 function t1_DYonDy(s)

        implicit none

        real*8, intent(in)    ::s

        t1_DYonDy=dtanh(beta/2.d0)*(dcosh(beta*(s-0.5d0))**2)/beta

    end function t1_DYonDy

    real*8 function t1_D2YonDy2(s)
        implicit none
        real*8  :: A, dfds, dsdy

        real*8, intent(in)          ::s

        A=beta*(s-0.5d0)
        dfds=2.d0*dtanh(beta/2.d0)*dtanh(A)*(dcosh(A))**2
        dsdy=t1_DYonDy(s)

        t1_D2YonDy2=dfds*dsdy

    end function t1_D2YonDy2


end module following_arctan_transfo


module following_lamballais_transfo
    use mathematical_constants
    use decomp_2d
    implicit none
    integer     :: istret=2
    real*8      :: alpha, beta
    real*8      :: yly

contains


    subroutine t2_init_transfo(beta1)
        implicit none
        real*8      :: beta1
        real*8 :: yinf,den,xnum,xcx,den1,den2,den3,den4,xnum1,cst

        beta=beta1

        yly=2.d0


        yinf=-yly/2.
        den=2.*beta*yinf
        xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
        alpha=abs(xnum/den)

    end subroutine t2_init_transfo

    real*8 function t2_get_computational_coords(i, n)

        use following_mesh
        use boundaries
        implicit none

        real*8, intent(in)      :: i
        integer, intent(in)     :: n

        real*8  :: yeta, yp

        if (i==1.d0) then
            if (istret.eq.0) yeta=0.
            if (istret.eq.1) yeta=0.
            if (istret.eq.2) yeta=-1./2.
            if (istret.eq.3) yeta=-1./2.

        else
            if (istret==0) then
                yeta=(i-1.)*(1./n)
            endif
            if (istret==1) then
                if (BC2.eq.UNBOUNDED) yeta=(i-1.)*(1./n)
                if ((BC2.eq.FREESLIP).or.(BC2.eq.NOSLIP)) yeta=(i-1.)*(1./(n-1.))
            endif
            if (istret==2) then
                if (BC2.eq.UNBOUNDED) yeta=(i-1.)*(1./n)-0.5
                if ((BC2.eq.FREESLIP).or.(BC2.eq.NOSLIP)) yeta=(i-1.)*(1./(n-1.))-0.5
            endif
            if (istret==3) then
                if (BC2.eq.UNBOUNDED) yeta=((i-1.)*(1./2./n)-0.5)
                if ((BC2.eq.FREESLIP).or.(BC2.eq.NOSLIP)) yeta=((i-1.)*(1./2./(n-1.))-0.5)
            endif

        endif

        t2_get_computational_coords=yeta

    end function t2_get_computational_coords


    real*8 function t2_get_computational_coords_from_physical(y_physical)

        use following_mesh
        use boundaries
        implicit none

        real*8, intent(in)      :: y_physical
        real*8  :: yeta

        yeta=0.5d0 + (1.d0/beta) * datanh( dtanh(beta*0.5d0)*(y_physical-1.d0) )
        t2_get_computational_coords_from_physical=yeta
        write(*,*) 't2_get_computational_coords_from_physical for Lamballais transfo NOT IMPLEMENTED YET!'

    end function t2_get_computational_coords_from_physical


    real*8 function t2_get_physical_coords(yeta)

        use following_mesh
        use boundaries
        implicit none

        real*8  :: yeta, yp

        real*8 :: yinf,den,xnum,xcx,den1,den2,den3,den4,xnum1,cst

        yinf=-yly/2.
        den=2.*beta*yinf
        xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
        xcx=1./beta/alpha

        if (alpha.ne.0.) then

            if (yeta==-0.5d0) then

                if (istret.eq.1) yp=0.
                if (istret.eq.2) yp=0.
                if (istret.eq.3) yp=0.

            else

                den1=sqrt(alpha*beta+1.)
                xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
                den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
                den3=((sin(pi*yeta))*(sin(pi*yeta))/beta/pi)+alpha/pi
                den4=2.*alpha*beta-cos(2.*pi*yeta)+1.

                if ((BC2.ne.UNBOUNDED).and.(yeta==0.5d0).and.(istret==1)) then
                    xnum1=0.
                else
                    xnum1=(atan(xnum*tan(pi*yeta)))*den4/den1/den3/den
                endif

                cst=sqrt(beta)*pi/(2.*sqrt(alpha)*sqrt(alpha*beta+1.))

                if (istret==0) then
                    yp=yeta*2.d0
                endif

                if (istret==1) then
                    if (yeta.lt.0.5) yp=xnum1-cst-yinf
                    if (yeta.eq.0.5) yp=0.-yinf
                    if (yeta.gt.0.5) yp=xnum1+cst-yinf
                endif
                if (istret==2) then
                    if (yeta.lt.0.5) yp=xnum1-cst+yly
                    if (yeta.eq.0.5) yp=0.+yly
                    if (yeta.gt.0.5) yp=xnum1+cst+yly
                endif
                if (istret==3) then
                    if (yeta.lt.0.5) yp=(xnum1-cst+yly)*2.
                    if (yeta.eq.0.5) yp=(0.+yly)*2.
                    if (yeta.gt.0.5) yp=(xnum1+cst+yly)*2.
                endif
            end if


        endif

        t2_get_physical_coords=yp

    end function t2_get_physical_coords


    real*8 function t2_DYonDy(yeta)

        implicit none

        real*8, intent(in)    ::yeta

        if (istret==2) then
            t2_DYonDy=(alpha/pi+(1./pi/beta)*sin(pi*yeta)*sin(pi*yeta))
        end if

        if (istret==0) then
            t2_DYonDy=0.5d0
        end if

    end function t2_DYonDy

    real*8 function t2_D2YonDy2(yeta)
        implicit none

        real*8, intent(in)          ::yeta
        real*8  :: dfds, dsdy

        if (istret==2) then
            dfds=(2./beta*cos(pi*yeta)*sin(pi*yeta))
            dsdy=t2_DYonDy(yeta)

            t2_D2YonDy2=dfds*dsdy
            !t2_D2YonDy2=(2./beta*cos(pi*yeta)*sin(pi*yeta))/yly
        end if

        if (istret==0) then
            t2_D2YonDy2=0.d0
        end if

    end function t2_D2YonDy2


end module following_lamballais_transfo

module following_mesh_interface
    use following_lamballais_transfo
    use following_arctan_transfo

    procedure(t1_init_transfo), pointer                             :: init_transfo
    procedure(t1_get_computational_coords), pointer                 :: get_computational_coords
    procedure(t1_get_computational_coords_from_physical), pointer   :: get_computational_coords_from_physical
    procedure(t1_get_physical_coords), pointer                      :: get_physical_coords
    procedure(t1_DYonDy), pointer                                   :: DYonDy
    procedure(t1_D2YonDy2), pointer                                 :: D2YonDy2

end module following_mesh_interface

module following_mesh_loader

    use following_mesh_interface
    use following_lamballais_transfo
    use following_arctan_transfo
    use following_mesh

    contains

        subroutine configure_mesh_transfo(mesh_type)

        if (mesh_type==ORLANDI_MESH) then
            init_transfo                            =>t1_init_transfo
            get_computational_coords                =>t1_get_computational_coords
            get_computational_coords_from_physical  =>t1_get_computational_coords_from_physical
            get_physical_coords                     =>t1_get_physical_coords
            DYonDy                                  =>t1_DYonDy
            D2YonDy2                                =>t1_D2YonDy2
        end if

        if (mesh_type==LAMBALLAIS_MESH) then
            init_transfo                            =>t2_init_transfo
            get_computational_coords                =>t2_get_computational_coords
            get_computational_coords_from_physical  =>t2_get_computational_coords_from_physical
            get_physical_coords                     =>t2_get_physical_coords
            DYonDy                                  =>t2_DYonDy
            D2YonDy2                                =>t2_D2YonDy2
        end if

        end subroutine configure_mesh_transfo

end module following_mesh_loader

module following_mesh_generator

    ! Data modules
    use following_irregular_derivative_coefficients
    use following_mesh

    ! Processing modules
    use following_mesh_interface
    use following_mesh_loader

    use decomp_2d
    use mpi

    implicit none

contains

    ! Define all the mesh attributes from the mesh description given by arguments (n1,n2...)
    !   ni : number of point along the i th direction
    !   Li : Lenght of the computational domain in the i th direction
    !   str: The stretching applied along the Y direction
    subroutine generate_mesh(n1, n2, n3, L1, L2, L3, str)

        use workspace_view, only: log_path
        use decomp_2d

        use following_data

        implicit none
        ! Entry arguments
        integer, intent(in) :: n1, n2, n3
        real*8, intent(in)  :: L1, L2, L3

        ! Local variables
        integer                     :: i, j, k
        real*8                      :: str, dy, yeta
        integer, parameter          :: YMESH_UNIT=47

        integer                     :: writer_proc=0
        character(100)              :: mesh_generator_path, mesh_generator_path_y
        integer                     :: cell_center_file_id, cell_Yface_file_id, zone_id, var_id
        integer, dimension(n2m)     :: j_array
        real*8, dimension(n2m)      :: y_reg_c
        real*8, dimension(n2)       :: y_reg

        allocate(Z(n1))
        allocate(Zc(n1-1))

        allocate(Y(n2))
        allocate(Yc(n2-1))

        allocate(X(n3))
        allocate(Xc(n3-1))

        allocate(cell_size_Y(n2-1))

        yeta=0.d0

        dx2=1.d0/dfloat(n2m)    ! half height of a regular cell
        dx1=L1/dfloat(n1m)      ! spanwise size of a cell
        dx3=L3/dfloat(n3m)      ! streamwise size of a cell

        if (mesh_type/=NO_TRANSFO_MESH) then
            call configure_mesh_transfo(mesh_type)
            call init_transfo(str)
        endif

        ! Point positions in X, Z uniform directions -------------------
        do i=1,n1
            Z(i)=dfloat(i-1)*dx1
        enddo
        do i=1,n1-1
            Zc(i)=Z(i+1)-0.5d0*dx1
        enddo

        do k=1,n3
            X(k)=dfloat(k-1)*dx3
        enddo
        do k=1,n3-1
            Xc(k)=X(k+1)-0.5d0*dx3
        enddo

        if (mesh_type/=NO_TRANSFO_MESH) then

            ! Point positions in Y non-uniform directions ------------------------
            yeta=get_computational_coords(1.d0*n2, n2)
            Y(n2)=get_physical_coords(yeta)

            do j=1,n2-1

                yeta=get_computational_coords(1.d0*j, n2)
                Y(j)=get_physical_coords(yeta)

                yeta=get_computational_coords(j+0.5d0, n2)
                Yc(j)=get_physical_coords(yeta)

                if (j>1) cell_size_Y(j-1)=Y(j)-Y(j-1)
            enddo

            cell_size_Y(n2-1)=Y(n2)-Y(n2-1)

            !open(1542, file="Y_field", position="append")
            !do i = 1,n1
             !   do k = 1,n3
              !      do j = 1,n2
               !         Y_field(i,j,k) = Y(j)
                !        write(1542,*)Y(j)
                 !   end do
                !end do
            !end do
            !close(1542)

            ! *************************************************************************
            ! Defining coefficient for schemes in irregular directions-----------------
            ! *************************************************************************

            ! Y direction--------------------------------------------------------------
            allocate(Y_to_YTr_for_D1(n2))
            Y_to_YTr_for_D1=0.d0
            allocate(Yc_to_YcTr_for_D1(n2-1))
            Yc_to_YcTr_for_D1=0.d0

            allocate(Y_to_YTr_for_D2(n2, 2))
            Y_to_YTr_for_D2=0.d0
            allocate(Yc_to_YcTr_for_D2(n2-1, 2))
            Yc_to_YcTr_for_D2=0.d0

            ! 1st derivative
            do j = 1, n2
                Y_to_YTr_for_D1(j)=DYonDy(get_computational_coords(j*1.d0, n2))
            end do

            do j = 1, n2m
                Yc_to_YcTr_for_D1(j)=DYonDy(get_computational_coords(j+0.5d0, n2))
            end do

                ! 2nd derivative

            do j = 2, n2m
                Y_to_YTr_for_D2(j,1)=D2YonDy2(get_computational_coords(j*1.d0, n2))
                Y_to_YTr_for_D2(j,2)=Y_to_YTr_for_D1(j)**2
            end do

            do j = 1, n2m
                Yc_to_YcTr_for_D2(j,1)=D2YonDy2(get_computational_coords(j+0.5d0, n2))
                Yc_to_YcTr_for_D2(j,2)=Yc_to_YcTr_for_D1(j)**2
            end do

        else

            ! Point positions in Y uniform directions -------------------
            ! /!\ Caution, here we use half height of the cells in y-direction
            do j=1,n2
                ! Y(j)=dfloat(j-1)*dx2
                Y(j)=dfloat(j-1)*2.d0*dx2
            enddo
            do j=1,n2-1
                ! Yc(j)=Y(j+1)-0.5d0*dx2
                Yc(j)=Y(j+1)-dx2

                if (j>1) cell_size_Y(j-1)=Y(j)-Y(j-1)
            enddo

            allocate(Y_to_YTr_for_D1(n2))
            Y_to_YTr_for_D1=1.d0
            allocate(Yc_to_YcTr_for_D1(n2-1))
            Yc_to_YcTr_for_D1=1.d0

            allocate(Y_to_YTr_for_D2(n2, 2))
            Y_to_YTr_for_D2(:,1)=0.d0
            Y_to_YTr_for_D2(:,2)=1.d0
            allocate(Yc_to_YcTr_for_D2(n2-1, 2))
            Yc_to_YcTr_for_D2(:,1)=0.d0
            Yc_to_YcTr_for_D2(:,2)=1.d0

        endif


        ! Write mesh in Log path ____________________________________________________

        if (nrank==0) then

            mesh_generator_path='Log/Mesh_generator/following/'

            mesh_generator_path=trim(log_path)//"Mesh_generator/"
            mesh_generator_path_y=trim(mesh_generator_path)//"following/"

            open(YMESH_UNIT,file=trim(mesh_generator_path_y)//'ymesh.out')

            write(YMESH_UNIT,*) 'j, y, dy, yc, a, b, c, a_c, b_c, c_c'
            do j=1,n2m
                dy=Y(j+1)-Y(j)
                write(YMESH_UNIT,*) j, Y(j), dy, Yc(j), Y_to_YTr_for_D1(j), Y_to_YTr_for_D2(j,2), Y_to_YTr_for_D2(j,1), Yc_to_YcTr_for_D1(j), Yc_to_YcTr_for_D2(j,2), Yc_to_YcTr_for_D2(j,1)
            end do

            close(YMESH_UNIT)

        end if


    end subroutine generate_mesh

end module following_mesh_generator
