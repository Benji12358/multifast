
module arctan_transfo
    implicit none
    real*8      :: beta

contains

    subroutine t1_init_transfo(beta1)
        implicit none
        real*8      :: beta1

        beta=beta1
    end subroutine t1_init_transfo


    real*8 function t1_get_computational_coords(i, n)

        use mesh
        use boundaries
        implicit none

        real*8, intent(in)      :: i
        integer, intent(in)     :: n

        real*8  :: yeta

        yeta=(i-1)/dfloat(n-1)
        t1_get_computational_coords=yeta

    end function t1_get_computational_coords


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


end module arctan_transfo


module lamballais_transfo
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

        use mesh
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


    real*8 function t2_get_physical_coords(yeta)

        use mesh
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


end module lamballais_transfo

module transfo_tester
    implicit none

contains

    subroutine test_transfo(y, a, b, c, dy)
        use mesh, only: n2
        use d1c_schemes
        use d2c_schemes
        use workspace_view, only: log_path

        implicit none

        real*8, dimension(n2)   :: fy, dfy, ddfy, dfy_th, ddfy_th
        real*8, dimension(:)    :: y, a, b, c
        integer                 :: j
        real*8                  :: dy
        integer                 :: TEST_UNIT=152
        character(100)          :: mesh_generator_path, mesh_generator_path_y

        fy=0.d0
        dfy=0.d0
        ddfy=0.d0
        dfy_th=0.d0
        ddfy_th=0.d0

        mesh_generator_path='Log/Mesh_generator/Y/'

        mesh_generator_path=trim(log_path)//"Mesh_generator/"
        mesh_generator_path_y=trim(mesh_generator_path)//"Y/"

        do j = 1, n2
            fy(j)=3*y(j)**4 + 5*y(j)**3+ 2*y(j)**2
            dfy_th(j)=12*y(j)**3+15*y(j)**2+4*y(j)
            ddfy_th(j)=36*y(j)**2+30*y(j)+4
        end do


        call D1_ExpCtr_O0Fp7_MULT(fy(:), dfy(:), n2, dy, .false., Dirichlet, a)

        open(TEST_UNIT,file=trim(mesh_generator_path_y)//'test_d1.csv')

        do j = 1, n2-1
            write(TEST_UNIT,*)j, dfy(j), dfy_th(j)
        end do

        close(TEST_UNIT)

        call D1_ExpCtr_O0Fp7_MULT(fy(:), ddfy(:), n2, dy, .false., Dirichlet, b)
        call D2_ExpCtr_O0Fp6_MULT_ACC(fy(:), ddfy(:), n2, dy, .false., Dirichlet, c)

        open(TEST_UNIT,file=trim(mesh_generator_path_y)//'test_d2.csv')

        do j = 1, n2-1
            write(TEST_UNIT,*)j, ddfy(j), ddfy_th(j)
        end do

        close(TEST_UNIT)


    end subroutine test_transfo

end module transfo_tester

module mesh_generator

    ! Data modules
    use irregular_derivative_coefficients
    use mesh

    ! Processing modules
    use lamballais_transfo
    use arctan_transfo

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
        use transfo_tester

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


        procedure(t1_init_transfo), pointer             :: init_transfo
        procedure(t1_get_computational_coords), pointer :: get_computational_coords
        procedure(t1_get_physical_coords), pointer      :: get_physical_coords
        procedure(t1_DYonDy), pointer                   :: DYonDy
        procedure(t1_D2YonDy2), pointer                 :: D2YonDy2

        if (mesh_type==ORLANDI_MESH) then
            init_transfo            =>t1_init_transfo
            get_computational_coords=>t1_get_computational_coords
            get_physical_coords     =>t1_get_physical_coords
            DYonDy                  =>t1_DYonDy
            D2YonDy2                =>t1_D2YonDy2
        end if

        if (mesh_type==LAMBALLAIS_MESH) then
            init_transfo            =>t2_init_transfo
            get_computational_coords=>t2_get_computational_coords
            get_physical_coords     =>t2_get_physical_coords
            DYonDy                  =>t2_DYonDy
            D2YonDy2                =>t2_D2YonDy2
        end if

        allocate(Z(n1))
        allocate(Zc(n1-1))

        allocate(Y(n2))
        allocate(Yc(n2-1))

        allocate(X(n3))
        allocate(Xc(n3-1))

        allocate(cell_size_Y(n2-1))

        allocate(Y_field(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

        yeta=0.d0

        dx2=1.d0/dfloat(n2m)    ! half height of a regular cell
        dx1=L1/dfloat(n1m)      ! spanwise size of a cell
        dx3=L3/dfloat(n3m)      ! streamwise size of a cell

        call init_transfo(str)

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

        do k = zstart(3), zend(3)
                do j = zstart(2), zend(2)
                    do i = zstart(1), zend(1)
                        Y_field(i,j,k) = Y(j)
                    end do
                end do
            end do

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

        ! Write mesh in Log path ____________________________________________________

        if (nrank==0) then

            mesh_generator_path='Log/Mesh_generator/Y/'

            mesh_generator_path=trim(log_path)//"Mesh_generator/"
            mesh_generator_path_y=trim(mesh_generator_path)//"Y/"

            open(YMESH_UNIT,file=trim(mesh_generator_path_y)//'ymesh.out')

            write(YMESH_UNIT,*) 'j, y, dy, yc, a, b, c'
            do j=1,n2m
                dy=Y(j+1)-Y(j)
                write(YMESH_UNIT,*) j, Y(j), dy, Yc(j), Y_to_YTr_for_D1(j), Y_to_YTr_for_D2(j,2), Y_to_YTr_for_D2(j,1)
            end do

            close(YMESH_UNIT)

        end if

        if (nrank==0) call test_transfo(Y, Y_to_YTr_for_D1, Y_to_YTr_for_D2(:,1), Y_to_YTr_for_D2(:,2), dx2)


    end subroutine generate_mesh

end module mesh_generator
