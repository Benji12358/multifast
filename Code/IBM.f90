module object_drawer
    implicit none

    public :: perform_mask
    private

contains



    subroutine sort_array(T, nb, verbose)
        implicit none
        real*8, dimension(nb)    :: T
        integer :: i,j, nb
        integer :: imin
        real*8  :: tmp
        integer, optional :: verbose


        if (present(verbose)) write(*,*)'minval(T)', minval(T)
        do i = 1, nb
            imin=-1
            do j = i, nb
                if (T(j)==minval(T(i:nb))) then
                    imin=j
                endif

            end do

            if (present(verbose)) write(*,*)'#####################'
            if (present(verbose)) write(*,*)'i=', i
            if (present(verbose)) write(*,*)'Minimum en imin=', imin
            if (present(verbose)) write(*,*)'T1', T

            if(imin>0) then
                tmp=T(i)
                T(i)=T(imin)
                T(imin)=tmp
            endif

            if (present(verbose)) write(*,*)'T2', T

        end do

    end subroutine sort_array

    subroutine get_objet_area(vertex, nb_vertex, xs, xe, ys, ye, zs, ze)
        implicit none
        integer                                     :: nb_vertex
        real*8, dimension(nb_vertex, 3)             :: vertex
        real*8                                      :: xs, xe, ys, ye, zs, ze

        integer :: i
        real*8  :: sommet(3)

        xs=10000.d0; ys=10000.d0; zs=10000.d0
        xe=0.d0; ye=0.d0; ze=0.d0

        do i = 1, nb_vertex
            sommet=vertex(i,:)
            xs=min(xs, sommet(1))
            ys=min(ys, sommet(2))
            zs=min(zs, sommet(3))

            xe=max(xe, sommet(1))
            ye=max(ye, sommet(2))
            ze=max(ze, sommet(3))
        end do

    end subroutine get_objet_area

    subroutine inTriangle(A1, B1, C1, P, isInTriangle, tolerance, verbose)
        implicit none
        real*8, dimension(3)        :: P, A1, B1, C1
        logical                     :: isInTriangle
        real*8                      :: tolerance
        integer, optional           :: verbose
        real*8  :: sum_angle, maxangle
        real*8, dimension(3)        :: vec01, vec02, vec03
        real*8                          :: vec01_norm, vec02_norm, vec03_norm
        real*8                          :: angle12, angle13, angle23
        real*8                          :: cosangle12, cosangle13, cosangle23
        real*8                          :: tmp
        real*8                          :: prod12, prod13, prod23

        real*8  :: normAP, normBP, normCP
        real*8 , dimension(3)   :: AP, BP, CP, AB, BC, AC
        real*8, dimension(3)    :: A, B, C
        real*8                  :: ca, cb, cc, cd,ce,cf, det1, det2, det3, alpha, beta
        logical                     :: isInTriangle1,isInTriangle2, isInTriangle3


        real*8 , dimension(3)   :: AP_sur_AB, AP_sur_AC, P_dans_triangle
        real*8                      :: tolerance2=1.d-14

        isInTriangle=.true.
        isInTriangle1=.true.
        isInTriangle2=.true.
        isInTriangle3=.true.

        if (present(verbose)) write(*,*)'TRIANGLE3'

        A=A1
        B=B1
        C=C1

        AP=P-A
        normAP=dsqrt(AP(1)**2+AP(2)**2+AP(3)**2)
        BP=P-B
        normBP=dsqrt(BP(1)**2+BP(2)**2+BP(3)**2)
        CP=P-C
        normCP=dsqrt(CP(1)**2+CP(2)**2+CP(3)**2)

        AB=B-A
        BC=C-B
        AC=C-A

        if(normAP<min(normBP,normCP)) then
            if (present(verbose)) write(*,*)'ROTATION'
            if (present(verbose)) write(*,*)'A', A
            if (present(verbose)) write(*,*)'B', B
            if (present(verbose)) write(*,*)'C', C

            A=B1
            B=C1
            C=A1

            AP=P-A
            BP=P-B
            CP=P-C

            AB=B-A
            BC=C-B
            AC=C-A

            if (present(verbose)) write(*,*)'A', A
            if (present(verbose)) write(*,*)'B', B
            if (present(verbose)) write(*,*)'C', C
        endif



        if (.true.) then !if(normCP>max(normAP,normBP)) then
            ca=BC(1)
            cb=-AP(1)
            ce=-AB(1)
            cc=BC(2)
            cd=-AP(2)
            cf=-AB(2)

            det1=ca*cd-cb*cc

            alpha=(ce*cd-cb*cf)/det1
            beta=(ca*cf-ce*cc)/det1


            if (((alpha>1.d0).or.(alpha<0.d0).or.(beta<1.d0)).and.(abs(det1)>tolerance2)) then
                isInTriangle1=.false.
                if (present(verbose)) write(*,*)'POINT1'
                if (present(verbose)) write(*,*)'alpha', alpha
                if (present(verbose)) write(*,*)'beta', beta
                if (present(verbose)) write(*,*)'det', det1

            else if (abs(det1)>tolerance2) then
                if (present(verbose)) write(*,*)'POINT1'
                if (present(verbose)) write(*,*)'alpha', alpha
                if (present(verbose)) write(*,*)'beta', beta
            endif

            if (abs(det1)<tolerance2) then

                ! VOIR SI EQ(1) et EQ(2) sont incompatibles ou identiques
                ! Si incompatible isInTriangle1=.false.

                if (present(verbose)) write(*,*)'POINT1 NE SAIT PAS', det1
                if(abs(ce-ca*cf/cc)<tolerance) then
                    if (present(verbose))write(*,*)'EQUATIONS 1 et 2 IDENTIQUES', det1
                endif
                if(abs(ce-ca*cf/cc)>tolerance) then
                    if (present(verbose)) write(*,*)'EQUATIONS 1 et 2 INCOMPATIBLES', det1
                    isInTriangle1=.false.
                endif

            endif

            ! Si Eq(1) compatible avec Eq(2) on verifie Eq(1) avec Eq(3).
            ! Si Eq(1) et Eq(2) identiques,  on utilise Eq(3)


            if (isInTriangle1) then

                cc=BC(3)
                cd=-AP(3)
                cf=-AB(3)
                det2=ca*cd-cb*cc

                alpha=(ce*cd-cb*cf)/det2
                beta=(ca*cf-ce*cc)/det2

                if (((alpha>1.d0).or.(alpha<0.d0).or.(beta<1.d0)).and.(abs(det2)>tolerance2)) then
                    isInTriangle2=.false.
                    if (present(verbose)) write(*,*)'POINT2'
                    if (present(verbose)) write(*,*)'alpha', alpha
                    if (present(verbose)) write(*,*)'beta', beta


                else if (abs(det2)>tolerance2) then
                    if (present(verbose)) write(*,*)'POINT2'
                    if (present(verbose)) write(*,*)'alpha', alpha
                    if (present(verbose)) write(*,*)'beta', beta

                endif

                if (abs(det2)<tolerance2) then

                    ! VOIR SI EQ(1) et EQ(3) sont incompatibles ou identiques
                    ! Si incompatible isInTriangle2=.false.

                    if (present(verbose)) write(*,*)'POINT2 NE SAIT PAS', det2
                    if(abs(ce-ca*cf/cc)<tolerance) then
                        if (present(verbose))write(*,*)'EQUATIONS 1 et 3 IDENTIQUES car'
                        if (present(verbose))write(*,*)abs(ce-ca*cf/cc),'<',tolerance
                    endif
                    if(abs(ce-ca*cf/cc)>tolerance) then
                        if (present(verbose)) write(*,*)'EQUATIONS 1 et 3 INCOMPATIBLES'
                        if (present(verbose))write(*,*)abs(ce-ca*cf/cc),'>',tolerance
                        isInTriangle2=.false.
                    endif

                endif



                if (isInTriangle2) then
                    ca=BC(2)
                    cb=-AP(2)
                    ce=-AB(2)
                    det3=ca*cd-cb*cc

                    alpha=(ce*cd-cb*cf)/det3
                    beta=(ca*cf-ce*cc)/det3

                    if (((alpha>1.d0).or.(alpha<0.d0).or.(beta<1.d0)).and.(abs(det3)>tolerance2)) then
                        isInTriangle3=.false.
                        if (present(verbose))  write(*,*)'POINT3', alpha, beta

                    else if (abs(det3)>tolerance2) then
                        if (present(verbose)) write(*,*)'POINT3'
                        if (present(verbose)) write(*,*)'alpha', alpha
                        if (present(verbose)) write(*,*)'beta', beta
                    endif

                    if (abs(det3)<tolerance2) then

                            ! VOIR SI EQ(1) et EQ(3) sont incompatibles ou identiques
                            ! Si incompatible isInTriangle2=.false.

                        if (present(verbose)) write(*,*)'POINT2 NE SAIT PAS', det3
                        if(abs(ce-ca*cf/cc)<tolerance) then
                            if (present(verbose))write(*,*)'EQUATIONS 2 et 3 IDENTIQUES'
                        endif
                        if(abs(ce-ca*cf/cc)>tolerance) then
                            if (present(verbose)) write(*,*)'EQUATIONS 2 et 3 INCOMPATIBLES'
                            isInTriangle3=.false.
                        endif
                    endif

                endif




            end if

            isInTriangle=isInTriangle1.and.isInTriangle2.and.isInTriangle3

            if (present(verbose)) then

                AP_sur_AB=AB(1)*AP(1)+AB(2)*AP(2)+AB(3)*AP(3)
                AP_sur_AC=AC(1)*AP(1)+AC(2)*AP(2)+AC(3)*AP(3)
                P_dans_triangle=A+(AP_sur_AB*AB+AP_sur_AC*AC)

                write(*,*)
                write(*,*)
                write(*,*)'POINT P               : ', P
                write(*,*)'POINT DANS TRIANGLE P : ', P_dans_triangle
                write(*,*) 'TRIANGLE:'
                write(*,*)'A', A1
                write(*,*)'B', B1
                write(*,*)'C', C1
                if (isInTriangle) write(*,*)'Le point EST dans le triangle'
                if (.not. isInTriangle) write(*,*)'Le point NEST PAS dans le triangle'
            endif


        endif





    end subroutine inTriangle

    subroutine find_equation_x(face1, face2, face3, a1, a2, a3, det,verbose)
        implicit none
        real*8, dimension(3)    :: face1, face2, face3
        real*8                      :: a1,a2,a3
        real*8                      :: a,b,c,d,e,f
        real*8                      :: x1,x2,x3,y1,y2,y3,z1,z2,z3
        integer                     :: orientation
        real*8                      :: det
        integer, optional           :: verbose

        ! eq: y2-y1= a*(x2-x1)+ b*(z2-z1)

        if (present(verbose)) write(*,*) "Surface x"

        x1=face1(1)
        x2=face2(1)
        x3=face3(1)

        y1=face1(2)
        y2=face2(2)
        y3=face3(2)

        z1=face1(3)
        z2=face2(3)
        z3=face3(3)

        det=(y2-y1)*(z3-z1)-(z2-z1)*(y3-y1)

        if (det.ne.(0.d0)) then

            a=(y2-y1)
            b=(z2-z1)
            e=(x2-x1)

            c=(y3-y1)
            d=(z3-z1)
            f=(x3-x1)

            a1=(e*d-b*f)/det
            a2=(a*f-e*c)/det
            a3=x1-a1*y1-a2*z1
            if (present(verbose)) then
                write(*,*)
                write(*,*) "a", a1, e,d,b,f
                write(*,*) "b", a2, a,f,e,c
                write(*,*) "c", a3

                write(*,*)"___________"
            endif

        else
            if (present(verbose)) write(*,*)'Surface oriente selon y ou z'
        endif

    end subroutine find_equation_x

    subroutine intersection(face1, face2, face3, y, z, x, error, verbose)
        implicit none

        real*8, dimension(3)      :: face1, face2, face3
        real*8                      :: y,z, x
        real*8                      :: det, a,b,c
        integer                     :: error
        integer, optional           :: verbose


        call find_equation_x(face1, face2, face3, a,b,c, det)
        error=0

        if (det==0) then
            error=1
        else
            if (present(verbose)) write(*,*)'a,b,c', a,b,c
            x=a*y+b*z+c
        endif

    end subroutine intersection

    subroutine perform_mask(vertex, faces, nb_vertex, nb_faces, mask_field, n1, n2, n2s,n2e,n3, n3s,n3e, X, Y, Z, value, verbose)
        implicit none
        integer                                     :: n1, n2, n2s,n2e,n3, n3s,n3e, nb_vertex, nb_faces
        real*8, dimension(n1, n2s:n2e,n3s:n3e)      :: mask_field
        integer, dimension(nb_faces, 3)             :: faces
        real*8, dimension(nb_vertex, 3)             :: vertex
        real*8, dimension(:)                        :: X, Y, Z
        real*8                                      :: value
        integer, optional                           :: verbose

        real*8, dimension(3)                        :: face1, face2, face3, point
        logical :: isInFace=.false.

        real*8  :: x_intersec(40), x_intersec_clean(40), x_i
        integer :: face_intersec(40)
        integer :: nb_xintersec
        integer :: pc

        integer     :: i,j,k,f

        real*8  :: xs, xe, ys, ye, zs, ze
        integer :: i0,i1, j0,j1, k0,k1, ierr
        real*8  :: tolerance=1.d-13



        mask_field=0.d0

        call get_objet_area(vertex, nb_vertex, xs, xe, ys, ye, zs, ze)

        i0=1; i1=n1; j0=1; j1=n2; k0=1; k1=n3;

        do i = 1, n1
            if (X(i)>=xs) then
                i0=i
                exit
            endif
        end do

        do i = n1, 1, -1
            if (X(i)<=xe) then
                i1=i
                exit
            endif
        end do

        do j = 1, n2
            if (Y(j)>=ys) then
                j0=j
                exit
            endif
        end do

        do j = n2, 1, -1
            if (Y(j)<=ye) then
                j1=j
                exit
            endif
        end do

        do k = 1, n3
            if (Z(k)>=zs) then
                k0=k
                exit
            endif
        end do

        do k = n3, 1, -1
            if (Z(k)<=ze) then
                k1=k
                exit
            endif
        end do

        write(*,*)'objet ext:', xs, xe, ys, ye, zs, ze
        write(*,*)'i', i0,':',i1
        write(*,*)'j', j0,':',j1
        write(*,*)'k', k0,':',k1
        write(*,*)'value', value
        call sleep(0)

        ! is=(xs-Ox)/(Lx-Ox)

        do j = n2s, n2e
            do k = n3s, n3e

                nb_xintersec=0
                x_intersec=0.d0

                do f = 1, nb_faces
                    face1=vertex(faces(f,1), :)
                    face2=vertex(faces(f,2), :)
                    face3=vertex(faces(f,3), :)
                    call intersection(face1, face2, face3, Y(j),Z(k), x_i, ierr)

                    if (ierr==0) then

                        point=(/x_i,Y(j),Z(k)/)

                        call inTriangle(face1, face2, face3, point, isInFace, tolerance)

                        if (isInFace) then

                            call inTriangle(face1, face2, face3, point, isInFace, tolerance)
                            nb_xintersec=nb_xintersec+1
                            x_intersec(nb_xintersec)=x_i
                            face_intersec(nb_xintersec)=f

                            if (present(verbose)) then
                                write(*,*)
                                write(*,*)'Face', f
                                write(*,*)'sommet1', face1
                                write(*,*)'sommet2', face2
                                write(*,*)'sommet3', face3

                                write(*,*) 'Ligne y=', Y(j), 'z=', Z(k)
                                write(*,*)'Point dintersection', x_i
                                write(*,*)
                            endif

                        endif

                    endif
                end do

                if(present(verbose)) write(*,*)'j,k', j,k
                if(present(verbose)) write(*,*) 'x_intersec', x_intersec(1:nb_xintersec)
                call sort_array(x_intersec, nb_xintersec)

                if(mod(nb_xintersec,2)==1) then

                    if(present(verbose)) then
                        write(*,*)
                        write(*,*)
                        write(*,*) 'Nombre dintersection impair ', nb_xintersec
                        write(*,*) 'j,k', j,k
                        write(*,*)'Y',Y(j)
                        write(*,*)'Z',Z(k)

                        write(*,*) 'Intersections en:', x_intersec(1:nb_xintersec)

                        write(*,*)
                        write(*,*) 'Recherche de doublons...'
                    endif

                    x_intersec_clean(1)=x_intersec(1)
                    pc=1
                    do i = 2, nb_xintersec
                        if(abs(x_intersec(i)-x_intersec_clean(pc))>1.d-9) then
                            pc=pc+1
                            x_intersec_clean(pc)=x_intersec(i)
                        endif
                    end do
                    nb_xintersec=pc
                    x_intersec=0.d0
                    x_intersec(1:nb_xintersec)=x_intersec_clean(1:pc)

                    if(present(verbose)) write(*,*)
                    if(present(verbose)) write(*,*) 'Tableaux des intersections nettoye'
                    if(present(verbose)) write(*,*) 'Intersections en:', x_intersec(1:nb_xintersec)
                    if(present(verbose)) call sleep(0)
                endif

                if(mod(nb_xintersec,2)==1) then
                    write(*,*) 'ERROR!!!', j,k, nb_xintersec
                    write(*,*) 'x_intersec', x_intersec(1:nb_xintersec)
                    write(*,*)'face_intersec', face_intersec(1:nb_xintersec)
                endif

                pc=1

                if ((nb_xintersec>0).and.(mod(nb_xintersec,2)==0)) then

                    if(present(verbose)) write(*,*) 'x_intersec', x_intersec(1:nb_xintersec)

                    do i = 1, n1
                        if ((X(i)>x_intersec(pc)).and.(X(i)<x_intersec(pc+1))) then
                            mask_field(i,j,k)=value

                            if (i>10000) then   ! POUR DEBUGGER CHANGER LA VALEUR DU SEUIL
                                write(*,*) 'i',i
                                write(*,*) 'pc',pc, '/', nb_xintersec
                                write(*,*) 'x_intersec', x_intersec
                                write(*,*) 'j,k', j,k
                                call sleep(2)
                            endif

                        else if((X(i)>x_intersec(pc)).and.(X(i)>x_intersec(pc+1))) then
                            if (nb_xintersec>=pc+3) pc=pc+2
                        else
                            mask_field(i,j,k)=0.d0
                        endif
                    end do

                    if(present(verbose)) write(*,*)'+++++++++++++++++++++++++++++++++++++++'

                endif

            end do
        end do

    end subroutine perform_mask

end module object_drawer


module IBM
    use IBM_data
    use IBM_settings
    use decomp_2d
    implicit none

contains

    subroutine IBM_setup()
        use snapshot_writer
        use object_drawer
        use object_file_reader
        use mesh
        use COMMON_workspace_view
        implicit none

        real*8, dimension(:,:), allocatable                 :: vertex
        integer, dimension(:,:), allocatable                :: faces
        integer                                             :: nb_faces, nb_vertex
        real*8  :: X1(n1), X2(n2), X3(n3)
        real*8  :: X1c(n1), X2c(n2), X3c(n3)
        integer     :: i,j,k
        character(200)                  :: snapshot_file

        ! In this context, the mesh array must be n1,n2 or n3-array even for the staggered ones
        ! (because of subroutine perform_mask)
        do i = 1, n1
            X1(i)=Z(i)
        end do

        do j = 1, n2
            X2(j)=Y(j)
        end do

        do k = 1, n3
            X3(k)=X(k)
        end do

        do i = 1, n1-1
            X1c(i)=Zc(i)
        end do

        do j = 1, n2-1
            X2c(j)=Yc(j)
        end do

        do k = 1, n3-1
            X3c(k)=Xc(k)
        end do

        X1c(n1)=X1c(n1-1)+dx1
        X2c(n2)=X2c(n2-1)+Y(n2)-Y(n2-1)
        X3c(n3)=X3c(n3-1)+dx3


        ! READING OBJECT FILE
        call read_object_size(trim(obj_file_path), nb_vertex, nb_faces)
        call read_object(trim(obj_file_path), vertex, faces, nb_vertex, nb_faces)

        ! Object is positioned according the user requirements
        vertex(:,1)=vertex(:,1)*body_scale_x1+body_x1*L1
        vertex(:,2)=vertex(:,2)*body_scale_x2+body_x2*2.d0
        vertex(:,3)=vertex(:,3)*body_scale_x3+body_x3*L3

        ! PERFORMING THE MASK
        call perform_mask(vertex, faces, nb_vertex, nb_faces, IBM_mask1, n1, &
        n2,xstart(2),xend(2), n3,xstart(3),xend(3), X1, X2c, X3c, 1.d0)

        call perform_mask(vertex, faces, nb_vertex, nb_faces, IBM_mask2, n1, &
        n2,xstart(2),xend(2), n3,xstart(3),xend(3), X1c, X2, X3c, 1.d0)

        call perform_mask(vertex, faces, nb_vertex, nb_faces, IBM_mask3, n1, &
        n2,xstart(2),xend(2), n3,xstart(3),xend(3), X1c, X2c, X3, 1.d0)

        ! FOR DEBUGGING, export the mask array in snapshot directory
        call create_stretch_snapshot(COMMON_snapshot_path, "IBM", IBM_mask1, "mask1", 1, X1,X2c,X3c)
        call create_stretch_snapshot(COMMON_snapshot_path, "IBM", IBM_mask2, "mask2", 1, X1c,X2,X3c)
        call create_stretch_snapshot(COMMON_snapshot_path, "IBM", IBM_mask3, "mask3", 1, X1c,X2c,X3)

        ! Default: qbound is zero everywhere at the interface body/flow
        qbound=0.d0

    end subroutine IBM_setup

    subroutine force_velocity(q1, q2, q3, sz, vel_term1, vel_term2, vel_term3)

        implicit none
        integer, dimension(3), intent(in)           :: sz
        real*8, dimension(:,:,:), intent(in)        :: q1,q2,q3
        real*8, dimension(:,:,:), intent(out)       :: vel_term1, vel_term2, vel_term3

        select case (interpol_type)

            case (IBM_INTERPOL_NONE)
                call no_interpol

            case default

        end select

    contains

        subroutine no_interpol() !vel_term correspond au terme (V^(n+1) - u^n )/dt

            implicit none
            integer ::   i,j,k

            do k=1, sz(3)
                do j=1, sz(2)
                    do i=1, sz(1)

                        vel_term1(i,j,k) = ( qbound(1) - q1(i,j,k) )  !* solid_cell(j,k)
                        vel_term2(i,j,k) = ( qbound(2) - q2(i,j,k) )  !* solid_cell(j,k)
                        vel_term3(i,j,k) = ( qbound(3) - q3(i,j,k) )  !* solid_cell(j,k)

                    enddo
                enddo
            enddo

        end subroutine no_interpol

    end subroutine force_velocity


end module IBM
