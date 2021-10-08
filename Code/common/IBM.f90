module object_drawer
    implicit none

    public :: perform_mask, perform_modulation_function, get_objet_area
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

    subroutine get_objet_area(vertex, nb_vertex)
        use IBM_data
        implicit none
        integer                                     :: nb_vertex
        real*8, dimension(nb_vertex, 3)             :: vertex

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
        use IBM_data
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
        integer :: i0, i1, j0, j1, k0, k1, ierr
        real*8  :: tolerance=1.d-13

        mask_field=0.d0

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

    subroutine perform_modulation_function(mask_field, n1, n2, n2s,n2e,n3, n3s,n3e, X, Y, Z, L1, L2, L3)
        use IBM_data
        use mathematical_constants
        use mpi
        implicit none
        integer                                     :: n1, n2, n2s,n2e,n3, n3s,n3e
        real*8, dimension(n1, n2s:n2e,n3s:n3e)      :: mask_field
        real*8, dimension(:)                        :: X, Y, Z
        real*8                                      :: L1, L2, L3

        integer     :: i,j,k,f

        integer :: i0, i1, j0, j1, k0, k1, ierr, mpi_err
        real*8  :: x_c, y_c, z_c
        real*8  :: dimension_x, dimension_y, dimension_z
        real*8  :: f_max, f_max_glob

        if ((xs<0).and.(xe>L1)) then
            dimension_x = L1
            x_c = L1/2.d0
        elseif (xs<0) then
            dimension_x = xe
            x_c = xe/2.d0
        elseif (xe>L1) then
            dimension_x = L1 - xs
            x_c = xs + dimension_x/2.d0
        else
            dimension_x = (xe - xs)
            x_c = xs + dimension_x/2.d0
        endif

        if ((ys<0).and.(ye>L2)) then
            dimension_y = L2
            y_c = L2/2.d0
        elseif (ys<0) then
            dimension_y = ye
            y_c = ye/2.d0
        elseif (ye>L2) then
            dimension_y = L2 - ys
            y_c = ys + dimension_y/2.d0
        else
            dimension_y = (ye - ys)
            y_c = ys + dimension_y/2.d0
        endif

        if ((zs<0).and.(ze>L3)) then
            dimension_z = L3
            z_c = L3/2.d0
        elseif (zs<0) then
            dimension_z = ze
            z_c = ze/2.d0
        elseif (ze>L3) then
            dimension_z = L3 - zs
            z_c = zs + dimension_z/2.d0
        else
            dimension_z = (ze - zs)
            z_c = zs + dimension_z/2.d0
        endif

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

        ! At this stage we have the volumic center of the object,
        ! the position of the bounds and the indexes of the domain
        write(*,*)'perform_modulation_function:', x_c, y_c, z_c
        write(*,*)dimension_x, dimension_y, dimension_z
        call sleep(0)

        mask_field = 0d0

        do i=i0,i1
            do j = n2s, n2e
                do k = n3s, n3e

                    if ((j.gt.j0).and.(j.lt.j1).and.(k.gt.k0).and.(k.lt.k1)) then

                        ! mask_field(i,j,k) = ( 0.5 - abs((X(i)-x_c)) / dimension_x) * & 
                        ! ( 0.5 - abs((Y(j)-y_c)) / dimension_y) * &
                        ! ( 0.5 - abs((Z(k)-z_c)) / dimension_z)

                        mask_field(i,j,k) = (1 - abs((X(i)-x_c)/dimension_x + (Y(j)-y_c)/ dimension_y) - abs((Y(j)-y_c)/ dimension_y - (X(i)-x_c)/dimension_x)) * & 
                        (1 - abs((Y(j)-y_c)/dimension_y + (Z(k)-z_c)/ dimension_z) - abs((Z(k)-z_c)/ dimension_z - (Y(j)-y_c)/dimension_y)) * &
                        (1 - abs((Z(k)-z_c)/dimension_z + (X(i)-x_c)/ dimension_x) - abs((X(i)-x_c)/ dimension_x - (Z(k)-z_c)/dimension_z))

                    endif

                enddo
            enddo
        enddo

        f_max = maxval(mask_field)

        call MPI_ALLREDUCE (f_max, f_max_glob, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)

        mask_field = mask_field * (f_max_glob - mask_field)
        ! mask_field = (f_max_glob - mask_field)

        f_max = maxval(mask_field)

        call MPI_ALLREDUCE (f_max, f_max_glob, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)

        mask_field = mask_field / f_max_glob

    end subroutine perform_modulation_function

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
        use DNS_settings
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
        ! Shift in direction 1 with L_y so that it is easier
        ! vertex(:,1)=vertex(:,1)*body_scale_x1+body_x1*2.d0
        vertex(:,2)=vertex(:,2)*body_scale_x2+body_x2*2.d0
        vertex(:,3)=vertex(:,3)*body_scale_x3+body_x3*L3

        call get_objet_area(vertex, nb_vertex)

        ! PERFORMING THE MASK
        call perform_mask(vertex, faces, nb_vertex, nb_faces, IBM_mask1_x, n1, &
        n2,xstart(2),xend(2), n3,xstart(3),xend(3), X1, X2c, X3c, 1.d0)

        call perform_mask(vertex, faces, nb_vertex, nb_faces, IBM_mask2_x, n1, &
        n2,xstart(2),xend(2), n3,xstart(3),xend(3), X1c, X2, X3c, 1.d0)

        call perform_mask(vertex, faces, nb_vertex, nb_faces, IBM_mask3_x, n1, &
        n2,xstart(2),xend(2), n3,xstart(3),xend(3), X1c, X2c, X3, 1.d0)

        ! if (interpol_type == ANTISYMMETRIC_INTERPOL) then
        !     call perform_modulation_function(IBM_modulation_x, n1, &
        !         n2, xstart(2),xend(2), n3,xstart(3),xend(3), X1, X2c, X3c, L1, L2, L3)

        !     call perform_modulation_function(IBM_modulation_y, n1, &
        !         n2, xstart(2),xend(2), n3,xstart(3),xend(3), X1c, X2, X3c, L1, L2, L3)

        !     call perform_modulation_function(IBM_modulation_z, n1, &
        !         n2, xstart(2),xend(2), n3,xstart(3),xend(3), X1c, X2c, X3, L1, L2, L3)
        ! endif

        ! FOR DEBUGGING, export the mask array in snapshot directory
        call create_stretch_snapshot(COMMON_snapshot_path, "IBM", IBM_mask1_x, "mask1", 1, X1,X2c,X3c)
        call create_stretch_snapshot(COMMON_snapshot_path, "IBM", IBM_mask2_x, "mask2", 1, X1c,X2,X3c)
        call create_stretch_snapshot(COMMON_snapshot_path, "IBM", IBM_mask3_x, "mask3", 1, X1c,X2c,X3)

        if (interpol_type == ANTISYMMETRIC_INTERPOL) then
            call create_stretch_snapshot(COMMON_snapshot_path, "IBM", IBM_modulation_x, "modulation1", 1, X1,X2c,X3c)
            call create_stretch_snapshot(COMMON_snapshot_path, "IBM", IBM_modulation_y, "modulation2", 1, X1c,X2,X3c)
            call create_stretch_snapshot(COMMON_snapshot_path, "IBM", IBM_modulation_z, "modulation3", 1, X1c,X2c,X3)
        endif

        ! if (streamwise==3) then
        !     call transpose_mask_x_to_z(IBM_mask1_x, IBM_mask2_x, IBM_mask3_x)
        ! endif

        ! Default: qbound is zero everywhere at the interface body/flow
        qbound=0.d0

    contains

        ! subroutine transpose_mask_x_to_z(IBM_mask1x, IBM_mask2x, IBM_mask3x)

        !     implicit none

        !     real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: IBM_mask1x
        !     real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: IBM_mask2x
        !     real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: IBM_mask3x

        !     real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))                :: IBM_mask_tmp_y


        !     call transpose_x_to_y(IBM_mask1x, IBM_mask_tmp_y)
        !     call transpose_y_to_z(IBM_mask_tmp_y, IBM_mask1_z)
        !     ! q2
        !     call transpose_x_to_y(IBM_mask2x, IBM_mask_tmp_y)
        !     call transpose_y_to_z(IBM_mask_tmp_y, IBM_mask2_z)
        !     ! q3
        !     call transpose_x_to_y(IBM_mask3x, IBM_mask_tmp_y)
        !     call transpose_y_to_z(IBM_mask_tmp_y, IBM_mask3_z)

        ! end subroutine transpose_mask_x_to_z

    end subroutine IBM_setup

    subroutine force_velocity(q1, q2, q3, sz, vel_term1, vel_term2, vel_term3)

        implicit none
        integer, dimension(3), intent(in)           :: sz
        real*8, dimension(:,:,:), intent(in)        :: q1,q2,q3
        real*8, dimension(:,:,:), intent(out)       :: vel_term1, vel_term2, vel_term3

        select case (interpol_type)

            case (IBM_INTERPOL_NONE, SECOND_ORDER_INTERPOL)
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

    subroutine set_antisymmetric_velocity_in_object(X,Y,Z,q,configuration)
        use mathematical_constants
        use mesh, only: n1, n2, n3, L1, L2, L3
        use mpi

        implicit none
        integer                                                         :: i,j,k
        integer                                                         :: it,jt,kt
        real*8                                                          :: x_c, y_c, z_c
        real*8                                                          :: dimension_x, dimension_y, dimension_z
        real*8, dimension(:)                                            :: X, Y, Z
        real*8, dimension(:,:,:)                                        :: q
        character*1                                                     :: configuration
        integer                                                         :: i0, i1, j0, j1, k0, k1
        integer                                                         :: corresponding_i, corresponding_j, corresponding_k

        i0=1; i1=n1; j0=1; j1=n2; k0=1; k1=n3;

        do i = 1, size(X)
            if (X(i)>=xs) then
                i0=i
                exit
            endif
        end do

        do i = size(X), 1, -1
            if (X(i)<=xe) then
                i1=i
                exit
            endif
        end do

        do j = 1, size(Y)
            if (Y(j)>=ys) then
                j0=j
                exit
            endif
        end do

        do j = size(Y), 1, -1
            if (Y(j)<=ye) then
                j1=j
                exit
            endif
        end do

        do k = 1, size(Z)
            if (Z(k)>=zs) then
                k0=k
                exit
            endif
        end do

        do k = size(Z), 1, -1
            if (Z(k)<=ze) then
                k1=k
                exit
            endif
        end do

        if ((xs<0).and.(xe>L1)) then
            dimension_x = L1
            x_c = L1/2.d0
        elseif (xs<0) then
            dimension_x = xe
            x_c = xe/2.d0
        elseif (xe>L1) then
            dimension_x = L1 - xs
            x_c = xs + dimension_x/2.d0
        else
            dimension_x = (xe - xs)
            x_c = xs + dimension_x/2.d0
        endif

        if ((ys<0).and.(ye>L2)) then
            dimension_y = L2
            y_c = L2/2.d0
        elseif (ys<0) then
            dimension_y = ye
            y_c = ye/2.d0
        elseif (ye>L2) then
            dimension_y = L2 - ys
            y_c = ys + dimension_y/2.d0
        else
            dimension_y = (ye - ys)
            y_c = ys + dimension_y/2.d0
        endif

        if ((zs<0).and.(ze>L3)) then
            dimension_z = L3
            z_c = L3/2.d0
        elseif (zs<0) then
            dimension_z = ze
            z_c = ze/2.d0
        elseif (ze>L3) then
            dimension_z = L3 - zs
            z_c = zs + dimension_z/2.d0
        else
            dimension_z = (ze - zs)
            z_c = zs + dimension_z/2.d0
        endif

        if (configuration=='X') then

            do i= 1, size(X)
                do j = 1, xsize(2)
                    do k = 1, xsize(3)

                        jt = min(j+xstart(2)-1,n2-1)
                        kt = min(k+xstart(3)-1,n3-1)

                        if ((X(i).gt.xs).and.(X(i).lt.xe).and.(Y(jt).gt.ys).and.(Y(jt).lt.ye).and.(Z(kt).gt.zs).and.(Z(kt).lt.ze)) then

                            corresponding_i = 0

                            if ((X(i)-xs).eq.0) then 
                                corresponding_i = 2*i0 - i
                            else if ((X(i)-xe).eq.0) then
                                corresponding_i = 2*i1 - i
                            else if ((X(i)-x_c).lt.0) then
                                corresponding_i = 2*i0 - i - 1
                            else 
                                corresponding_i = 2*i1 - i + 1
                            endif

                            if ((corresponding_i.ge.1).and.(corresponding_i.le.(n1-1))) then
                                q(i,j,k) = - sin(2*pi*(X(i)-x_c)**2.d0 / dimension_x**2.d0) * q(corresponding_i,j,k)
                            endif

                        endif
                    enddo
                enddo
            enddo

        endif


        if (configuration=='Y') then

            do i= 1, ysize(1)
                do j = 1, size(Y)
                    do k = 1, ysize(3)

                        it = min(i+ystart(1)-1,n1-1)
                        kt = min(k+ystart(3)-1,n3-1)

                        if ((X(it).gt.xs).and.(X(it).lt.xe).and.(Y(j).gt.ys).and.(Y(j).lt.ye).and.(Z(kt).gt.zs).and.(Z(kt).lt.ze)) then

                            corresponding_j = 0

                            if ((Y(j)-ys).eq.0) then 
                                corresponding_j = 2*j0 - j
                            else if ((Y(j)-ye).eq.0) then
                                corresponding_j = 2*j1 - j
                            else if ((Y(j)-y_c).lt.0) then
                                corresponding_j = 2*j0 - j - 1
                            else 
                                corresponding_j = 2*j1 - j + 1
                            endif

                            if ((corresponding_j.ge.1).and.(corresponding_j.le.(n2-1))) then
                                q(i,j,k) = - sin(2*pi*(Y(j)-y_c)**2.d0 / dimension_y**2.d0) * q(i,corresponding_j,k)
                            endif

                        endif
                    enddo
                enddo
            enddo

        endif


        if (configuration=='Z') then

            do i= 1, zsize(1)
                do j = 1, zsize(2)
                    do k = 1, size(Z)

                        it = min(i+zstart(1)-1,n1-1)
                        jt = min(j+zstart(2)-1,n2-1)

                        if ((X(it).gt.xs).and.(X(it).lt.xe).and.(Y(jt).gt.ys).and.(Y(jt).lt.ye).and.(Z(k).gt.zs).and.(Z(k).lt.ze)) then

                            corresponding_k = 0

                            if ((Z(k)-zs).eq.0) then 
                                corresponding_k = 2*k0 - k
                            else if ((Z(k)-ze).eq.0) then
                                corresponding_k = 2*k1 - k
                            else if ((Z(k)-z_c).lt.0) then
                                corresponding_k = 2*k0 - k - 1
                            else 
                                corresponding_k = 2*k1 - k + 1
                            endif

                            if ((corresponding_k.ge.1).and.(corresponding_k.le.(n3-1))) then
                                q(i,j,k) = - sin(2*pi*(Z(k)-z_c)**2.d0 / dimension_z**2.d0) * q(i,j,corresponding_k)
                            endif

                        endif
                    enddo
                enddo
            enddo

        endif

    end subroutine set_antisymmetric_velocity_in_object

    subroutine get_antisymmetric_velocity_in_object(X,Y,Z,q,q_out,configuration)
        use mathematical_constants
        use mesh, only: n1, n2, n3, L1, L2, L3
        use mpi

        implicit none
        integer                                                         :: i,j,k
        integer                                                         :: it,jt,kt
        real*8                                                          :: x_c, y_c, z_c
        real*8                                                          :: dimension_x, dimension_y, dimension_z
        real*8, dimension(:)                                            :: X, Y, Z
        real*8, dimension(:,:,:)                                        :: q, q_out
        character*1                                                     :: configuration
        integer                                                         :: i0, i1, j0, j1, k0, k1
        integer                                                         :: corresponding_i, corresponding_j, corresponding_k

        i0=1; i1=n1; j0=1; j1=n2; k0=1; k1=n3;

        do i = 1, size(X)
            if (X(i)>=xs) then
                i0=i
                exit
            endif
        end do

        do i = size(X), 1, -1
            if (X(i)<=xe) then
                i1=i
                exit
            endif
        end do

        do j = 1, size(Y)
            if (Y(j)>=ys) then
                j0=j
                exit
            endif
        end do

        do j = size(Y), 1, -1
            if (Y(j)<=ye) then
                j1=j
                exit
            endif
        end do

        do k = 1, size(Z)
            if (Z(k)>=zs) then
                k0=k
                exit
            endif
        end do

        do k = size(Z), 1, -1
            if (Z(k)<=ze) then
                k1=k
                exit
            endif
        end do

        if ((xs<0).and.(xe>L1)) then
            dimension_x = L1
            x_c = L1/2.d0
        elseif (xs<0) then
            dimension_x = xe
            x_c = xe/2.d0
        elseif (xe>L1) then
            dimension_x = L1 - xs
            x_c = xs + dimension_x/2.d0
        else
            dimension_x = (xe - xs)
            x_c = xs + dimension_x/2.d0
        endif

        if ((ys<0).and.(ye>L2)) then
            dimension_y = L2
            y_c = L2/2.d0
        elseif (ys<0) then
            dimension_y = ye
            y_c = ye/2.d0
        elseif (ye>L2) then
            dimension_y = L2 - ys
            y_c = ys + dimension_y/2.d0
        else
            dimension_y = (ye - ys)
            y_c = ys + dimension_y/2.d0
        endif

        if ((zs<0).and.(ze>L3)) then
            dimension_z = L3
            z_c = L3/2.d0
        elseif (zs<0) then
            dimension_z = ze
            z_c = ze/2.d0
        elseif (ze>L3) then
            dimension_z = L3 - zs
            z_c = zs + dimension_z/2.d0
        else
            dimension_z = (ze - zs)
            z_c = zs + dimension_z/2.d0
        endif

        if (configuration=='X') then

            do i= 1, size(X)
                do j = 1, xsize(2)
                    do k = 1, xsize(3)

                        jt = min(j+xstart(2)-1,n2-1)
                        kt = min(k+xstart(3)-1,n3-1)

                        if ((X(i).gt.xs).and.(X(i).lt.xe).and.(Y(jt).gt.ys).and.(Y(jt).lt.ye).and.(Z(kt).gt.zs).and.(Z(kt).lt.ze)) then

                            corresponding_i = 0

                            if ((X(i)-xs).eq.0) then 
                                corresponding_i = 2*i0 - i
                            else if ((X(i)-xe).eq.0) then
                                corresponding_i = 2*i1 - i
                            else if ((X(i)-x_c).lt.0) then
                                corresponding_i = 2*i0 - i - 1
                            else 
                                corresponding_i = 2*i1 - i + 1
                            endif

                            if ((corresponding_i.ge.1).or.(corresponding_i.le.(n1-1))) then
                                q_out(i,j,k) = - sin(2*pi*(X(i)-x_c)**2.d0 / dimension_x**2.d0) * sin(2*pi*(Y(j)-y_c)**2.d0 / dimension_y**2.d0) * sin(2*pi*(Z(k)-z_c)**2.d0 / dimension_z**2.d0) * q(corresponding_i,j,k)
                                ! q_out(i,j,k) = - IBM_modulation_x(i,jt,kt) * q(corresponding_i,j,k)
                            endif

                        endif
                    enddo
                enddo
            enddo

        endif


        if (configuration=='Y') then

            do i= 1, ysize(1)
                do j = 1, size(Y)
                    do k = 1, ysize(3)

                        it = min(i+ystart(1)-1,n1-1)
                        kt = min(k+ystart(3)-1,n3-1)

                        if ((X(it).gt.xs).and.(X(it).lt.xe).and.(Y(j).gt.ys).and.(Y(j).lt.ye).and.(Z(kt).gt.zs).and.(Z(kt).lt.ze)) then

                            corresponding_j = 0

                            if ((Y(j)-ys).eq.0) then 
                                corresponding_j = 2*j0 - j
                            else if ((Y(j)-ye).eq.0) then
                                corresponding_j = 2*j1 - j
                            else if ((Y(j)-y_c).lt.0) then
                                corresponding_j = 2*j0 - j - 1
                            else 
                                corresponding_j = 2*j1 - j + 1
                            endif

                            if ((corresponding_j.ge.1).or.(corresponding_j.le.(n2-1))) then
                                q_out(i,j,k) = - sin(2*pi*(Y(j)-y_c)**2.d0 / dimension_y**2.d0) * q(i,corresponding_j,k)
                            endif

                        endif
                    enddo
                enddo
            enddo

        endif


        if (configuration=='Z') then

            do i= 1, zsize(1)
                do j = 1, zsize(2)
                    do k = 1, size(Z)

                        it = min(i+zstart(1)-1,n1-1)
                        jt = min(j+zstart(2)-1,n2-1)

                        if ((X(it).gt.xs).and.(X(it).lt.xe).and.(Y(jt).gt.ys).and.(Y(jt).lt.ye).and.(Z(k).gt.zs).and.(Z(k).lt.ze)) then

                            corresponding_k = 0

                            if ((Z(k)-zs).eq.0) then 
                                corresponding_k = 2*k0 - k
                            else if ((Z(k)-ze).eq.0) then
                                corresponding_k = 2*k1 - k
                            else if ((Z(k)-z_c).lt.0) then
                                corresponding_k = 2*k0 - k - 1
                            else 
                                corresponding_k = 2*k1 - k + 1
                            endif

                            if ((corresponding_k.ge.1).or.(corresponding_k.le.(n3-1))) then
                                q_out(i,j,k) = - sin(2*pi*(Z(k)-z_c)**2.d0 / dimension_z**2.d0) * q(i,j,corresponding_k)
                            endif

                        endif
                    enddo
                enddo
            enddo

        endif

    end subroutine get_antisymmetric_velocity_in_object

    subroutine reset_velocity(X,Y,Z,q,configuration)
        use mathematical_constants
        use mesh, only: n1, n2, n3, L1, L2, L3

        implicit none
        integer                                                         :: i,j,k
        real*8                                                          :: x_c, y_c, z_c
        real*8                                                          :: dimension_x, dimension_y, dimension_z
        real*8, dimension(:)                                            :: X, Y, Z
        real*8, dimension(:,:,:)                                        :: q
        character*1                                                     :: configuration
        integer                                                         :: i0, i1, j0, j1, k0, k1
        integer                                                         :: corresponding_i, corresponding_j, corresponding_k

        i0=1; i1=n1; j0=1; j1=n2; k0=1; k1=n3;

        do i = 1, size(X)
            if (X(i)>=xs) then
                i0=i
                exit
            endif
        end do

        do i = size(X), 1, -1
            if (X(i)<=xe) then
                i1=i
                exit
            endif
        end do

        do j = 1, size(Y)
            if (Y(j)>=ys) then
                j0=j
                exit
            endif
        end do

        do j = size(Y), 1, -1
            if (Y(j)<=ye) then
                j1=j
                exit
            endif
        end do

        do k = 1, size(Z)
            if (Z(k)>=zs) then
                k0=k
                exit
            endif
        end do

        do k = size(Z), 1, -1
            if (Z(k)<=ze) then
                k1=k
                exit
            endif
        end do

        if (configuration=='X') then

            do i=i0,i1
                ! do j = xstart(2), min(xend(2),n2-1)
                !     do k = xstart(3), min(xend(3),n3-1)
                do j = 1, xsize(2)
                    do k = 1, xsize(3)

                        if (((j+xstart(2)-1).ge.j0).and.((j+xstart(2)-1).le.j1).and.((k+xstart(3)-1).ge.k0).and.((k+xstart(3)-1).le.k1)) then

                            q(i,j,k) = 0.d0

                        endif
                    enddo
                enddo
            enddo

        endif

        if (configuration=='Y') then

            ! do i=ystart(1), min(yend(1),n1-1)
            do i= 1, ysize(1)
                do j = j0,j1
                    do k = 1, ysize(3)
                    ! do k = ystart(3), min(yend(3),n3-1)

                        if (((i+ystart(1)-1).ge.i0).and.((i+ystart(1)-1).le.i1).and.((k+ystart(3)-1).ge.k0).and.((k+ystart(3)-1).le.k1)) then

                            q(i,j,k) = 0.d0

                        endif
                    enddo
                enddo
            enddo

        endif

        if (configuration=='Z') then

            ! do i=zstart(1), min(zend(1),n1-1)
            !     do j = zstart(2), min(zend(2),n2-1)
            do i= 1, zsize(1)
                do j = 1, zsize(2)
                    do k = k0,k1

                        if (((i+zstart(1)-1).ge.i0).and.((i+zstart(1)-1).le.i1).and.((j+zstart(2)-1).ge.j0).and.((j+zstart(2)-1).le.j1)) then

                            q(i,j,k) = 0.d0

                        endif
                    enddo
                enddo
            enddo

        endif

    end subroutine reset_velocity


    subroutine compute_antisymmetric_velocity(q1, q2, q3, sz, vel_term1, vel_term2, vel_term3)

        use mesh
        use snapshot_writer
        use COMMON_workspace_view

        implicit none
        integer, dimension(3), intent(in)           :: sz
        real*8, dimension(:,:,:), intent(in)        :: q1,q2,q3
        real*8, dimension(:,:,:), intent(out)       :: vel_term1, vel_term2, vel_term3
        real*8  :: X1(n1), X2(n2), X3(n3)
        real*8  :: X1c(n1), X2c(n2), X3c(n3)
        integer     :: i,j,k

        ! Z => X
        ! Y => Y
        ! X => Z
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

        vel_term1 = 0.d0
        vel_term2 = 0.d0
        vel_term3 = 0.d0

        call antisymmetric_interpol(vel_term1, IBM_modulation_x, q1, n1, n2, n3, X1, X2c, X3c, .true., .true., .true.)

        call antisymmetric_interpol(vel_term2, IBM_modulation_y, q2, n1, n2, n3, X1c, X2, X3c, .true., .true., .true.)

        call antisymmetric_interpol(vel_term3, IBM_modulation_z, q3, n1, n2, n3, X1c, X2c, X3, .true., .true., .true.)


        call create_stretch_snapshot(COMMON_snapshot_path, "IBM", vel_term1, "vel_term1", 1, X1,X2c,X3c)
        call create_stretch_snapshot(COMMON_snapshot_path, "IBM", vel_term2, "vel_term2", 1, X1c,X2,X3c)
        call create_stretch_snapshot(COMMON_snapshot_path, "IBM", vel_term3, "vel_term3", 1, X1c,X2c,X3)


        do k=1, sz(3)
            do j=1, sz(2)
                do i=1, sz(1)

                    vel_term1(i,j,k) = ( vel_term1(i,j,k) - q1(i,j,k) )  !* solid_cell(j,k)
                    vel_term2(i,j,k) = ( vel_term2(i,j,k) - q2(i,j,k) )  !* solid_cell(j,k)
                    vel_term3(i,j,k) = ( vel_term3(i,j,k) - q3(i,j,k) )  !* solid_cell(j,k)

                enddo
            enddo
        enddo

        contains

            subroutine antisymmetric_interpol(vel_term, IBM_modulation, q, n1, n2, n3, X, Y, Z, ax, ay, az) !vel_term correspond au terme (V^(n+1) - u^n )/dt
                use mpi

                implicit none
                integer                                                         :: i,j,k
                integer                                                         :: n1, n2, n3
                real*8, dimension(n1, xstart(2):xend(2), xstart(3):xend(3))     :: IBM_modulation, q, vel_term
                real*8, dimension(ystart(1):yend(1), n2, ystart(3):yend(3))     :: q_y
                real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), n3)     :: q_z
                ! real*8, dimension(:,:,:), allocatable                           :: q_buffer
                real*8                                                          :: x_c, y_c, z_c
                real*8                                                          :: dimension_x, dimension_y, dimension_z
                real*8, dimension(:)                                            :: X, Y, Z
                integer                                                         :: i0, i1, j0, j1, k0, k1, mpi_err
                real*8, dimension(:,:,:), allocatable                           :: velocity_x, velocity_y, velocity_z, velocity_y_x, velocity_z_x
                real*8, dimension(:,:,:), allocatable                           :: corresponding_i_x, corresponding_i_y, corresponding_i_z
                real*8, dimension(:,:,:), allocatable                           :: corresponding_j_x, corresponding_j_y, corresponding_j_z
                real*8, dimension(:,:,:), allocatable                           :: corresponding_k_x, corresponding_k_y, corresponding_k_z
                real*8, dimension(:,:,:), allocatable                           :: denom_x, denom_y, denom_z
                logical                                                         :: extrapolate
                logical                                                         :: ax, ay, az
                ! integer                                                         :: size_q_buffer, size_q, errclass, resultlen
                ! integer, dimension(3)                                           :: displs, recv_buffer
                ! character(200)                  :: err_buffer

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

                if ((xs<0).and.(xe>L1)) then
                    dimension_x = L1
                    x_c = L1/2.d0
                elseif (xs<0) then
                    dimension_x = xe
                    x_c = xe/2.d0
                elseif (xe>L1) then
                    dimension_x = L1 - xs
                    x_c = xs + dimension_x/2.d0
                else
                    dimension_x = (xe - xs)
                    x_c = xs + dimension_x/2.d0
                endif

                if ((ys<0).and.(ye>L2)) then
                    dimension_y = L2
                    y_c = L2/2.d0
                elseif (ys<0) then
                    dimension_y = ye
                    y_c = ye/2.d0
                elseif (ye>L2) then
                    dimension_y = L2 - ys
                    y_c = ys + dimension_y/2.d0
                else
                    dimension_y = (ye - ys)
                    y_c = ys + dimension_y/2.d0
                endif

                if ((zs<0).and.(ze>L3)) then
                    dimension_z = L3
                    z_c = L3/2.d0
                elseif (zs<0) then
                    dimension_z = ze
                    z_c = ze/2.d0
                elseif (ze>L3) then
                    dimension_z = L3 - zs
                    z_c = zs + dimension_z/2.d0
                else
                    dimension_z = (ze - zs)
                    z_c = zs + dimension_z/2.d0
                endif

                allocate(corresponding_i_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
                allocate(corresponding_i_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
                allocate(corresponding_i_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

                allocate(corresponding_j_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
                allocate(corresponding_j_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
                allocate(corresponding_j_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

                allocate(corresponding_k_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
                allocate(corresponding_k_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
                allocate(corresponding_k_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))



                ! allocate(corresponding_j(i0:i1, max(j0,xstart(2)):min(j1,xend(2)), max(k0,xstart(3)):min(k1,xend(3))))
                ! allocate(corresponding_k(i0:i1, max(j0,xstart(2)):min(j1,xend(2)), max(k0,xstart(3)):min(k1,xend(3))))

                ! allocate(velocity_x(i0:i1, max(j0,xstart(2)):min(j1,xend(2)), max(k0,xstart(3)):min(k1,xend(3))))
                ! allocate(velocity_y(i0:i1, max(j0,xstart(2)):min(j1,xend(2)), max(k0,xstart(3)):min(k1,xend(3))))
                ! allocate(velocity_z(i0:i1, max(j0,xstart(2)):min(j1,xend(2)), max(k0,xstart(3)):min(k1,xend(3))))

                allocate(velocity_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
                allocate(velocity_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
                allocate(velocity_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

                allocate(velocity_y_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
                allocate(velocity_z_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))

                allocate(denom_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
                allocate(denom_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
                allocate(denom_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

                corresponding_i_x=0.d0
                corresponding_i_y=0.d0
                corresponding_i_z=0.d0

                corresponding_j_x=0.d0
                corresponding_j_y=0.d0
                corresponding_j_z=0.d0

                corresponding_k_x=0.d0
                corresponding_k_y=0.d0
                corresponding_k_z=0.d0

                velocity_x=0.d0
                velocity_y=0.d0
                velocity_z=0.d0
                
                denom_x=0.d0
                denom_y=0.d0
                denom_z=0.d0


                ! allocate(q_buffer(n1,n2,n3))
                ! q_buffer = 0.d0
                ! allocate(q_buffer(n1,n2,n3))
                ! q_buffer = 0.d0

                ! size_q = size(q)
                ! size_q_buffer = size(q_buffer)
                ! displs(1) = 0
                ! displs(2) = xstart(2) - 1
                ! displs(3) = xstart(3) - 1

                ! recv_buffer(1) = n1
                ! recv_buffer(2) = xend(2) - xstart(2) + 1
                ! recv_buffer(3) = xend(3) - xstart(3) + 1



                ! write(*,*) 'nrank=', nrank, 'displs', displs
                ! write(*,*) 'nrank=', nrank, 'recv_buffer', recv_buffer
                ! write(*,*) 'nrank=', nrank, 'xsize', xsize
                ! write(*,*) 'nrank=', nrank, recv_buffer(1)*recv_buffer(2)*recv_buffer(3), size_q

                ! ! call MPI_ALLGATHER(q(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), size_q, MPI_REAL8, & 
                ! !     q_buffer(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), size_q_buffer, MPI_REAL8, MPI_COMM_WORLD, mpi_err)

                ! call MPI_ALLGATHERV(q(:,:,:), size_q, MPI_REAL8, & 
                !     q_buffer(:,:,:), recv_buffer, displs, MPI_REAL8, MPI_COMM_WORLD, mpi_err)

                ! ! call MPI_ALLGATHER(q(xstart(1):min(xend(1),n1m),xstart(2):min(xend(2),n2m),xstart(3):min(xend(3),n3m), size_q, MPI_REAL8, & 
                ! !     q_buffer(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), size_q_buffer, MPI_REAL8, MPI_COMM_WORLD, mpi_err)

                ! write(*,*) 'after aie aie aie'

                ! write(*,*) sum(abs(q_buffer))

                ! Get the velocities from the x direction

                if (ax) then

                    do i=i0,i1
                        do j = xstart(2), xend(2)
                            do k = xstart(3), xend(3)

                                if ((j.gt.j0).and.(j.lt.j1).and.(k.gt.k0).and.(k.lt.k1)) then

                                    if ((X(i)-x_c).gt.0) then 
                                        corresponding_i_x(i,j,k) = i1
                                    else
                                        corresponding_i_x(i,j,k) = i0
                                    endif

                                    if (((2*corresponding_i_x(i,j,k)-i).le.0).or.((2*corresponding_i_x(i,j,k)-i).ge.(n1-1))) then
                                        velocity_x(i,j,k) = 0.d0
                                    else
                                        velocity_x(i,j,k) = q(2*int(corresponding_i_x(i,j,k))-i,j,k)
                                        denom_x(i,j,k) = denom_x(i,j,k) + real(abs(corresponding_i_x(i,j,k)-i),8)
                                    endif

                                endif
                            enddo
                        enddo
                    enddo

                endif

                ! Get the velocities from the y direction
                call transpose_x_to_y(q, q_y)
                call transpose_x_to_y(denom_x, denom_y)

                if (ay) then

                    do i=ystart(1), yend(1)
                        do j = j0,j1
                            do k = ystart(3), yend(3)

                                if ((i.gt.i0).and.(i.lt.i1).and.(k.gt.k0).and.(k.lt.k1)) then

                                    if ((Y(j)-y_c).gt.0) then 
                                        corresponding_j_y(i,j,k) = j1
                                    else
                                        corresponding_j_y(i,j,k) = j0
                                    endif

                                    if (((2*corresponding_j_y(i,j,k)-j).le.0).or.((2*corresponding_j_y(i,j,k)-j).ge.(n2-1))) then
                                        velocity_y(i,j,k) = 0.d0
                                    else
                                        velocity_y(i,j,k) = q_y(i, 2*int(corresponding_j_y(i,j,k))-j,k)
                                        denom_y(i,j,k) = denom_y(i,j,k) + real(abs(corresponding_j_y(i,j,k)-j),8)
                                    endif

                                endif
                            enddo
                        enddo
                    enddo

                endif

                ! Get the velocities from the z direction
                call transpose_y_to_z(q_y, q_z)
                call transpose_y_to_z(denom_y, denom_z)

                if (az) then

                    do i=zstart(1), zend(1)
                        do j = zstart(2), zend(2)
                            do k = k0,k1

                                if ((i.gt.i0).and.(i.lt.i1).and.(j.gt.j0).and.(j.lt.j1)) then

                                    if ((Z(k)-z_c).gt.0) then 
                                        corresponding_k_z(i,j,k) = k1
                                    else
                                        corresponding_k_z(i,j,k) = k0
                                    endif

                                    if (((2*corresponding_k_z(i,j,k)-k).le.0).or.((2*corresponding_k_z(i,j,k)-k).ge.(n3-1))) then
                                        velocity_z(i,j,k) = 0.d0
                                    else
                                        velocity_z(i,j,k) = q_z(i,j, 2*int(corresponding_k_z(i,j,k))-k)
                                        denom_z(i,j,k) = denom_z(i,j,k) + real(abs(corresponding_k_z(i,j,k)-k),8)
                                    endif

                                endif
                            enddo
                        enddo
                    enddo

                endif

                call transpose_z_to_y(denom_z, denom_y)
                call transpose_y_to_x(denom_y, denom_x)

                call transpose_z_to_y(corresponding_k_z, corresponding_k_y)
                call transpose_y_to_x(corresponding_k_y, corresponding_k_x)

                call transpose_y_to_x(corresponding_j_y, corresponding_j_x)

                ! Extrapolate the velocity inside the object
                do i=i0,i1
                    do j = xstart(2), xend(2)
                        do k = xstart(3), xend(3)

                            if ((j.gt.j0).and.(j.lt.j1).and.(k.gt.k0).and.(k.lt.k1)) then

                                if (denom_x(i,j,k).eq.0) then
                                    vel_term(i,j,k) = 0.d0
                                else

                                    vel_term(i,j,k) = - IBM_modulation(i,j,k) * ( &
                                        (denom_x(i,j,k) - real(abs(corresponding_i_x(i,j,k)-i),8))/(2*denom_x(i,j,k)) * velocity_x(i,j,k) + &
                                        (denom_x(i,j,k) - real(abs(corresponding_j_x(i,j,k)-j),8))/(2*denom_x(i,j,k)) * velocity_y_x(i,j,k) + &
                                        (denom_x(i,j,k) - real(abs(corresponding_k_x(i,j,k)-k),8))/(2*denom_x(i,j,k)) * velocity_z_x(i,j,k) )
                                endif

                            endif
                        enddo
                    enddo
                enddo

                ! call create_stretch_snapshot(COMMON_snapshot_path, "IBM", denom_x, "denom_x", 1, X,Y,Z)

                ! call create_stretch_snapshot(COMMON_snapshot_path, "IBM", velocity_x, "velocity_x", 1, X1,X2c,X3c)
                ! call create_stretch_snapshot(COMMON_snapshot_path, "IBM", velocity_y_x, "velocity_y_x", 1, X1,X2c,X3c)
                ! call create_stretch_snapshot(COMMON_snapshot_path, "IBM", velocity_z_x, "velocity_z_x", 1, X1,X2c,X3c)

                ! call create_stretch_snapshot(COMMON_snapshot_path, "IBM", corresponding_i_x, "corresponding_i_x", 1, X1,X2c,X3c)
                ! call create_stretch_snapshot(COMMON_snapshot_path, "IBM", corresponding_j_x, "corresponding_j_x", 1, X1,X2c,X3c)
                ! call create_stretch_snapshot(COMMON_snapshot_path, "IBM", corresponding_k_x, "corresponding_k_x", 1, X1,X2c,X3c)

                ! deallocate everything
                deallocate(corresponding_i_x)
                deallocate(corresponding_i_y)
                deallocate(corresponding_i_z)

                deallocate(corresponding_j_x)
                deallocate(corresponding_j_y)
                deallocate(corresponding_j_z)

                deallocate(corresponding_k_x)
                deallocate(corresponding_k_y)
                deallocate(corresponding_k_z)

                deallocate(velocity_x)
                deallocate(velocity_y)
                deallocate(velocity_z)

                deallocate(denom_x)
                deallocate(denom_y)
                deallocate(denom_z)

            end subroutine antisymmetric_interpol

    end subroutine compute_antisymmetric_velocity

    ! return bounds of the object, +/- 6 nodes
    subroutine get_ijk_IBM(X,Y,Z,i_IBM_start,i_IBM_end,j_IBM_start,j_IBM_end,k_IBM_start,k_IBM_end)

        use mesh, only:n1,n2,n3
        use snapshot_writer
        use COMMON_workspace_view

        implicit none
        real*8, dimension(:), intent(in)    :: X,Y,Z
        integer                             :: i,j,k
        integer, intent(inout)              :: i_IBM_start,i_IBM_end,j_IBM_start,j_IBM_end,k_IBM_start,k_IBM_end

        i_IBM_start=1; i_IBM_end=n1; j_IBM_start=1; j_IBM_end=n2; k_IBM_start=1; k_IBM_end=n3;

        do i = 1, n1-1
            if (X(i)>=xs) then
                i_IBM_start=i
                exit
            endif
        end do

        do i = n1-1, 1, -1
            if (X(i)<=xe) then
                i_IBM_end=i
                exit
            endif
        end do

        do j = 1, n2-1
            if (Y(j)>=ys) then
                j_IBM_start=j
                exit
            endif
        end do

        do j = n2-1, 1, -1
            if (Y(j)<=ye) then
                j_IBM_end=j
                exit
            endif
        end do

        do k = 1, n3-1
            if (Z(k)>=zs) then
                k_IBM_start=k
                exit
            endif
        end do

        do k = n3-1, 1, -1
            if (Z(k)<=ze) then
                k_IBM_end=k
                exit
            endif
        end do

        i_IBM_start = i_IBM_start - 6
        i_IBM_end = i_IBM_end + 6

        j_IBM_start = j_IBM_start - 6
        j_IBM_end = j_IBM_end + 6

        k_IBM_start = k_IBM_start - 6
        k_IBM_end = k_IBM_end + 6

    end subroutine get_ijk_IBM


    subroutine activate_O2_for_pressure_correction()

        use schemes_loader
        implicit none

        call configure_schemes(O2_SCHEMES)

    end subroutine activate_O2_for_pressure_correction

    subroutine deactivate_O2_for_intermediate_velocity()

        use schemes_loader
        use numerical_methods_settings
        implicit none

        call configure_schemes(schemes_configuration)

    end subroutine deactivate_O2_for_intermediate_velocity

end module IBM
