module object_drawer
    implicit none

    public :: perform_mask, perform_modulation_function, generate_mesh_ibm, get_objet_area, extrapolate_mask_fine_mesh_ibm
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
        use mesh, only:L1,L2,L3
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

        ! mask_field=0.d0
        ! initialized outside of the subroutine to take into account multi object

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

        ibm_volume = ibm_volume + (min(xe,L1) - max(xs,0.d0)) * (min(ye,L2) - max(ys,0.d0)) * (min(ze,L3) - max(zs,0.d0))

        do i=i0,i1
            do j = n2s, n2e
                do k = n3s, n3e

                    if ((j.ge.j0).and.(j.le.j1).and.(k.ge.k0).and.(k.le.k1)) mask_field(i,j,k)=1.d0

                enddo
            enddo
        enddo

        ! do j = n2s, n2e
        !     do k = n3s, n3e

        !         nb_xintersec=0
        !         x_intersec=0.d0

        !         do f = 1, nb_faces
        !             face1=vertex(faces(f,1), :)
        !             face2=vertex(faces(f,2), :)
        !             face3=vertex(faces(f,3), :)
        !             call intersection(face1, face2, face3, Y(j),Z(k), x_i, ierr)

        !             if (ierr==0) then

        !                 point=(/x_i,Y(j),Z(k)/)

        !                 call inTriangle(face1, face2, face3, point, isInFace, tolerance)

        !                 if (isInFace) then

        !                     call inTriangle(face1, face2, face3, point, isInFace, tolerance)
        !                     nb_xintersec=nb_xintersec+1
        !                     x_intersec(nb_xintersec)=x_i
        !                     face_intersec(nb_xintersec)=f

        !                     if (present(verbose)) then
        !                         write(*,*)
        !                         write(*,*)'Face', f
        !                         write(*,*)'sommet1', face1
        !                         write(*,*)'sommet2', face2
        !                         write(*,*)'sommet3', face3

        !                         write(*,*) 'Ligne y=', Y(j), 'z=', Z(k)
        !                         write(*,*)'Point dintersection', x_i
        !                         write(*,*)
        !                     endif

        !                 endif

        !             endif
        !         end do

        !         if(present(verbose)) write(*,*)'j,k', j,k
        !         if(present(verbose)) write(*,*) 'x_intersec', x_intersec(1:nb_xintersec)
        !         call sort_array(x_intersec, nb_xintersec)

        !         if(mod(nb_xintersec,2)==1) then

        !             if(present(verbose)) then
        !                 write(*,*)
        !                 write(*,*)
        !                 write(*,*) 'Nombre dintersection impair ', nb_xintersec
        !                 write(*,*) 'j,k', j,k
        !                 write(*,*)'Y',Y(j)
        !                 write(*,*)'Z',Z(k)

        !                 write(*,*) 'Intersections en:', x_intersec(1:nb_xintersec)

        !                 write(*,*)
        !                 write(*,*) 'Recherche de doublons...'
        !             endif

        !             x_intersec_clean(1)=x_intersec(1)
        !             pc=1
        !             do i = 2, nb_xintersec
        !                 if(abs(x_intersec(i)-x_intersec_clean(pc))>1.d-9) then
        !                     pc=pc+1
        !                     x_intersec_clean(pc)=x_intersec(i)
        !                 endif
        !             end do
        !             nb_xintersec=pc
        !             x_intersec=0.d0
        !             x_intersec(1:nb_xintersec)=x_intersec_clean(1:pc)

        !             if(present(verbose)) write(*,*)
        !             if(present(verbose)) write(*,*) 'Tableaux des intersections nettoye'
        !             if(present(verbose)) write(*,*) 'Intersections en:', x_intersec(1:nb_xintersec)
        !             if(present(verbose)) call sleep(0)
        !         endif

        !         if(mod(nb_xintersec,2)==1) then
        !             write(*,*) 'ERROR!!!', j,k, nb_xintersec
        !             write(*,*) 'x_intersec', x_intersec(1:nb_xintersec)
        !             write(*,*)'face_intersec', face_intersec(1:nb_xintersec)
        !         endif

        !         pc=1

        !         if ((nb_xintersec>0).and.(mod(nb_xintersec,2)==0)) then

        !             if(present(verbose)) write(*,*) 'x_intersec', x_intersec(1:nb_xintersec)

        !             do i = 1, n1
        !                 if ((X(i)>x_intersec(pc)).and.(X(i)<x_intersec(pc+1))) then
        !                     mask_field(i,j,k)=value

        !                     if (i>10000) then   ! POUR DEBUGGER CHANGER LA VALEUR DU SEUIL
        !                         write(*,*) 'i',i
        !                         write(*,*) 'pc',pc, '/', nb_xintersec
        !                         write(*,*) 'x_intersec', x_intersec
        !                         write(*,*) 'j,k', j,k
        !                         call sleep(2)
        !                     endif

        !                 else if((X(i)>x_intersec(pc)).and.(X(i)>x_intersec(pc+1))) then
        !                     if (nb_xintersec>=pc+3) pc=pc+2
        !                 else
        !                     mask_field(i,j,k)=0.d0
        !                 endif
        !             end do

        !             if(present(verbose)) write(*,*)'+++++++++++++++++++++++++++++++++++++++'

        !         endif

        !     end do
        ! end do

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

    ! Define all the mesh attributes from the mesh description given by arguments (n1,n2...)
    !   ni : number of point along the i th direction
    !   Li : Lenght of the computational domain in the i th direction
    !   str: The stretching applied along the Y direction
    ! THIS IS ONLY FOR THE REFINED GRID USED WITH THE IBM
    subroutine generate_mesh_ibm()

        ! Data modules
        use irregular_derivative_coefficients
        use mesh
        use IBM_data

        ! Processing modules
        use lamballais_transfo
        use arctan_transfo

        use workspace_view, only: log_path
        use transfo_tester

        implicit none

        ! Local variables
        integer                     :: i, j, k
        real*8                      :: dy, yeta
        integer, parameter          :: YMESH_UNIT=47

        integer                     :: writer_proc=0
        character(100)              :: mesh_generator_path, mesh_generator_path_y
        integer                     :: cell_center_file_id, cell_Yface_file_id, zone_id, var_id


        procedure(t1_init_transfo), pointer             :: init_transfo
        procedure(t1_get_computational_coords), pointer :: get_computational_coords
        procedure(t1_get_physical_coords), pointer      :: get_physical_coords
        procedure(t1_DYonDy), pointer                   :: DYonDy
        procedure(t1_D2YonDy2), pointer                 :: D2YonDy2

        init_transfo            =>t1_init_transfo
        get_computational_coords=>t1_get_computational_coords
        get_physical_coords     =>t1_get_physical_coords
        DYonDy                  =>t1_DYonDy
        D2YonDy2                =>t1_D2YonDy2

        n1_ibm = 2*(n1-1)+1
        n2_ibm = 2*(n2-1)+1
        n3_ibm = 2*(n3-1)+1

        ! n1_ibm = 2*n1
        ! n2_ibm = 2*n2
        ! n3_ibm = 2*n3

        allocate(Z_ibm(n1_ibm))
        allocate(Zc_ibm(n1_ibm-1))

        allocate(Y_ibm(n2_ibm))
        allocate(Yc_ibm(n2_ibm-1))

        allocate(X_ibm(n3_ibm))
        allocate(Xc_ibm(n3_ibm-1))

        allocate(cell_size_Y_ibm(n2_ibm-1))

        yeta=0.d0

        dx2_ibm=dx2/2    ! half height of a regular cell
        dx1_ibm=dx1/2     ! spanwise size of a cell
        dx3_ibm=dx3/2      ! streamwise size of a cell

        call init_transfo(stretch_Y)

        ! Point positions in X, Z uniform directions -------------------
        do i=1,n1_ibm
            Z_ibm(i)=dfloat(i-1)*dx1_ibm
        enddo
        do i=1,n1_ibm-1
            Zc_ibm(i)=Z_ibm(i+1)-0.5d0*dx1_ibm
        enddo

        do k=1,n3_ibm
            X_ibm(k)=dfloat(k-1)*dx3_ibm
        enddo
        do k=1,n3_ibm-1
            Xc_ibm(k)=X_ibm(k+1)-0.5d0*dx3_ibm
        enddo

        ! Point positions in Y non-uniform directions ------------------------
        yeta=get_computational_coords(1.d0*n2_ibm, n2_ibm)
        Y_ibm(n2_ibm)=get_physical_coords(yeta)

        do j=1,n2_ibm-1

            yeta=get_computational_coords(1.d0*j, n2_ibm)
            Y_ibm(j)=get_physical_coords(yeta)

            yeta=get_computational_coords(j+0.5d0, n2_ibm)
            Yc_ibm(j)=get_physical_coords(yeta)

            if (j>1) cell_size_Y_ibm(j-1)=Y_ibm(j)-Y_ibm(j-1)
        enddo

        cell_size_Y_ibm(n2_ibm-1)=Y_ibm(n2_ibm)-Y_ibm(n2_ibm-1)

        ! *************************************************************************
        ! Defining coefficient for schemes in irregular directions-----------------
        ! *************************************************************************

        ! Y direction--------------------------------------------------------------
        allocate(Y_to_YTr_for_D1_ibm(n2_ibm))
        Y_to_YTr_for_D1_ibm=0.d0
        allocate(Yc_to_YcTr_for_D1_ibm(n2_ibm-1))
        Yc_to_YcTr_for_D1_ibm=0.d0

        allocate(Y_to_YTr_for_D2_ibm(n2_ibm, 2))
        Y_to_YTr_for_D2_ibm=0.d0
        allocate(Yc_to_YcTr_for_D2_ibm(n2_ibm-1, 2))
        Yc_to_YcTr_for_D2_ibm=0.d0

        ! 1st derivative
        do j = 1, n2_ibm
            Y_to_YTr_for_D1_ibm(j)=DYonDy(get_computational_coords(j*1.d0, n2_ibm))
        end do

        do j = 1, n2_ibm-1
            Yc_to_YcTr_for_D1_ibm(j)=DYonDy(get_computational_coords(j+0.5d0, n2_ibm))
        end do

            ! 2nd derivative

        do j = 2, n2_ibm-1
            Y_to_YTr_for_D2_ibm(j,1)=D2YonDy2(get_computational_coords(j*1.d0, n2_ibm))
            Y_to_YTr_for_D2_ibm(j,2)=Y_to_YTr_for_D1_ibm(j)**2
        end do

        do j = 1, n2_ibm-1
            Yc_to_YcTr_for_D2_ibm(j,1)=D2YonDy2(get_computational_coords(j+0.5d0, n2_ibm))
            Yc_to_YcTr_for_D2_ibm(j,2)=Yc_to_YcTr_for_D1_ibm(j)**2
        end do

        ! Write mesh in Log path ____________________________________________________

        if (nrank==0) then

            mesh_generator_path='Log/Mesh_generator/Y/'

            mesh_generator_path=trim(log_path)//"Mesh_generator/"
            mesh_generator_path_y=trim(mesh_generator_path)//"Y/"

            open(YMESH_UNIT,file=trim(mesh_generator_path_y)//'ymesh_ibm.out')

            write(YMESH_UNIT,*) 'j, y, dy, yc, a, b, c, a_c, b_c, c_c'
            do j=1,n2_ibm-1
                dy=Y_ibm(j+1)-Y_ibm(j)
                write(YMESH_UNIT,*) j, Y_ibm(j), dy, Yc_ibm(j), Y_to_YTr_for_D1_ibm(j), Y_to_YTr_for_D2_ibm(j,2), Y_to_YTr_for_D2_ibm(j,1), Yc_to_YcTr_for_D1_ibm(j), Yc_to_YcTr_for_D2_ibm(j,2), Yc_to_YcTr_for_D2_ibm(j,1)
            end do

            close(YMESH_UNIT)

        end if

        if (nrank.eq.0) write(*,*) 'size fine mesh for second order IBM', n1_ibm, n2_ibm, n3_ibm

    end subroutine generate_mesh_ibm

    subroutine extrapolate_mask_fine_mesh_ibm()

        use IBM_data
        use decomp_2d

        implicit none
        integer :: i,j,k

        do i=xstart(1), xend(1)
            do j=xstart(2), xend(2)
                do k=xstart(3), xend(3)

                    IBM_mask1_x_ibm(2*i-1,2*j-1,2*k-1)                                  = IBM_mask1_x(i,j,k)
                    IBM_mask1_x_ibm(min(2*i,n1_ibm),2*j-1,2*k-1)                        = IBM_mask1_x(i,j,k)
                    IBM_mask1_x_ibm(2*i-1,min(2*j,n2_ibm),2*k-1)                        = IBM_mask1_x(i,j,k)
                    IBM_mask1_x_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),2*k-1)              = IBM_mask1_x(i,j,k)
                    IBM_mask1_x_ibm(2*i-1,2*j-1,min(2*k,n3_ibm))                        = IBM_mask1_x(i,j,k)
                    IBM_mask1_x_ibm(2*i-1,min(2*j,n2_ibm),min(2*k,n3_ibm))              = IBM_mask1_x(i,j,k)
                    IBM_mask1_x_ibm(min(2*i,n1_ibm),2*j-1,min(2*k,n3_ibm))              = IBM_mask1_x(i,j,k)
                    IBM_mask1_x_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),min(2*k,n3_ibm))    = IBM_mask1_x(i,j,k)

                    IBM_mask2_x_ibm(2*i-1,2*j-1,2*k-1)                                  = IBM_mask2_x(i,j,k)
                    IBM_mask2_x_ibm(min(2*i,n1_ibm),2*j-1,2*k-1)                        = IBM_mask2_x(i,j,k)
                    IBM_mask2_x_ibm(2*i-1,min(2*j,n2_ibm),2*k-1)                        = IBM_mask2_x(i,j,k)
                    IBM_mask2_x_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),2*k-1)              = IBM_mask2_x(i,j,k)
                    IBM_mask2_x_ibm(2*i-1,2*j-1,min(2*k,n3_ibm))                        = IBM_mask2_x(i,j,k)
                    IBM_mask2_x_ibm(2*i-1,min(2*j,n2_ibm),min(2*k,n3_ibm))              = IBM_mask2_x(i,j,k)
                    IBM_mask2_x_ibm(min(2*i,n1_ibm),2*j-1,min(2*k,n3_ibm))              = IBM_mask2_x(i,j,k)
                    IBM_mask2_x_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),min(2*k,n3_ibm))    = IBM_mask2_x(i,j,k)

                    IBM_mask3_x_ibm(2*i-1,2*j-1,2*k-1)                                  = IBM_mask3_x(i,j,k)
                    IBM_mask3_x_ibm(min(2*i,n1_ibm),2*j-1,2*k-1)                        = IBM_mask3_x(i,j,k)
                    IBM_mask3_x_ibm(2*i-1,min(2*j,n2_ibm),2*k-1)                        = IBM_mask3_x(i,j,k)
                    IBM_mask3_x_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),2*k-1)              = IBM_mask3_x(i,j,k)
                    IBM_mask3_x_ibm(2*i-1,2*j-1,min(2*k,n3_ibm))                        = IBM_mask3_x(i,j,k)
                    IBM_mask3_x_ibm(2*i-1,min(2*j,n2_ibm),min(2*k,n3_ibm))              = IBM_mask3_x(i,j,k)
                    IBM_mask3_x_ibm(min(2*i,n1_ibm),2*j-1,min(2*k,n3_ibm))              = IBM_mask3_x(i,j,k)
                    IBM_mask3_x_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),min(2*k,n3_ibm))    = IBM_mask3_x(i,j,k)

                enddo
            enddo
        enddo


    end subroutine extrapolate_mask_fine_mesh_ibm

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
        real*8, dimension(:), allocatable                   :: X1_ibm, X2_ibm, X3_ibm, X1c_ibm, X2c_ibm, X3c_ibm
        integer     :: i,j,k
        character(200)                  :: snapshot_file

        ! ALLOCATE GLOBAL DATA
        allocate(i_start_obj_q1(number_of_objects))
        i_start_obj_q1=0
        allocate(i_end_obj_q1(number_of_objects))
        i_end_obj_q1=0
        allocate(j_start_obj_q1(number_of_objects))
        j_start_obj_q1=0
        allocate(j_end_obj_q1(number_of_objects))
        j_end_obj_q1=0
        allocate(k_start_obj_q1(number_of_objects))
        k_start_obj_q1=0
        allocate(k_end_obj_q1(number_of_objects))
        k_end_obj_q1=0

        allocate(i_start_obj_q2(number_of_objects))
        i_start_obj_q2=0
        allocate(i_end_obj_q2(number_of_objects))
        i_end_obj_q2=0
        allocate(j_start_obj_q2(number_of_objects))
        j_start_obj_q2=0
        allocate(j_end_obj_q2(number_of_objects))
        j_end_obj_q2=0
        allocate(k_start_obj_q2(number_of_objects))
        k_start_obj_q2=0
        allocate(k_end_obj_q2(number_of_objects))
        k_end_obj_q2=0

        allocate(i_start_obj_q3(number_of_objects))
        i_start_obj_q3=0
        allocate(i_end_obj_q3(number_of_objects))
        i_end_obj_q3=0
        allocate(j_start_obj_q3(number_of_objects))
        j_start_obj_q3=0
        allocate(j_end_obj_q3(number_of_objects))
        j_end_obj_q3=0
        allocate(k_start_obj_q3(number_of_objects))
        k_start_obj_q3=0
        allocate(k_end_obj_q3(number_of_objects))
        k_end_obj_q3=0

        allocate(xstart_ibm_q1(number_of_objects,3))
        xstart_ibm_q1=0
        allocate(xend_ibm_q1(number_of_objects,3))
        xend_ibm_q1=0
        allocate(xsize_ibm_q1(number_of_objects,3))
        xsize_ibm_q1=0
        allocate(ystart_ibm_q1(number_of_objects,3))
        ystart_ibm_q1=0
        allocate(yend_ibm_q1(number_of_objects,3))
        yend_ibm_q1=0
        allocate(ysize_ibm_q1(number_of_objects,3))
        ysize_ibm_q1=0
        allocate(zstart_ibm_q1(number_of_objects,3))
        zstart_ibm_q1=0
        allocate(zend_ibm_q1(number_of_objects,3))
        zend_ibm_q1=0
        allocate(zsize_ibm_q1(number_of_objects,3))
        zsize_ibm_q1=0

        allocate(object_in_current_proc_x_q1(number_of_objects))
        allocate(object_in_current_proc_y_q1(number_of_objects))
        allocate(object_in_current_proc_z_q1(number_of_objects))

        allocate(xstart_ibm_q2(number_of_objects,3))
        xstart_ibm_q2=0
        allocate(xend_ibm_q2(number_of_objects,3))
        xend_ibm_q2=0
        allocate(xsize_ibm_q2(number_of_objects,3))
        xsize_ibm_q2=0
        allocate(ystart_ibm_q2(number_of_objects,3))
        ystart_ibm_q2=0
        allocate(yend_ibm_q2(number_of_objects,3))
        yend_ibm_q2=0
        allocate(ysize_ibm_q2(number_of_objects,3))
        ysize_ibm_q2=0
        allocate(zstart_ibm_q2(number_of_objects,3))
        zstart_ibm_q2=0
        allocate(zend_ibm_q2(number_of_objects,3))
        zend_ibm_q2=0
        allocate(zsize_ibm_q2(number_of_objects,3))
        zsize_ibm_q2=0

        allocate(object_in_current_proc_x_q2(number_of_objects))
        allocate(object_in_current_proc_y_q2(number_of_objects))
        allocate(object_in_current_proc_z_q2(number_of_objects))

        allocate(xstart_ibm_q3(number_of_objects,3))
        xstart_ibm_q3=0
        allocate(xend_ibm_q3(number_of_objects,3))
        xend_ibm_q3=0
        allocate(xsize_ibm_q3(number_of_objects,3))
        xsize_ibm_q3=0
        allocate(ystart_ibm_q3(number_of_objects,3))
        ystart_ibm_q3=0
        allocate(yend_ibm_q3(number_of_objects,3))
        yend_ibm_q3=0
        allocate(ysize_ibm_q3(number_of_objects,3))
        ysize_ibm_q3=0
        allocate(zstart_ibm_q3(number_of_objects,3))
        zstart_ibm_q3=0
        allocate(zend_ibm_q3(number_of_objects,3))
        zend_ibm_q3=0
        allocate(zsize_ibm_q3(number_of_objects,3))
        zsize_ibm_q3=0

        allocate(object_in_current_proc_x_q3(number_of_objects))
        allocate(object_in_current_proc_y_q3(number_of_objects))
        allocate(object_in_current_proc_z_q3(number_of_objects))


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

        ! initialized masks
        IBM_mask1_x = 0.d0
        IBM_mask2_x = 0.d0
        IBM_mask3_x = 0.d0
        ibm_volume = 0.d0

        do i=1,number_of_objects

            ! READING OBJECT FILE
            call read_object_size(trim(obj_file_path), nb_vertex, nb_faces)
            call read_object(trim(obj_file_path), vertex, faces, nb_vertex, nb_faces)

            ! Object is positioned according the user requirements
            ! vertex(:,1)=vertex(:,1)*body_scale_x1(i)+body_x1(i)*L1
            ! ! Shift in direction 1 with L_y so that it is easier
            ! ! vertex(:,1)=vertex(:,1)*body_scale_x1+body_x1*2.d0
            ! vertex(:,2)=vertex(:,2)*body_scale_x2(i)+body_x2(i)*2.d0
            ! vertex(:,3)=vertex(:,3)*body_scale_x3(i)+body_x3(i)*L3

            ! Scale with H
            vertex(:,1)=vertex(:,1)*body_scale_x1(i)+body_x1(i)
            ! Shift in direction 1 with L_y so that it is easier
            ! vertex(:,1)=vertex(:,1)*body_scale_x1+body_x1*2.d0
            vertex(:,2)=vertex(:,2)*body_scale_x2(i)+body_x2(i)
            vertex(:,3)=vertex(:,3)*body_scale_x3(i)+body_x3(i)

            ! ! Shift in directions with H so that it is easier
            ! vertex(:,1)=vertex(:,1)*body_scale_x1(i)+body_x1(i)*1.d0
            ! vertex(:,2)=vertex(:,2)*body_scale_x2(i)+body_x2(i)*1.d0
            ! vertex(:,3)=vertex(:,3)*body_scale_x3(i)+body_x3(i)*1.d0

            call get_objet_area(vertex, nb_vertex)

            ! PERFORMING THE MASK
            call perform_mask(vertex, faces, nb_vertex, nb_faces, IBM_mask1_x, n1, &
            n2,xstart(2),xend(2), n3,xstart(3),xend(3), X1, X2c, X3c, 1.d0)

            call perform_mask(vertex, faces, nb_vertex, nb_faces, IBM_mask2_x, n1, &
            n2,xstart(2),xend(2), n3,xstart(3),xend(3), X1c, X2, X3c, 1.d0)

            call perform_mask(vertex, faces, nb_vertex, nb_faces, IBM_mask3_x, n1, &
            n2,xstart(2),xend(2), n3,xstart(3),xend(3), X1c, X2c, X3, 1.d0)

            call perform_mask(vertex, faces, nb_vertex, nb_faces, IBM_maskcc_x, n1, &
            n2,xstart(2),xend(2), n3,xstart(3),xend(3), X1c, X2c, X3c, 1.d0)

            call get_ijk_IBM(X1, X2c, X3c,i_start_obj_q1(i),i_end_obj_q1(i),j_start_obj_q1(i),j_end_obj_q1(i),k_start_obj_q1(i),k_end_obj_q1(i))
            call get_ijk_IBM(X1c, X2, X3c,i_start_obj_q2(i),i_end_obj_q2(i),j_start_obj_q2(i),j_end_obj_q2(i),k_start_obj_q2(i),k_end_obj_q2(i))
            call get_ijk_IBM(X1c, X2c, X3,i_start_obj_q3(i),i_end_obj_q3(i),j_start_obj_q3(i),j_end_obj_q3(i),k_start_obj_q3(i),k_end_obj_q3(i))

            call decomp_objet_among_procs(xstart_ibm_q1(i,:), xend_ibm_q1(i,:), &
                                      ystart_ibm_q1(i,:), yend_ibm_q1(i,:), &
                                      zstart_ibm_q1(i,:), zend_ibm_q1(i,:), &
                                      xsize_ibm_q1(i,:), ysize_ibm_q1(i,:), zsize_ibm_q1(i,:), &
                                      xstart, xend, &
                                      ystart, yend, &
                                      zstart, zend, &
                                      i_start_obj_q1(i), i_end_obj_q1(i), &
                                      j_start_obj_q1(i), j_end_obj_q1(i), &
                                      k_start_obj_q1(i), k_end_obj_q1(i), &
                                      object_in_current_proc_x_q1(i), object_in_current_proc_y_q1(i), object_in_current_proc_z_q1(i))

            call decomp_objet_among_procs(xstart_ibm_q2(i,:), xend_ibm_q2(i,:), &
                                      ystart_ibm_q2(i,:), yend_ibm_q2(i,:), &
                                      zstart_ibm_q2(i,:), zend_ibm_q2(i,:), &
                                      xsize_ibm_q2(i,:), ysize_ibm_q2(i,:), zsize_ibm_q2(i,:), &
                                      xstart, xend, &
                                      ystart, yend, &
                                      zstart, zend, &
                                      i_start_obj_q2(i), i_end_obj_q2(i), &
                                      j_start_obj_q2(i), j_end_obj_q2(i), &
                                      k_start_obj_q2(i), k_end_obj_q2(i), &
                                      object_in_current_proc_x_q2(i), object_in_current_proc_y_q2(i), object_in_current_proc_z_q2(i))

            call decomp_objet_among_procs(xstart_ibm_q3(i,:), xend_ibm_q3(i,:), &
                                      ystart_ibm_q3(i,:), yend_ibm_q3(i,:), &
                                      zstart_ibm_q3(i,:), zend_ibm_q3(i,:), &
                                      xsize_ibm_q3(i,:), ysize_ibm_q3(i,:), zsize_ibm_q3(i,:), &
                                      xstart, xend, &
                                      ystart, yend, &
                                      zstart, zend, &
                                      i_start_obj_q3(i), i_end_obj_q3(i), &
                                      j_start_obj_q3(i), j_end_obj_q3(i), &
                                      k_start_obj_q3(i), k_end_obj_q3(i), &
                                      object_in_current_proc_x_q3(i), object_in_current_proc_y_q3(i), object_in_current_proc_z_q3(i))

        enddo

        if (maxval(IBM_mask1_x).gt.1) call exit
        if (maxval(IBM_mask2_x).gt.1) call exit
        if (maxval(IBM_mask3_x).gt.1) call exit

        ibm_volume = ibm_volume / 3.d0
        if (nrank.eq.0) write(*,*) 'ibm_volume', ibm_volume
        if (nrank.eq.0) write(*,*) 100*ibm_volume/(2.d0*L1*L3), ' % of total domain'

        if (interpol_type == SECOND_ORDER_INTERPOL) then
            call generate_mesh_ibm

            ! Call the function that initialize the decomposition on the fine mesh
            call decomp2d_info_ibm

            ! Allocate all data on the fine mesh but only on the subdomain
            call allocate_data_ibm

            ! PERFORMING THE MASK
            call extrapolate_mask_fine_mesh_ibm

        endif

        call transpose_x_to_y(IBM_maskcc_x,IBM_maskcc_y)

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

    end subroutine IBM_setup

    subroutine from_coarse_to_fine_velocity(q1_save2_x,q2_save2_x,q3_save2_x)
        use IBM_data
        use physical_fields
        use DRP_IBM
        use boundaries
        use schemes3D_interface
        use mesh, only:n1,n2,n3
        use mesh, only:X,Y,Z
        use mesh, only:Xc,Yc,Zc

        implicit none

        integer :: i,j,k
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in) :: q1_save2_x,q2_save2_x,q3_save2_x
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: q1_save2_y,q2_save2_y,q3_save2_y
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: q1_save2_z,q2_save2_z,q3_save2_z

        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: q1_z_fcx_cy_ez
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: q1_y_fcx_cy_ez, q1_y_fcx_ey_cz, q1_y_fcx_ey_ez
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: q1_x_fcx_ey_cz, q1_x_fcx_ey_ez, q1_x_fcx_cy_ez, q1_x_ccx_cy_cz, q1_x_ccx_cy_ez, q1_x_ccx_ey_cz, q1_x_ccx_ey_ez

        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: q2_z_cx_fcy_ez,q2_z_ex_fcy_ez,q2_z_ex_fcy_cz
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: q2_y_ex_fcy_cz, q2_y_cx_fcy_ez, q2_y_ex_fcy_ez, q2_y_cx_ccy_cz, q2_y_ex_ccy_cz, q2_y_cx_ccy_ez, q2_y_ex_ccy_ez
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: q2_x_ex_fcy_cz

        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: q3_z_ex_cy_fcz, q3_z_cx_ey_fcz, q3_z_ex_ey_fcz, q3_z_cx_cy_ccz, q3_z_ex_cy_ccz, q3_z_cx_ey_ccz, q3_z_ex_ey_ccz
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: q3_y_ex_cy_fcz, q3_y_cx_ey_fcz, q3_y_ex_ey_fcz
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: q3_x_ex_cy_fcz

        integer :: MPI_COMM_WORLD, mpi_err


        ! Inner values (in x,y,z decomposition configuration)-----------------
        
        call transpose_x_to_y(q1_save2_x, q1_save2_y)
        call transpose_x_to_y(q2_save2_x, q2_save2_y)
        call transpose_x_to_y(q3_save2_x, q3_save2_y)

        call transpose_y_to_z(q1_save2_y, q1_save2_z)
        call transpose_y_to_z(q2_save2_y, q2_save2_z)
        call transpose_y_to_z(q3_save2_y, q3_save2_z)

        ! q1_location
        q1_x_ibm = 0.d0
        q2_y_ibm = 0.d0
        q3_z_ibm = 0.d0

        ! prepare data

        ! q1

        ! first in z-direction
        call D0s_3Dz(q1_save2_z, q1_z_fcx_cy_ez, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_Q1_BC3)

        call transpose_z_to_y(q1_z_fcx_cy_ez, q1_y_fcx_cy_ez)

        ! then, y-direction
        call D0s_3Dy(q1_save2_y, q1_y_fcx_ey_cz, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_Q1_BC2)
        call D0s_3Dy(q1_y_fcx_cy_ez, q1_y_fcx_ey_ez, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_Q1_BC2)

        call transpose_y_to_x(q1_y_fcx_ey_cz, q1_x_fcx_ey_cz)
        call transpose_y_to_x(q1_y_fcx_ey_ez, q1_x_fcx_ey_ez)
        call transpose_y_to_x(q1_y_fcx_cy_ez, q1_x_fcx_cy_ez)

        ! finally, x_direction
        call D0s_3Dx(q1_save2_x, q1_x_ccx_cy_cz, n1,xsize(2),xsize(2),xsize(3),xsize(3), NS_Q1_BC1)
        call D0s_3Dx(q1_x_fcx_cy_ez, q1_x_ccx_cy_ez, n1,xsize(2),xsize(2),xsize(3),xsize(3), NS_Q1_BC1)
        call D0s_3Dx(q1_x_fcx_ey_cz, q1_x_ccx_ey_cz, n1,xsize(2),xsize(2),xsize(3),xsize(3), NS_Q1_BC1)
        call D0s_3Dx(q1_x_fcx_ey_ez, q1_x_ccx_ey_ez, n1,xsize(2),xsize(2),xsize(3),xsize(3), NS_Q1_BC1)

        if (object_in_current_proc_ccm_x) then

            do i=xstart_ibm_ccm(1), min(n1-1,xend_ibm_ccm(1))
                do j=xstart_ibm_ccm(2), min(n2-1,xend_ibm_ccm(2))
                    do k=xstart_ibm_ccm(3), min(n3-1,xend_ibm_ccm(3))

                        ! Initialize all 8 fine cells in 1 coarse cell
                        q1_x_ibm(2*i-1,2*j-1,2*k-1) = q1_save2_x(i,j,k)
                        q1_x_ibm(2*i-1,min(2*j,n2_ibm),2*k-1) = q1_x_fcx_ey_cz(i,j,k)
                        q1_x_ibm(2*i-1,2*j-1,min(2*k,n3_ibm)) = q1_x_fcx_cy_ez(i,j,k)
                        q1_x_ibm(2*i-1,min(2*j,n2_ibm),min(2*k,n3_ibm)) = q1_x_fcx_ey_ez(i,j,k)

                        q1_x_ibm(min(2*i,n1_ibm),2*j-1,2*k-1) = q1_x_ccx_cy_cz(i,j,k)
                        q1_x_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),2*k-1) = q1_x_ccx_cy_ez(i,j,k)
                        q1_x_ibm(min(2*i,n1_ibm),2*j-1,min(2*k,n3_ibm)) = q1_x_ccx_ey_cz(i,j,k)
                        q1_x_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),min(2*k,n3_ibm)) = q1_x_ccx_ey_ez(i,j,k)

                    enddo
                enddo
            enddo

        endif

        ! q2

        ! first in x-direction
        call D0s_3Dx(q2_save2_x, q2_x_ex_fcy_cz, n1,xsize(2),xsize(2),xsize(3),xsize(3), NS_Q2_BC1)

        call transpose_x_to_y(q2_x_ex_fcy_cz, q2_y_ex_fcy_cz)
        call transpose_y_to_z(q2_y_ex_fcy_cz, q2_z_ex_fcy_cz)

        ! then, z-direction
        call D0s_3Dz(q2_save2_z, q2_z_cx_fcy_ez, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_Q2_BC3)
        call D0s_3Dz(q2_z_ex_fcy_cz, q2_z_ex_fcy_ez, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_Q2_BC3)

        call transpose_z_to_y(q2_z_cx_fcy_ez, q2_y_cx_fcy_ez)
        call transpose_z_to_y(q2_z_ex_fcy_ez, q2_y_ex_fcy_ez)

        ! finally, y_direction
        call D0s_3Dy(q2_save2_y, q2_y_cx_ccy_cz, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_Q1_BC2)
        call D0s_3Dy(q2_y_ex_fcy_cz, q2_y_ex_ccy_cz, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_Q1_BC2)
        call D0s_3Dy(q2_y_cx_fcy_ez, q2_y_cx_ccy_ez, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_Q1_BC2)
        call D0s_3Dy(q2_y_ex_fcy_ez, q2_y_ex_ccy_ez, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_Q1_BC2)

        if (object_in_current_proc_ccm_y) then

            do i=ystart_ibm_ccm(1), min(n1-1,yend_ibm_ccm(1))
                do j=ystart_ibm_ccm(2), min(n2-1,yend_ibm_ccm(2))
                    do k=ystart_ibm_ccm(3), min(n3-1,yend_ibm_ccm(3))

                        ! Initialize all 8 fine cells in 1 coarse cell
                        q2_y_ibm(2*i-1,2*j-1,2*k-1) = q2_save2_y(i,j,k)
                        q2_y_ibm(min(2*i,n1_ibm),2*j-1,2*k-1) = q2_y_ex_fcy_cz(i,j,k)
                        q2_y_ibm(2*i-1,2*j-1,min(2*k,n3_ibm)) = q2_y_cx_fcy_ez(i,j,k)
                        q2_y_ibm(min(2*i,n1_ibm),2*j-1,min(2*k,n3_ibm)) = q2_y_ex_fcy_ez(i,j,k)

                        q2_y_ibm(2*i-1,min(2*j,n2_ibm),2*k-1) = q2_y_cx_ccy_cz(i,j,k)
                        q2_y_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),2*k-1) = q2_y_ex_ccy_cz(i,j,k)
                        q2_y_ibm(2*i-1,min(2*j,n2_ibm),min(2*k,n3_ibm)) = q2_y_cx_ccy_ez(i,j,k)
                        q2_y_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),min(2*k,n3_ibm)) = q2_y_ex_ccy_ez(i,j,k)

                    enddo
                enddo
            enddo

        endif

        ! q3

        ! first in x-direction
        call D0s_3Dx(q3_save2_x, q3_x_ex_cy_fcz, n1,xsize(2),xsize(2),xsize(3),xsize(3), NS_Q3_BC1)

        call transpose_x_to_y(q3_x_ex_cy_fcz, q3_y_ex_cy_fcz)

        ! then, y-direction
        call D0s_3Dy(q3_save2_y, q3_y_cx_ey_fcz, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_Q3_BC2)
        call D0s_3Dy(q3_y_ex_cy_fcz, q3_y_ex_ey_fcz, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_Q3_BC2)

        call transpose_y_to_z(q3_y_ex_cy_fcz, q3_z_ex_cy_fcz)
        call transpose_y_to_z(q3_y_cx_ey_fcz, q3_z_cx_ey_fcz)
        call transpose_y_to_z(q3_y_ex_ey_fcz, q3_z_ex_ey_fcz)

        ! finally, z_direction
        call D0s_3Dz(q3_save2_z, q3_z_cx_cy_ccz, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_Q3_BC3)
        call D0s_3Dz(q3_z_ex_cy_fcz, q3_z_ex_cy_ccz, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_Q3_BC3)
        call D0s_3Dz(q3_z_cx_ey_fcz, q3_z_cx_ey_ccz, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_Q3_BC3)
        call D0s_3Dz(q3_z_ex_ey_fcz, q3_z_ex_ey_ccz, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_Q3_BC3)

        if (object_in_current_proc_ccm_z) then

            do i=zstart_ibm_ccm(1), min(n1-1,zend_ibm_ccm(1))
                do j=zstart_ibm_ccm(2), min(n2-1,zend_ibm_ccm(2))
                    do k=zstart_ibm_ccm(3), min(n3-1,zend_ibm_ccm(3))

                        ! Initialize all 8 fine cells in 1 coarse cell
                        q3_z_ibm(2*i-1,2*j-1,2*k-1) = q3_save2_z(i,j,k)
                        q3_z_ibm(2*i-1,min(2*j,n2_ibm),2*k-1) = q3_z_cx_ey_fcz(i,j,k)
                        q3_z_ibm(min(2*i,n1_ibm),2*j-1,2*k-1) = q3_z_ex_cy_fcz(i,j,k)
                        q3_z_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),2*k-1) = q3_z_ex_ey_fcz(i,j,k)

                        q3_z_ibm(min(2*i,n1_ibm),2*j-1,min(2*k,n3_ibm)) = q3_z_ex_cy_ccz(i,j,k)
                        q3_z_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),min(2*k,n3_ibm)) = q3_z_ex_ey_ccz(i,j,k)
                        q3_z_ibm(2*i-1,2*j-1,min(2*k,n3_ibm)) = q3_z_cx_cy_ccz(i,j,k)
                        q3_z_ibm(2*i-1,min(2*j,n2_ibm),min(2*k,n3_ibm)) = q3_z_cx_ey_ccz(i,j,k)

                    enddo
                enddo
            enddo

        endif

        ! q2_location
        q2_x_ibm = 0.d0
        q2_z_ibm = 0.d0

        q1_y_ibm = 0.d0
        q1_z_ibm = 0.d0

        call transpose_y_to_x(q2_y_ibm, q2_x_ibm, decomp_fine)
        call transpose_y_to_z(q2_y_ibm, q2_z_ibm, decomp_fine)

        call transpose_x_to_y(q1_x_ibm, q1_y_ibm, decomp_fine)
        call transpose_y_to_z(q1_y_ibm, q1_z_ibm, decomp_fine)

        ! q3_location
        q3_y_ibm = 0.d0
        q3_x_ibm = 0.d0

        call transpose_z_to_y(q3_z_ibm, q3_y_ibm, decomp_fine)
        call transpose_y_to_x(q3_y_ibm, q3_x_ibm, decomp_fine)

    end subroutine from_coarse_to_fine_velocity

    subroutine force_temperature(sca_y, sz, sca_term, deltaT)

        implicit none
        integer, dimension(3), intent(in)           :: sz
        real*8, intent(in)                          :: deltaT
        real*8, dimension(:,:,:), intent(in)        :: sca_y
        real*8, dimension(:,:,:), intent(out)       :: sca_term

        select case (interpol_type)

            case (IBM_INTERPOL_NONE, SECOND_ORDER_INTERPOL)
                call no_interpol

            case default

        end select

    contains

        subroutine no_interpol() !sca_term correspond au terme (sca^(n+1) - sca^n )/dt
            use mesh, only: n2
            implicit none
            integer ::   i,j,k

            do k=1, sz(3)
                do i=1, sz(1)
                    do j=1, n2/2
                        sca_term(i,j,k) = ( -deltaT - sca_y(i,j,k) )
                    enddo


                    do j=n2/2+1,n2
                        sca_term(i,j,k) = ( deltaT - sca_y(i,j,k) )
                    enddo
                enddo
            enddo

        end subroutine no_interpol

    end subroutine force_temperature

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

    ! return bounds of the object, depending on the location of the variable
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

    end subroutine get_ijk_IBM

    ! return bounds of the object, depending on the location of the variable
    subroutine get_n_start_IBM(mask_ibm,n,n_IBM_start,n_objects)

        use mesh, only:n1,n2,n3
        use snapshot_writer
        use COMMON_workspace_view

        implicit none
        real*8, dimension(:), intent(in)                        :: mask_ibm
        integer                                                 :: index
        integer, intent(in)                                     :: n
        integer, intent(inout)                                  :: n_objects
        integer, dimension(number_of_objects), intent(inout)    :: n_IBM_start

        !if (dir=='X') then
        !    allocate(array_1D(n1))
        !    array_1D = mask_ibm(:,j0,k0)
        !elseif (dir=='Z') then
        !    allocate(array_1D(n3))
        !    array_1D = mask_ibm(i0,j0,:)
        !endif

        ! The idea is to go through the probe, and check the rising edges of an object
        ! To do so, we check the difference of two consecutives nodes for the ibm_mask
        ! This difference is equal to 0 in the flowfield and in the object
        ! For a rising edge, it is equal to 1 and for a falling edge, it is equal to -1

        n_objects = 1

        if (mask_ibm(1).gt.(0.d0)) then
            n_IBM_start(n_objects) = 1
            n_objects = n_objects + 1
        endif

        do index = 2, n
            if ((mask_ibm(index)-mask_ibm(index-1)).gt.(0.d0)) then
                n_IBM_start(n_objects) = index
                n_objects = n_objects + 1
            endif
        end do

        ! because n_objects was initialized to 1, we have to remove 1
        n_objects = n_objects - 1

    end subroutine get_n_start_IBM

    ! return the local decomposition of the IBM object with its own meshes
    subroutine decomp_objet_among_procs(xstart_ibm, xend_ibm, ystart_ibm, yend_ibm, zstart_ibm, zend_ibm, & 
                                        xsize_ibm, ysize_ibm, zsize_ibm, &
                                        xstart_global, xend_global, ystart_global, yend_global, zstart_global, zend_global, &
                                        i_IBM_start, i_IBM_end, j_IBM_start, j_IBM_end, k_IBM_start, k_IBM_end, &
                                        is_object_in_proc_x, is_object_in_proc_y, is_object_in_proc_z)

        implicit none

        integer, intent(in) :: i_IBM_start, i_IBM_end, j_IBM_start, j_IBM_end, k_IBM_start, k_IBM_end
        integer, dimension(3), intent(in) :: xstart_global, xend_global, ystart_global, yend_global, zstart_global, zend_global
        integer, dimension(3), intent(out) :: xstart_ibm, xend_ibm, ystart_ibm, yend_ibm, zstart_ibm, zend_ibm
        integer, dimension(3), intent(out) :: xsize_ibm, ysize_ibm, zsize_ibm

        logical :: is_object_in_proc_x, is_object_in_proc_y, is_object_in_proc_z

        is_object_in_proc_x=.true.
        is_object_in_proc_y=.true.
        is_object_in_proc_z=.true.

        xstart_ibm(1) = i_IBM_start
        xend_ibm(1) = i_IBM_end

        ystart_ibm(2) = j_IBM_start
        yend_ibm(2) = j_IBM_end

        zstart_ibm(3) = k_IBM_start
        zend_ibm(3) = k_IBM_end

        ! x - direction
        if ((i_IBM_end-i_IBM_start).gt.0) then
            ! the objet is placed on "to the inlet" x_direction

            if ((yend_global(1).gt.i_IBM_start).and.(ystart_global(1).lt.i_IBM_end)) then 
                ystart_ibm(1) = max(i_IBM_start,ystart_global(1))
                yend_ibm(1) = min(i_IBM_end,yend_global(1))
            else
                ystart_ibm(1) = 0
                yend_ibm(1) = 0
                is_object_in_proc_y=.false.
            endif

            if ((zend_global(1).gt.i_IBM_start).and.(zstart_global(1).lt.i_IBM_end)) then 
                zstart_ibm(1) = max(i_IBM_start,zstart_global(1))
                zend_ibm(1) = min(i_IBM_end,zend_global(1))
            else
                zstart_ibm(1) = 0
                zend_ibm(1) = 0
                is_object_in_proc_z=.false.
            endif

        ! else
        !     ! the objet is placed on "to the outlet" x_direction

        !     if (yend_global(1).lt.i_IBM_start) then 
        !         ystart_ibm(1) = 0
        !         yend_ibm(1) = 0
        !         is_object_in_proc=.false.
        !     else
        !         ! xstart_ibm(2) = max(j_IBM_start,xstart_global(2))
        !         ! xend_ibm(2) = min(j_IBM_end,xend_global(2))
        !     endif

        !     if (zend_global(1).lt.j_IBM_start) then 
        !         zstart_ibm(1) = 0
        !         zend_ibm(1) = 0
        !         is_object_in_proc=.false.
        !     else
        !         ! xstart_ibm(2) = max(j_IBM_start,xstart_global(2))
        !         ! xend_ibm(2) = min(j_IBM_end,xend_global(2))
        !     endif
        endif

        ! y - direction
        if ((j_IBM_end-j_IBM_start).gt.0) then
            ! the objet is placed on the lower wall

            if ((xend_global(2).gt.j_IBM_start).and.(xstart_global(2).lt.j_IBM_end)) then 
                xstart_ibm(2) = max(j_IBM_start,xstart_global(2))
                xend_ibm(2) = min(j_IBM_end,xend_global(2))
            else
                xstart_ibm(2) = 0
                xend_ibm(2) = 0
                is_object_in_proc_x=.false.
            endif

            if ((zend_global(2).gt.j_IBM_start).and.(zstart_global(2).lt.j_IBM_end)) then 
                zstart_ibm(2) = max(j_IBM_start,zstart_global(2))
                zend_ibm(2) = min(j_IBM_end,zend_global(2))
            else
                zstart_ibm(2) = 0
                zend_ibm(2) = 0
                is_object_in_proc_z=.false.
            endif

        ! else
        !     ! the object is placed on the upper wall

        !     if (xend_global(2).lt.j_IBM_start) then 
        !         xstart_ibm(2) = 0
        !         xend_ibm(2) = 0
        !         is_object_in_proc=.false.
        !     else
        !         ! xstart_ibm(2) = max(j_IBM_start,xstart_global(2))
        !         ! xend_ibm(2) = min(j_IBM_end,xend_global(2))
        !     endif

        !     if (zend_global(2).lt.j_IBM_start) then 
        !         zstart_ibm(2) = 0
        !         zend_ibm(2) = 0
        !         is_object_in_proc=.false.
        !     else
        !         ! xstart_ibm(2) = max(j_IBM_start,xstart_global(2))
        !         ! xend_ibm(2) = min(j_IBM_end,xend_global(2))
        !     endif
        endif

        ! z - direction
        if ((k_IBM_end-k_IBM_start).gt.0) then
            ! the objet is placed on "to the back" z_direction

            if ((xend_global(3).gt.k_IBM_start).and.(xstart_global(3).lt.k_IBM_end)) then 
                xstart_ibm(3) = max(k_IBM_start,xstart_global(3))
                xend_ibm(3) = min(k_IBM_end,xend_global(3))
            else
                xstart_ibm(3) = 0
                xend_ibm(3) = 0
                is_object_in_proc_x=.false.
            endif

            if ((yend_global(3).gt.k_IBM_start).and.(ystart_global(3).lt.k_IBM_end)) then 
                ystart_ibm(3) = max(k_IBM_start,ystart_global(3))
                yend_ibm(3) = min(k_IBM_end,yend_global(3))
            else
                ystart_ibm(3) = 0
                yend_ibm(3) = 0
                is_object_in_proc_y=.false.
            endif

        ! else
        !     ! the object is placed on "to the front" z_direction

        !     if (xend_global(3).lt.k_IBM_start) then 
        !         xstart_ibm(3) = 0
        !         xend_ibm(3) = 0
        !         is_object_in_proc=.false.
        !     else
        !         ! xstart_ibm(2) = max(k_IBM_start,xstart_global(2))
        !         ! xend_ibm(2) = min(k_IBM_end,xend_global(2))
        !     endif

        !     if (yend_global(3).lt.k_IBM_start) then 
        !         ystart_ibm(3) = 0
        !         yend_ibm(3) = 0
        !         is_object_in_proc=.false.
        !     else
        !         ! ystart_ibm(2) = max(k_IBM_start,ystart_global(2))
        !         ! yend_ibm(2) = min(k_IBM_end,yend_global(2))
        !     endif
        endif

        ! now, fill sizes ...
        xsize_ibm =  xend_ibm - xstart_ibm + 1
        ysize_ibm =  yend_ibm - ystart_ibm + 1
        zsize_ibm =  zend_ibm - zstart_ibm + 1

    end subroutine decomp_objet_among_procs 

    subroutine init_decomp_fine()

        use decomp_2d
        use mesh, only : n1, n2, n3

        implicit none

        integer :: prow, pcol
        TYPE(DECOMP_INFO) :: decomp_coarse, decomp_ref
        integer :: i

        ! initialize decomp fine object
        ! Then correct it by using part of the 2decompfft library
        ! The idea here, is to get twice the number of elements of any direction on all procs
        ! Somehow, using decomp_info_init is not doing exactly this
        ! So to control everything, it has to be done "by hand"
        call decomp_info_init(n1_ibm,n2_ibm,n3_ibm,decomp_fine)
        call decomp_info_init(n1,n2,n3,decomp_coarse)

        i=1

        ! x config
        decomp_fine%xst(i) = 2*xstart(i) - 1
        decomp_fine%xen(i) = min(2*xend(i),n1_ibm)
        decomp_fine%xsz(i) = decomp_fine%xen(i) - decomp_fine%xst(i) + 1

        ! y config
        decomp_fine%yst(i) = 2*ystart(i) - 1
        decomp_fine%yen(i) = min(2*yend(i),n1_ibm)
        decomp_fine%ysz(i) = decomp_fine%yen(i) - decomp_fine%yst(i) + 1

        ! z config
        decomp_fine%zst(i) = 2*zstart(i) - 1
        decomp_fine%zen(i) = min(2*zend(i),n1_ibm)
        decomp_fine%zsz(i) = decomp_fine%zen(i) - decomp_fine%zst(i) + 1


        i=2

        ! x config
        decomp_fine%xst(i) = 2*xstart(i) - 1
        decomp_fine%xen(i) = min(2*xend(i),n2_ibm)
        decomp_fine%xsz(i) = decomp_fine%xen(i) - decomp_fine%xst(i) + 1

        ! y config
        decomp_fine%yst(i) = 2*ystart(i) - 1
        decomp_fine%yen(i) = min(2*yend(i),n2_ibm)
        decomp_fine%ysz(i) = decomp_fine%yen(i) - decomp_fine%yst(i) + 1

        ! z config
        decomp_fine%zst(i) = 2*zstart(i) - 1
        decomp_fine%zen(i) = min(2*zend(i),n2_ibm)
        decomp_fine%zsz(i) = decomp_fine%zen(i) - decomp_fine%zst(i) + 1


        i=3

        ! x config
        decomp_fine%xst(i) = 2*xstart(i) - 1
        decomp_fine%xen(i) = min(2*xend(i),n3_ibm)
        decomp_fine%xsz(i) = decomp_fine%xen(i) - decomp_fine%xst(i) + 1

        ! y config
        decomp_fine%yst(i) = 2*ystart(i) - 1
        decomp_fine%yen(i) = min(2*yend(i),n3_ibm)
        decomp_fine%ysz(i) = decomp_fine%yen(i) - decomp_fine%yst(i) + 1

        ! z config
        decomp_fine%zst(i) = 2*zstart(i) - 1
        decomp_fine%zen(i) = min(2*zend(i),n3_ibm)
        decomp_fine%zsz(i) = decomp_fine%zen(i) - decomp_fine%zst(i) + 1

        prow = size(decomp_fine%x1dist)
        pcol = size(decomp_fine%y2dist)

        do i=0,prow-2

            decomp_fine%x1dist(i) = 2*decomp_coarse%x1dist(i)
            decomp_fine%y1dist(i) = 2*decomp_coarse%y1dist(i)

        enddo

        decomp_fine%x1dist(prow-1) = 2*decomp_coarse%x1dist(prow-1) - 1
        decomp_fine%y1dist(prow-1) = 2*decomp_coarse%y1dist(prow-1) - 1

        if (pcol.gt.1) then

            do i=0,pcol-2

                decomp_fine%y2dist(i) = 2*decomp_coarse%y2dist(i)
                decomp_fine%z2dist(i) = 2*decomp_coarse%z2dist(i)

            enddo

        endif

        decomp_fine%y2dist(pcol-1) = 2*decomp_coarse%y2dist(pcol-1) - 1
        decomp_fine%z2dist(pcol-1) = 2*decomp_coarse%z2dist(pcol-1) - 1

        do i=0, prow-1
           decomp_fine%x1cnts(i) = decomp_fine%x1dist(i)*decomp_fine%xsz(2)*decomp_fine%xsz(3)
           decomp_fine%y1cnts(i) = decomp_fine%ysz(1)*decomp_fine%y1dist(i)*decomp_fine%ysz(3)
           if (i==0) then
              decomp_fine%x1disp(i) = 0  ! displacement is 0-based index
              decomp_fine%y1disp(i) = 0
           else
              decomp_fine%x1disp(i) = decomp_fine%x1disp(i-1) + decomp_fine%x1cnts(i-1)
              decomp_fine%y1disp(i) = decomp_fine%y1disp(i-1) + decomp_fine%y1cnts(i-1)
           end if
        end do
        
        do i=0, pcol-1
           decomp_fine%y2cnts(i) = decomp_fine%ysz(1)*decomp_fine%y2dist(i)*decomp_fine%ysz(3)
           decomp_fine%z2cnts(i) = decomp_fine%zsz(1)*decomp_fine%zsz(2)*decomp_fine%z2dist(i)
           if (i==0) then
              decomp_fine%y2disp(i) = 0  ! displacement is 0-based index
              decomp_fine%z2disp(i) = 0
           else
              decomp_fine%y2disp(i) = decomp_fine%y2disp(i-1) + decomp_fine%y2cnts(i-1)
              decomp_fine%z2disp(i) = decomp_fine%z2disp(i-1) + decomp_fine%z2cnts(i-1)
           end if
        end do

        ! for X <=> Y transposes
        decomp_fine%x1count = maxval(decomp_fine%x1dist) * &
             maxval(decomp_fine%y1dist) * decomp_fine%xsz(3)
        decomp_fine%y1count = decomp_fine%x1count
        ! for Y <=> Z transposes
        decomp_fine%y2count = maxval(decomp_fine%y2dist) * &
             maxval(decomp_fine%z2dist) * decomp_fine%zsz(1)
        decomp_fine%z2count = decomp_fine%y2count

        ! make a copy of the decomposition information associated with the
        ! default global size in these global variables so applications can
        ! use them to create data structures 
        xstart_fine = decomp_fine%xst
        ystart_fine = decomp_fine%yst
        zstart_fine = decomp_fine%zst
        xend_fine   = decomp_fine%xen
        yend_fine   = decomp_fine%yen
        zend_fine   = decomp_fine%zen
        xsize_fine  = decomp_fine%xsz
        ysize_fine  = decomp_fine%ysz
        zsize_fine  = decomp_fine%zsz

    end subroutine init_decomp_fine

    ! actually generate all 2D decomposition information for the fine mesh
    ! This comes from the subroutine decomp_2d_init of the file decomp_2d.f90
    subroutine decomp2d_info_ibm()
        use decomp_2d
        use IBM_data
        use mesh, only:X,Xc,Y,Yc,Z,Zc
        use mesh, only:n1, n2, n3
        implicit none

        integer :: i_IBM_start, i_IBM_end, j_IBM_start, j_IBM_end, k_IBM_start, k_IBM_end
        integer :: i_IBM_start_tmp, i_IBM_end_tmp, j_IBM_start_tmp, j_IBM_end_tmp, k_IBM_start_tmp, k_IBM_end_tmp

        ! staring/ending index and size of data held by current processor
        ! duplicate 'decomp_fine', needed by apps to define data structure
        call init_decomp_fine

        ! get the bounding box, where the refinement is done
        i_IBM_start_tmp=n1_ibm; i_IBM_end_tmp=1;
        j_IBM_start_tmp=n2_ibm; j_IBM_end_tmp=1;
        k_IBM_start_tmp=n3_ibm; k_IBM_end_tmp=1;

        ! for q1 location
        call get_ijk_IBM(Z,Yc,Xc,i_IBM_start,i_IBM_end,j_IBM_start,j_IBM_end,k_IBM_start,k_IBM_end)
        if (i_IBM_start.lt.i_IBM_start_tmp) i_IBM_start_tmp=i_IBM_start
        if (j_IBM_start.lt.j_IBM_start_tmp) j_IBM_start_tmp=j_IBM_start
        if (k_IBM_start.lt.k_IBM_start_tmp) k_IBM_start_tmp=k_IBM_start
        if (i_IBM_end.gt.i_IBM_end_tmp) i_IBM_end_tmp=i_IBM_end
        if (j_IBM_end.gt.j_IBM_end_tmp) j_IBM_end_tmp=j_IBM_end
        if (k_IBM_end.gt.k_IBM_end_tmp) k_IBM_end_tmp=k_IBM_end

        ! for q2 location
        call get_ijk_IBM(Zc,Y,Xc,i_IBM_start,i_IBM_end,j_IBM_start,j_IBM_end,k_IBM_start,k_IBM_end)
        if (i_IBM_start.lt.i_IBM_start_tmp) i_IBM_start_tmp=i_IBM_start
        if (j_IBM_start.lt.j_IBM_start_tmp) j_IBM_start_tmp=j_IBM_start
        if (k_IBM_start.lt.k_IBM_start_tmp) k_IBM_start_tmp=k_IBM_start
        if (i_IBM_end.gt.i_IBM_end_tmp) i_IBM_end_tmp=i_IBM_end
        if (j_IBM_end.gt.j_IBM_end_tmp) j_IBM_end_tmp=j_IBM_end
        if (k_IBM_end.gt.k_IBM_end_tmp) k_IBM_end_tmp=k_IBM_end

        ! for q3 location
        call get_ijk_IBM(Zc,Yc,X,i_IBM_start,i_IBM_end,j_IBM_start,j_IBM_end,k_IBM_start,k_IBM_end)
        if (i_IBM_start.lt.i_IBM_start_tmp) i_IBM_start_tmp=i_IBM_start
        if (j_IBM_start.lt.j_IBM_start_tmp) j_IBM_start_tmp=j_IBM_start
        if (k_IBM_start.lt.k_IBM_start_tmp) k_IBM_start_tmp=k_IBM_start
        if (i_IBM_end.gt.i_IBM_end_tmp) i_IBM_end_tmp=i_IBM_end
        if (j_IBM_end.gt.j_IBM_end_tmp) j_IBM_end_tmp=j_IBM_end
        if (k_IBM_end.gt.k_IBM_end_tmp) k_IBM_end_tmp=k_IBM_end

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DEBUG
        ! physical bounds for the coarse mesh
        i_start_ibm_2nd_cp=max(1,i_IBM_start_tmp-7); i_end_ibm_2nd_cp=min(n1-1,i_IBM_end_tmp+7);
        j_start_ibm_2nd_cp=max(1,j_IBM_start_tmp-7); j_end_ibm_2nd_cp=min(n2-1,j_IBM_end_tmp+7);
        k_start_ibm_2nd_cp=max(1,k_IBM_start_tmp-7); k_end_ibm_2nd_cp=min(n3-1,k_IBM_end_tmp+7);

        ! computational bounds for the coarse mesh
        i_start_ibm_2nd_cc=max(1,i_start_ibm_2nd_cp-3); i_end_ibm_2nd_cc=min(n1-1,i_end_ibm_2nd_cp+3);
        j_start_ibm_2nd_cc=max(1,j_start_ibm_2nd_cp-3); j_end_ibm_2nd_cc=min(n2-1,j_end_ibm_2nd_cp+3);
        k_start_ibm_2nd_cc=max(1,k_start_ibm_2nd_cp-3); k_end_ibm_2nd_cc=min(n3-1,k_end_ibm_2nd_cp+3);

        ! physical bounds for the fine mesh
        i_start_ibm_2nd_fp=max(1,(2*i_IBM_start_tmp-1)-14); i_end_ibm_2nd_fp=min(n1_ibm-1,(2*i_IBM_end_tmp)+14);
        j_start_ibm_2nd_fp=max(1,(2*j_IBM_start_tmp-1)-14); j_end_ibm_2nd_fp=min(n2_ibm-1,(2*j_IBM_end_tmp)+14);
        k_start_ibm_2nd_fp=max(1,(2*k_IBM_start_tmp-1)-14); k_end_ibm_2nd_fp=min(n3_ibm-1,(2*k_IBM_end_tmp)+14);

        ! computational bounds for the fine mesh
        i_start_ibm_2nd_fc=max(1,i_start_ibm_2nd_fp-6); i_end_ibm_2nd_fc=min(n1_ibm-1,i_end_ibm_2nd_fp+6);
        j_start_ibm_2nd_fc=max(1,j_start_ibm_2nd_fp-6); j_end_ibm_2nd_fc=min(n2_ibm-1,j_end_ibm_2nd_fp+6);
        k_start_ibm_2nd_fc=max(1,k_start_ibm_2nd_fp-6); k_end_ibm_2nd_fc=min(n3_ibm-1,k_end_ibm_2nd_fp+6);

        if (nrank.eq.0) then

            write(*,*) '============================='
            write(*,*) 'COMPUTATIONAL COARSE MESH'
            write(*,*) i_start_ibm_2nd_cc, i_end_ibm_2nd_cc
            write(*,*) j_start_ibm_2nd_cc, j_end_ibm_2nd_cc
            write(*,*) k_start_ibm_2nd_cc, k_end_ibm_2nd_cc

            write(*,*) '============================='
            write(*,*) 'PHYSICAL COARSE MESH'
            write(*,*) i_start_ibm_2nd_cp, i_end_ibm_2nd_cp
            write(*,*) j_start_ibm_2nd_cp, j_end_ibm_2nd_cp
            write(*,*) k_start_ibm_2nd_cp, k_end_ibm_2nd_cp

            write(*,*) '============================='
            write(*,*) 'COMPUTATIONAL FINE MESH'
            write(*,*) i_start_ibm_2nd_fc, i_end_ibm_2nd_fc
            write(*,*) j_start_ibm_2nd_fc, j_end_ibm_2nd_fc
            write(*,*) k_start_ibm_2nd_fc, k_end_ibm_2nd_fc

            write(*,*) '============================='
            write(*,*) 'PHYSICAL FINE MESH'
            write(*,*) i_start_ibm_2nd_fp, i_end_ibm_2nd_fp
            write(*,*) j_start_ibm_2nd_fp, j_end_ibm_2nd_fp
            write(*,*) k_start_ibm_2nd_fp, k_end_ibm_2nd_fp

        endif

        call decomp_objet_among_procs(xstart_ibm_ccm, xend_ibm_ccm, &
                                      ystart_ibm_ccm, yend_ibm_ccm, &
                                      zstart_ibm_ccm, zend_ibm_ccm, &
                                      xsize_ibm_ccm, ysize_ibm_ccm, zsize_ibm_ccm, &
                                      xstart, xend, &
                                      ystart, yend, &
                                      zstart, zend, &
                                      i_start_ibm_2nd_cc, i_end_ibm_2nd_cc, &
                                      j_start_ibm_2nd_cc, j_end_ibm_2nd_cc, &
                                      k_start_ibm_2nd_cc, k_end_ibm_2nd_cc, &
                                      object_in_current_proc_ccm_x, object_in_current_proc_ccm_y, object_in_current_proc_ccm_z)

        call decomp_objet_among_procs(xstart_ibm_cpm, xend_ibm_cpm, &
                                      ystart_ibm_cpm, yend_ibm_cpm, &
                                      zstart_ibm_cpm, zend_ibm_cpm, &
                                      xsize_ibm_cpm, ysize_ibm_cpm, zsize_ibm_cpm, &
                                      xstart, xend, &
                                      ystart, yend, &
                                      zstart, zend, &
                                      i_start_ibm_2nd_cp, i_end_ibm_2nd_cp, &
                                      j_start_ibm_2nd_cp, j_end_ibm_2nd_cp, &
                                      k_start_ibm_2nd_cp, k_end_ibm_2nd_cp, &
                                      object_in_current_proc_cpm_x, object_in_current_proc_cpm_y, object_in_current_proc_cpm_z)

        call decomp_objet_among_procs(xstart_ibm_fcm, xend_ibm_fcm, &
                                      ystart_ibm_fcm, yend_ibm_fcm, &
                                      zstart_ibm_fcm, zend_ibm_fcm, &
                                      xsize_ibm_fcm, ysize_ibm_fcm, zsize_ibm_fcm, &
                                      xstart_fine, xend_fine, &
                                      ystart_fine, yend_fine, &
                                      zstart_fine, zend_fine, &
                                      i_start_ibm_2nd_fc, i_end_ibm_2nd_fc, &
                                      j_start_ibm_2nd_fc, j_end_ibm_2nd_fc, &
                                      k_start_ibm_2nd_fc, k_end_ibm_2nd_fc, &
                                      object_in_current_proc_fcm_x, object_in_current_proc_fcm_y, object_in_current_proc_fcm_z)

        call decomp_objet_among_procs(xstart_ibm_fpm, xend_ibm_fpm, &
                                      ystart_ibm_fpm, yend_ibm_fpm, &
                                      zstart_ibm_fpm, zend_ibm_fpm, &
                                      xsize_ibm_fpm, ysize_ibm_fpm, zsize_ibm_fpm, &
                                      xstart_fine, xend_fine, &
                                      ystart_fine, yend_fine, &
                                      zstart_fine, zend_fine, &
                                      i_start_ibm_2nd_fp, i_end_ibm_2nd_fp, &
                                      j_start_ibm_2nd_fp, j_end_ibm_2nd_fp, &
                                      k_start_ibm_2nd_fp, k_end_ibm_2nd_fp, &
                                      object_in_current_proc_fpm_x, object_in_current_proc_fpm_y, object_in_current_proc_fpm_z)


    end subroutine decomp2d_info_ibm

    ! actually allocate data for the all arrays related to the fine grid
    subroutine allocate_data_ibm()
        use IBM_data
        implicit none
        
        ! Inner values (in x,y,z decomposition configuration)-----------------
        
        ! q1_location

        allocate(q1_x_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
        q1_x_ibm=0.d0
        allocate(IBM_mask1_x_ibm(xstart_fine(1):xend_fine(1),xstart_fine(2):xend_fine(2),xstart_fine(3):xend_fine(3)))
        IBM_mask1_x_ibm=0.d0
        allocate(vel_term1_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
        vel_term1_ibm=0.d0

        allocate(q1_y_ibm(ystart_fine(1):yend_fine(1), ystart_fine(2):yend_fine(2), ystart_fine(3):yend_fine(3)))
        q1_y_ibm=0.d0

        allocate(q1_z_ibm(zstart_fine(1):zend_fine(1), zstart_fine(2):zend_fine(2), zstart_fine(3):zend_fine(3)))
        q1_z_ibm=0.d0

        ! Wall values allocation ---------------------------------------------
        !ATTENTION
        allocate(q1_wall10_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)))
        q1_wall10_ibm=0.d0
        allocate(q1_wall11_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)))
        q1_wall11_ibm=0.d0
        allocate(q1_wall20_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3)))
        q1_wall20_ibm=0.d0
        allocate(q1_wall21_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3)))
        q1_wall21_ibm=0.d0
        allocate(q1_wall30_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2)))
        q1_wall30_ibm=0.d0
        allocate(q1_wall31_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2)))
        q1_wall31_ibm=0.d0
        

        ! q2_location

        allocate(q2_x_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
        q2_x_ibm=0.d0
        allocate(IBM_mask2_x_ibm(xstart_fine(1):xend_fine(1),xstart_fine(2):xend_fine(2),xstart_fine(3):xend_fine(3)))
        IBM_mask2_x_ibm=0.d0
        allocate(vel_term2_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
        vel_term2_ibm=0.d0

        allocate(q2_y_ibm(ystart_fine(1):yend_fine(1), ystart_fine(2):yend_fine(2), ystart_fine(3):yend_fine(3)))
        q2_y_ibm=0.d0

        allocate(q2_z_ibm(zstart_fine(1):zend_fine(1), zstart_fine(2):zend_fine(2), zstart_fine(3):zend_fine(3)))
        q2_z_ibm=0.d0

        ! Wall values allocation ---------------------------------------------
        !ATTENTION
        allocate(q2_wall10_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)))
        q2_wall10_ibm=0.d0
        allocate(q2_wall11_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)))
        q2_wall11_ibm=0.d0
        allocate(q2_wall20_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3)))
        q2_wall20_ibm=0.d0
        allocate(q2_wall21_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3)))
        q2_wall21_ibm=0.d0
        allocate(q2_wall30_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2)))
        q2_wall30_ibm=0.d0
        allocate(q2_wall31_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2)))
        q2_wall31_ibm=0.d0
        

        ! q3_location

        allocate(q3_x_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
        q3_x_ibm=0.d0
        allocate(IBM_mask3_x_ibm(xstart_fine(1):xend_fine(1),xstart_fine(2):xend_fine(2),xstart_fine(3):xend_fine(3)))
        IBM_mask3_x_ibm=0.d0
        allocate(vel_term3_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
        vel_term3_ibm=0.d0


        allocate(q3_y_ibm(ystart_fine(1):yend_fine(1), ystart_fine(2):yend_fine(2), ystart_fine(3):yend_fine(3)))
        q3_y_ibm=0.d0

        allocate(q3_z_ibm(zstart_fine(1):zend_fine(1), zstart_fine(2):zend_fine(2), zstart_fine(3):zend_fine(3)))
        q3_z_ibm=0.d0

        ! Wall values allocation ---------------------------------------------
        !ATTENTION
        allocate(q3_wall10_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)))
        q3_wall10_ibm=0.d0
        allocate(q3_wall11_ibm(xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)))
        q3_wall11_ibm=0.d0
        allocate(q3_wall20_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3)))
        q3_wall20_ibm=0.d0
        allocate(q3_wall21_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(3):yend_ibm_fcm(3)))
        q3_wall21_ibm=0.d0
        allocate(q3_wall30_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2)))
        q3_wall30_ibm=0.d0
        allocate(q3_wall31_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2)))
        q3_wall31_ibm=0.d0

    end subroutine allocate_data_ibm

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
