module MHD_field_generator
    use decomp_2d
    use mpi
    use mesh
    use MHD_data
    use mathematical_constants
    use schemes3D_interface
    use HDF5_IO
    use irregular_derivative_coefficients

!
    implicit none
!
contains
!
!
    subroutine generate_magnets(num, B1_y, B2_y, B3_y)
        use boundaries

        implicit none

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))  :: B1_y, B2_y, B3_y
        real*8  xsize, zsize, xcenter, zcenter, ycenter, B_Factor, mu_0, pole_density, h_f, h_c
        integer :: num

        xsize   = mag_array(num)%magnet_size_x  *pi /2.d0    ! half size
        zsize   = mag_array(num)%magnet_size_z  *pi /2.d0    ! half_size
        pole_density = mag_array(num)%sigma ! -1 or +1

        !!!!! Generate defined magnet !!!!!
        xcenter = mag_array(num)%magnet_center_x *pi
        zcenter = mag_array(num)%magnet_center_z *pi
        ycenter = mag_array(num)%magnet_center_y

        call generate_1_magnet(xcenter,zcenter,ycenter,xsize,zsize,pole_density)

        !!!!! Generate closest periodic magnets !!!!
        if ((BC1==UNBOUNDED).and.(L1>zsize*2.d0)) then

            zcenter = zcenter + L1
            call generate_1_magnet(xcenter,zcenter,ycenter,xsize,zsize,pole_density)

            zcenter = zcenter - 2.d0*L1
            call generate_1_magnet(xcenter,zcenter,ycenter,xsize,zsize,pole_density)

            zcenter = zcenter + L1
        endif

        if ((BC3==UNBOUNDED).and.(L3>xsize*2.d0)) then

            xcenter = xcenter + L3
            call generate_1_magnet(xcenter,zcenter,ycenter,xsize,zsize,pole_density)

            xcenter = xcenter - 2.d0*L3
            call generate_1_magnet(xcenter,zcenter,ycenter,xsize,zsize,pole_density)

            xcenter = xcenter + L3
        endif


    contains
        subroutine generate_1_magnet(xcenter,zcenter,ycenter,xsize,zsize,pole_density)

            implicit none

            real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))  :: R1,R2,R3
            real*8, dimension(ystart(1):yend(1))  :: Ti, Tic
            real*8, dimension(ystart(3):yend(3))  :: Sk, Skc

            real*8  , intent(in)    :: xsize, zsize, xcenter, zcenter, ycenter, pole_density
            real*8                  :: B_Factor, mu_0, h_f, h_c
            integer                 :: i,j,k,m,n


            mu_0     = 1.d0
            B_Factor = pole_density/(4.d0*pi*mu_0)

    !      do k = ystart(3), min(n3-1, yend(3))
    !        do i = ystart(1), min(n1-1, yend(1))
    !           do j = 1, n2-1
            R1 = 0.d0
            R2 = 0.d0
            R3 = 0.d0

    !            enddo
    !         enddo
    !       enddo

    !      B1max=0.d0
    !      B2max=0.d0
    !      B3max=0.d0
    !      B1max_glob=0.d0
    !      B2max_glob=0.d0
    !      B3max_glob=0.d0


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! WARNING : Assuming that streamwise = 3
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !
    ! WARNING B1_y,B2_y,B3_y are calculated at each FACE CENTER.
    !
    !
    !
            do m = 0,1
                do n = 0,1

                    do i = ystart(1), yend(1)
                        do k = ystart(3), yend(3)
                            do j = 1, n2

                                h_f = Y(j) - ycenter
                                if (j/=n2) h_c = Yc(j) - ycenter

                                if(k/=n3) Skc(k) = (Xc(k)-xcenter) - (-1.d0)**dfloat(m) *xsize  !centered value
                                if(i/=n1) Tic(i) = (Zc(i)-zcenter) - (-1.d0)**dfloat(n) *zsize  !centered value
                                Sk(k) = (X(k)-xcenter) - (-1.d0)**dfloat(m) *xsize    !face value
                                Ti(i) = (Z(i)-zcenter) - (-1.d0)**dfloat(n) *zsize    !face value

                                if((j/=n2).and.(k/=n3)) R1(i,j,k) = ( Skc(k)**2.d0 + Ti(i)**2.d0 + h_c**2.d0)**0.5d0
                                if((k/=n3).and.(i/=n1)) R2(i,j,k) = ( Skc(k)**2.d0 + Tic(i)**2.d0 + h_f**2.d0)**0.5d0
                                if((j/=n2).and.(i/=n1)) R3(i,j,k) = ( Sk(k)**2.d0 + Tic(i)**2.d0 + h_c**2.d0)**0.5d0

                    !            B1_y(i,j,k) =  B1_y(i,j,k) + ( (-1.d0)**dfloat((m+n)) *dlog(R(i,j,k)-Ti(i)) )*B_Factor
                    !            B3_y(i,j,k) =  B3_y(i,j,k) + ( (-1.d0)**dfloat((m+n)) *dlog(R(i,j,k)-Sk(k)) )*B_Factor
                    !            B2_y(i,j,k) =  B2_y(i,j,k) + ( (-1.d0)**dfloat((m+n)) *datan(Sk(k)*Ti(i)/(R(i,j,k)*Yc(j))) )*B_Factor
                                if((j/=n2).and.(k/=n3)) B1_y(i,j,k) =  B1_y(i,j,k) + ( (-1.d0)**dfloat((m+n)) *dlog(R1(i,j,k)-Skc(k)) )*B_Factor
                                if((j/=n2).and.(i/=n1)) B3_y(i,j,k) =  B3_y(i,j,k) + ( (-1.d0)**dfloat((m+n)) *dlog(R3(i,j,k)-Tic(i)) )*B_Factor
                                if((k/=n3).and.(i/=n1)) B2_y(i,j,k) =  B2_y(i,j,k) + ( (-1.d0)**dfloat((m+n)) *datan(Skc(k)*Tic(i)/(R2(i,j,k)*h_f)) )*B_Factor


                            enddo !j
                        enddo !k
                    enddo !i

                enddo !n
            enddo !m

    !       do i = ystart(1), min(n1-1, yend(1))
    !        do k = ystart(3), min(n3-1, yend(3))
    !          do j = 2, n2-1
    !
    !             B1max=max(B1_y(i,j,k), B1max)
    !             B2max=max(B2_y(i,j,k), B2max)
    !             B3max=max(B3_y(i,j,k), B3max)
    !
    !
    !          enddo !j
    !         enddo !k
    !        enddo !i

    !       enddo !n
    !      enddo !m

    !
    !      call MPI_ALLREDUCE (B1max, B1max_glob, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
    !      call MPI_ALLREDUCE (B2max, B2max_glob, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
    !      call MPI_ALLREDUCE (B3max, B3max_glob, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
    !
    !       write(6,*) 'B1max =',B1max_glob
    !       write(6,*) 'B2max =',B2max_glob
    !       write(6,*) 'B3max =',B3max_glob

        end subroutine generate_1_magnet

!
    end subroutine generate_magnets


    subroutine generate_total_magnetic_field(B001_x, B002_y, B003_z, B001c_y, B002c_y, B003c_y)
        implicit none

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: B001_x, B001c_x, B1_x
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))  :: B002_y, B1_y, B2_y, B3_y
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))  :: B003_z, B003c_z, B3_z
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))  :: B001c_y, B002c_y, B003c_y
        !real*8 B01max,B02max, B03max,B01max_glob,B02max_glob, B03max_glob

        integer :: i,j,k,l,nb,mpi_err,n1e,n3e, n2e, n1s, n2s, n3s
!        real*8  :: NormB
        integer :: magnets_nb
!        character*1 coma
!        coma=char(44)
!

!      B01max=0.d0
!      B02max=0.d0
!      B03max=0.d0
!      B01max_glob=0.d0
!      B02max_glob=0.d0
!      B03max_glob=0.d0

      ! ************************************** CORR_JO **************************** !
      ! Le champ magnetique est défini sur les faces des cellules, comme la vitesse !
      ! L'initialisation est corrigée pour être effectuée jusqu'à n1, n2 et n3      !
      ! De plus chaque composante est défini pour une décomp2D dans sa direciton
      ! *************************************************************************** !

        if (layer_type==0) then
            call generate_uniform_field(B001_x, B002_y, B003_z ,         &
                                        Magnetic_field_unit_vector_x,    &
                                        Magnetic_field_unit_vector_y,    &
                                        Magnetic_field_unit_vector_z      )
        elseif(layer_type==2) then

            if ((Magnetic_field_unit_vector_x/=0.d0).or.(Magnetic_field_unit_vector_z/=0.d0))then
                call generate_layer_B1_B3(B001_x, B003_z ,              &
                                          Magnetic_field_unit_vector_x, &
                                          Magnetic_field_unit_vector_z   )
            endif

            if (Magnetic_field_unit_vector_y/=0.d0)then
                call generate_bipolar_layer_B2(B001_x, B002_y, B003_z ,     &
                                               Magnetic_field_unit_vector_y )
            endif

        elseif (layer_type==1) then

            call generate_sinusoidal_layer(B001_x, B002_y, B003_z ,     &
                                            Magnetic_field_unit_vector_x,  &
                                            Magnetic_field_unit_vector_y,  &
                                            Magnetic_field_unit_vector_z   )
        endif

        magnets_nb = 2*magnet_number_pairs  ! obligatoirement pair

      ! ****************************************** CORR_JO ***************************************** !
      ! ATTENTION modifier generate_fields_1_magnet pour calculer le champ dans chaque transposition !
      ! ******************************************************************************************** !
        B1_y = 0.d0
        B1_x = 0.d0
        B2_y = 0.d0
        B3_y = 0.d0
        B3_z = 0.d0

        do nb= 1, magnets_nb
            call generate_magnets(nb, B1_y, B2_y, B3_y)
!            do k = ystart(3), min(n3-1, yend(3))
!                do i = ystart(1), min(n1-1, yend(1))
!                    do j = 2, n2-1
!
!                          B001_y(i,j,k) =  B001_y(i,j,k) + B1_y(i,j,k)
!                          B002_y(i,j,k) =  B002_y(i,j,k) + B2_y(i,j,k)
!                          B003_y(i,j,k) =  B003_y(i,j,k) + B3_y(i,j,k)
!
!                          B01max=max(B001_y(i,j,k), B01max)
!                          B02max=max(B002_y(i,j,k), B02max)
!                          B03max=max(B003_y(i,j,k), B03max)
!
!                    enddo
!                enddo
!            enddo
        enddo

        call transpose_y_to_x(B1_y, B1_x)
        call transpose_y_to_z(B3_y, B3_z)

        !!!! ADDING magnets field to the general field !!!!!
        B001_x = B001_x + B1_x
        B002_y = B002_y + B2_y
        B003_z = B003_z + B3_z

!**************** Interpolation of B_face_centered to B_cell_centered *****************
!     Can be performed only at intialization since B (and then B_centered) is constant
!

         call D0s_3Dy(B002_y, B002c_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), Dirichlet)
         call D0s_3Dz(B003_z, B003c_z, zsize(1), zsize(1),zsize(2),zsize(2),n3, Dirichlet)
         call D0s_3Dx(B001_x, B001c_x, n1,xsize(2),xsize(2),xsize(3),xsize(3), Dirichlet)

         call transpose_x_to_y(B001c_x, B001c_y)
         call transpose_z_to_y(B003c_z, B003c_y)


!      call MPI_ALLREDUCE (B01max, B01max_glob, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
!      call MPI_ALLREDUCE (B02max, B02max_glob, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
!      call MPI_ALLREDUCE (B03max, B03max_glob, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
!
!       write(6,*) 'B01max =',B01max_glob
!       write(6,*) 'B02max =',B02max_glob
!       write(6,*) 'B03max =',B03max_glob

    !*** Check the zero divergence of the magnetic field in the doamain *** !
    !   call div
    !
    !
    !*** Output Magnetic Field *** !
    !
!************************************************************************
!******************** VISU :CHAMPS INTENSITY ****************************
!************************************************************************

    if(nrank==0)  call hdf_create_file_with_3Dmesh('Magnetic_field', Zc, Yc,Xc, "x1", "x2", "x3", n1m, n2m,n3m)

    if(nrank==0)  call hdf_addgroup('Magnetic_field', "B")

    n1s=ystart(1)
    n1e=min(yend(1),n1m)
    n2s=ystart(2)
    n2e=min(yend(2),n2m)
    n3s=ystart(3)
    n3e=min(yend(3),n3m)

    call hdf_add_3Dfield('Magnetic_field', B001c_y(n1s:n1e,n2s:n2e,n3s:n3e), "B/B1", n1m, n2m, n3m, n1s, n1e, n2s, n2e, n3s, n3e)
    call hdf_add_3Dfield('Magnetic_field', B002c_y(n1s:n1e,n2s:n2e,n3s:n3e), "B/B2", n1m, n2m, n3m, n1s, n1e, n2s, n2e, n3s, n3e)
    call hdf_add_3Dfield('Magnetic_field', B003c_y(n1s:n1e,n2s:n2e,n3s:n3e), "B/B3", n1m, n2m, n3m, n1s, n1e, n2s, n2e, n3s, n3e)


!            open(28,file='Magnetic_field.plt')
!
!
!             write(28,10)' TITLE = "Magnetic Field"'
!             write(28,10)' VARIABLES = "X","Y","Z","B1","B2","B3","NormB"'
!             write(28,11) 'ZONE T = "Mag"',coma,'I=',n3-1,coma,'J=',n2-1,coma,'K=',n1-1,coma,'DATAPACKING=POINT'
!
!    ! ************ Champs 3D  **************
!         n1e=min(n1-1, yend(1))
!         n3e=min(n3-1, yend(3))
!
!         do i = ystart(1), n1e
!          do j = 1, n2-1
!            do k = ystart(3), n3e
!
!             NormB = dsqrt(((B001c_y(i,j,k))**2 +(B002c_y(i,j,k))**2 + (B003c_y(i,j,k))**2))
!    !
!             write(28,*)Xc(k),Yc(j),Zc(i),B001c_y(i,j,k),B002c_y(i,j,k),B003c_y(i,j,k),NormB
!
!            enddo
!           enddo
!          enddo
!
!          close(28)
!
!       10 format(a)
!       11 format(a,a,a,i4,a,a,i3,a,a,i3,a,a)
!       12 format(a,a,a,i4,a,a,i3,a,a)

    contains
        subroutine smooth_heaviside(x,a,b,H)

            implicit none
            real*8  :: x,a,b, H

            H = 0.5d0 +0.5d0*tanh((x-a)/b)
        end subroutine

        subroutine generate_uniform_field(B1x,B2y,B3z,e1,e2,e3)

            implicit none

            real*8  :: e1,e2,e3, H1, H2
            integer :: n1e, n2e, n3e

            real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: B1x
            real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))  :: B2y
            real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))  :: B3z

            n3e = min(xend(3), n3-1)
            n2e = min(xend(2), n2-1)

            do j = xstart(2), n2e

                call smooth_heaviside(Yc(j),delta_B,smooth*delta_B,H1)
                call smooth_heaviside(Yc(j),L2-delta_B,smooth*delta_B,H2)

                do k = xstart(3), n3e
                      do i = 1, n1
                        B1x(i,j,k) = e1*(1-H1+H2)
                      enddo
                 enddo
            enddo

            n1e = min(zend(1), n1-1)
            n2e = min(zend(2), n2-1)

            do j = zstart(2), n2e

                 call smooth_heaviside(Yc(j),delta_B,smooth*delta_B,H1)
                 call smooth_heaviside(Yc(j),L2-delta_B,smooth*delta_B,H2)

                 do i = zstart(1), n1e
                      do k = 1, n3
                          B3z(i,j,k) = e3*(1-H1+H2)
                      enddo
                 enddo
            enddo

            n3e = min(yend(3), n3-1)
            n1e = min(yend(1), n1-1)

            do j = 1, n2

                call smooth_heaviside(Y(j),delta_B,smooth*delta_B,H1)
                call smooth_heaviside(Y(j),L2-delta_B,smooth*delta_B,H2)

                do k = ystart(3), n3e
                    do i = ystart(1), n1e
                        B2y(i,j,k) = e2*(1-H1+H2)

                        if((e1/=0.d0).and.Y(j)>1.d0)then
                            B2y(i,j,k) = -e2*(1-H1+H2)
                        endif
                    enddo
                enddo
            enddo

        end subroutine

        subroutine generate_layer_B1_B3(B1x,B3z,e1,e3)

            implicit none

            real*8  :: e1,e3, H1, H2
            integer :: n1e, n2e, n3e

            real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: B1x
            real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))  :: B3z

              n3e = min(xend(3), n3-1)
              n2e = min(xend(2), n2-1)

              do j = xstart(2), n2e
                   call smooth_heaviside(Yc(j),delta_B,smooth*delta_B,H1)
                   call smooth_heaviside(Yc(j),L2-delta_B,smooth*delta_B,H2)

                   do k = xstart(3), n3e
                        do i = 1, n1
                            B1x(i,j,k) = e1*(1-H1+H2)
                        enddo
                   enddo
              enddo

              n1e = min(zend(1), n1-1)
              n2e = min(zend(2), n2-1)

              do j = zstart(2), n2e
                   call smooth_heaviside(Yc(j),delta_B,smooth*delta_B,H1)
                   call smooth_heaviside(Yc(j),L2-delta_B,smooth*delta_B,H2)

                   do i = zstart(1), n1e
                        do k = 1, n3
                            B3z(i,j,k) = e3*(1-H1+H2)
                        enddo
                   enddo
              enddo


        end subroutine

        subroutine generate_bipolar_layer_B2(B1x,B2y,B3z,e2)

            use DNS_settings

            implicit none

            real*8  :: e2, H1, H2, H3, H4, H5
            integer :: n1e, n2e, n3e

            real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: B1x, dB1cx
            real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))  :: B2y, dB2cy
            real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))  :: B3z, dB3cz

              n3e = min(yend(3), n3-1)
              n1e = min(yend(1), n1-1)

              dB2cy = 0.d0
              dB3cz = 0.d0
              dB1cx = 0.d0

              do j = 1, n2
                  call smooth_heaviside(Y(j),delta_B,smooth*delta_B,H1)
                  call smooth_heaviside(Y(j),L2-delta_B,smooth*delta_B,H2)

                  do k = ystart(3), n3e
                    call smooth_heaviside(Xc(k),0.d0,L3*0.005d0,H4)
                    call smooth_heaviside(Xc(k),L3*0.5d0,L3*0.005d0,H3)
                    call smooth_heaviside(Xc(k),L3,L3*0.005d0,H5)
                       do i = ystart(1), n1e
                                B2y(i,j,k) = (2.d0*H4-1.d0)*(1.d0-2.d0*H3)*(1.d0-2.d0*H5)*e2*(1-H1+H2)
                        enddo
                  enddo
              enddo

                call D1s_MULT_3Dy(B2y, dB2cy, ysize(1),ysize(1),n2,ysize(3),ysize(3), dx2, .false., antisymetric, Yc_to_YcTr_for_D1)

                if(streamwise==3) then

                    call transpose_y_to_z(dB2cy,dB3cz)
                    dB3cz = - dB3cz

                    n1e = min(zend(1), n1-1)
                    n2e = min(zend(2), n2-1)

                    do j = zstart(2), n2e
                       do i = zstart(1), n1e
                            do k = 2, n3
                                B3z(i,j,k) = B3z(i,j,k-1) + dB3cz(i,j,k-1)*dx3
                            enddo
                       enddo
                    enddo

                elseif(streamwise==1) then

                    call transpose_y_to_x(dB2cy,dB1cx)
                    dB1cx = - dB1cx

                    n3e = min(xend(3), n3-1)
                    n2e = min(xend(2), n2-1)

                    do j = xstart(2), n2e
                       do k = xstart(3), n3e
                            do i = 2, n1
                                B1x(i,j,k) = B1x(i-1,j,k) + dB1cx(i-1,j,k)*dx1
                            enddo
                       enddo
                    enddo

                endif

        end subroutine

        subroutine generate_sinusoidal_layer(B1x,B2y,B3z,e1,e2,e3)

            use DNS_settings

            implicit none

            real*8  :: e1,e2,e3, H1, H2
            integer :: n1e, n2e, n3e

            real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: B1x, dB1cx
            real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))  :: B2y, dB2cy
            real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))  :: B3z, dB3cz


              dB2cy = 0.d0
              dB3cz = 0.d0
              dB1cx = 0.d0

              n3e = min(xend(3), n3-1)
              n2e = min(xend(2), n2-1)

              do j = xstart(2), n2e
                  call smooth_heaviside(Yc(j),delta_B,smooth*delta_B,H1)
                  call smooth_heaviside(Yc(j),L2-delta_B,smooth*delta_B,H2)

                  do k = xstart(3), n3e
                       do i = 1, n1
                            if(streamwise==3) then

                                B1x(i,j,k) = e1*sin(2.d0*pi*Xc(k)*num_period/L3)*(1-H1+H2)

                            elseif(streamwise==1) then

                                B1x(i,j,k) = e1*sin(2.d0*pi*Xc(k)*num_period*(dble(n3m)/dble(n1m))/L3)*(1-H1+H2)

                            endif
                        enddo
                  enddo
              enddo

              n2e = min(zend(2), n2-1)
              n1e = min(zend(1), n1-1)

              do j = zstart(2), n2e
                  call smooth_heaviside(Yc(j),delta_B,smooth*delta_B,H1)
                  call smooth_heaviside(Yc(j),L2-delta_B,smooth*delta_B,H2)

                  do k = 1, n3
                       do i = zstart(1), n1e
                            if(streamwise==3) then

                                B3z(i,j,k) = e3*sin(2.d0*pi*Zc(i)*num_period*(dble(n1m)/dble(n3m))/L1)*(1-H1+H2)

                            elseif(streamwise==1) then

                                B3z(i,j,k) = e3*sin(2.d0*pi*Zc(i)*num_period/L1)*(1-H1+H2)

                            endif
                        enddo
                  enddo
              enddo

              n3e = min(yend(3), n3-1)
              n1e = min(yend(1), n1-1)

              do j = 1, n2
                  call smooth_heaviside(Y(j),delta_B,smooth*delta_B,H1)
                  call smooth_heaviside(Y(j),L2-delta_B,smooth*delta_B,H2)

                  do k = ystart(3), n3e
                       do i = ystart(1), n1e
                            if(streamwise==3) then

                                B2y(i,j,k) = e2*sin(2.d0*pi*Xc(k)*num_period/L3)*(1-H1+H2)

                            elseif(streamwise==1) then

                                B2y(i,j,k) = e2*sin(2.d0*pi*Zc(i)*num_period/L1)*(1-H1+H2)

                            endif
                        enddo
                  enddo
              enddo

                call D1s_MULT_3Dy(B2y, dB2cy, ysize(1),ysize(1),n2,ysize(3),ysize(3), dx2, .false., antisymetric, Yc_to_YcTr_for_D1)

                if(streamwise==3) then

                    call transpose_y_to_z(dB2cy,dB3cz)
                    dB3cz = - dB3cz

                    n1e = min(zend(1), n1-1)
                    n2e = min(zend(2), n2-1)

                    do j = zstart(2), n2e
                       do i = zstart(1), n1e
                            do k = 2, n3
                                B3z(i,j,k) = B3z(i,j,k-1) + dB3cz(i,j,k-1)*dx3
                            enddo
                       enddo
                    enddo

                elseif(streamwise==1) then

                    call transpose_y_to_x(dB2cy,dB1cx)
                    dB1cx = - dB1cx

                    n3e = min(xend(3), n3-1)
                    n2e = min(xend(2), n2-1)

                    do j = xstart(2), n2e
                       do k = xstart(3), n3e
                            do i = 2, n1
                                B1x(i,j,k) = B1x(i-1,j,k) + dB1cx(i-1,j,k)*dx1
                            enddo
                       enddo
                    enddo

                endif

        end subroutine


    end subroutine generate_total_magnetic_field

!
end module MHD_field_generator
