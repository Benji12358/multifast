module Bubble_generator

    use decomp_2d
    use mesh
    use DNS_settings
    use physical_fields
    use subdomains_view
    use bubble_data
    use bubble_parallel_data
    use interpol
!    use VELOCITY_solver    ATTENTION BUBBLECHANGE

    implicit none

    real*8,dimension (:,:,:), allocatable       :: q1c_x, q2c_x, q3c_x
    real*8, dimension(:,:), allocatable         :: q1b20_x, q2b20_x, q3b20_x
    real*8, dimension(:,:), allocatable         :: q1b21_x, q2b21_x, q3b21_x

    real*8, dimension(:,:), allocatable         :: om1b20_x, om2b20_x, om3b20_x
    real*8, dimension(:,:), allocatable         :: om1b21_x, om2b21_x, om3b21_x

    real*8,dimension (:,:,:), allocatable       :: q1h2, q2h2, q3h2
    real*8,dimension (:,:,:), allocatable       :: om1h2, om2h2, om3h2
    real*8,dimension (:,:,:), allocatable       :: gradtau3Dh2



contains

    subroutine init_bubbles()
        implicit none

        allocate(q3c_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        q3c_x=0.d0
        allocate(q2c_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        q2c_x=0.d0
        allocate(q1c_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        q1c_x=0.d0
        allocate(om1_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        om1_x=0.d0
        allocate(om2_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        om2_x=0.d0
        allocate(om3_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        om3_x=0.d0

        allocate(gradtau3D_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        gradtau3D_x=0.d0


        !
        allocate(q1b20_x(xsize(1), xsize(3)))
        q1b20_x=0.d0
        allocate(q2b20_x(xsize(1), xsize(3)))
        q2b20_x=0.d0
        allocate(q3b20_x(xsize(1), xsize(3)))
        q3b20_x=0.d0
        allocate(q1b21_x(xsize(1), xsize(3)))
        q1b21_x=0.d0
        allocate(q2b21_x(xsize(1), xsize(3)))
        q2b21_x=0.d0
        allocate(q3b21_x(xsize(1), xsize(3)))
        q3b21_x=0.d0

        allocate(om1b20_x(xsize(1), xsize(3)))
        om1b20_x=0.d0
        allocate(om2b20_x(xsize(1), xsize(3)))
        om2b20_x=0.d0
        allocate(om3b20_x(xsize(1), xsize(3)))
        om3b20_x=0.d0

        allocate(om1b21_x(xsize(1), xsize(3)))
        om1b21_x=0.d0
        allocate(om2b21_x(xsize(1), xsize(3)))
        om2b21_x=0.d0
        allocate(om3b21_x(xsize(1), xsize(3)))
        om3b21_x=0.d0


        allocate(q1h2(0:n1, 0:xsize(2)+1, 0:xsize(3)+1))
        q1h2=0.d0
        allocate(q2h2(0:n1, 0:xsize(2)+1, 0:xsize(3)+1))
        q2h2=0.d0
        allocate(q3h2(0:n1, 0:xsize(2)+1, 0:xsize(3)+1))
        q3h2=0.d0
        allocate(om1h2(0:n1, 0:xsize(2)+1, 0:xsize(3)+1))
        om1h2=0.d0
        allocate(om2h2(0:n1, 0:xsize(2)+1, 0:xsize(3)+1))
        om2h2=0.d0
        allocate(om3h2(0:n1, 0:xsize(2)+1, 0:xsize(3)+1))
        om3h2=0.d0

        allocate(gradtau3Dh2(0:n1, 0:xsize(2)+1, 0:xsize(3)+1))
        gradtau3Dh2=0.d0

        allocate(fb1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        fb1=0.d0
        allocate(fb2(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        fb2=0.d0
        allocate(fb3(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        fb3=0.d0
        allocate(void_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        void_x=0.d0

        allocate(fb1_alpha_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        fb1_alpha_x=0.d0
        allocate(fb2_alpha_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        fb2_alpha_x=0.d0
        allocate(fb3_alpha_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        fb3_alpha_x=0.d0
        !allocate(fb1_alpha_y(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        !allocate(fb2_alpha_y(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        !allocate(fb1_alpha_z(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        !allocate(fb2_alpha_z(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))


    end subroutine init_bubbles

    ! BUBBLES GENERATION---------------------------------------------------------
    subroutine generate_bubbles(decomp, regen, growth_factor)
        implicit none
        integer, intent(in) :: decomp, regen, growth_factor
        integer             :: i, j, k, bg, nb_bubbles
        real*8              :: x1,x2,x3, v1,v2,v3, u01,u02,u03, u11,u12,u13, om01,om02,om03, om11,om12,om13, a1,a2,a3, b1,b2,b3, c1,c2,c3, d
        integer             :: proc_num
        type(bubble)        :: bb
        real*8              :: x1s, x2s, x3s, x1e, x2e, x3e

        if (nrank==0) then


            do bg = 1, bublgen_nb

                x1s=Z(max(1, bubl_gen(bg)%n1))
                x2s=Y(max(1, bubl_gen(bg)%n2))
                x3s=X(max(1, bubl_gen(bg)%n3))

                x1e=Z(min(n1, bubl_gen(bg)%n1+bubl_gen(bg)%s1))
                x2e=Y(min(n2, bubl_gen(bg)%n2+bubl_gen(bg)%s2))
                x3e=X(min(n3, bubl_gen(bg)%n3+bubl_gen(bg)%s3))

                if(regen==0) nb_bubbles=bubl_gen(bg)%nb_start
                if(regen==1)  nb_bubbles=bubl_gen(bg)%nb_regen*growth_factor


                do i = 1, nb_bubbles
                    call random_number(x1)
                    call random_number(x2)
                    call random_number(x3)

                    !                    call random_number(v1)
                    !                    call random_number(v2)
                    !                    call random_number(v3)

                    !x1=x1s+x1*(x1e-x1s)
                    x1=x1s
                    !x2=x2s+x2*(x2e-x2s)
                    x2=x2s
                    !x3=x3s+x3*(x3e-x3s)
                    x3=0.5d0*(x3s+x3e)

                    !x1=0.001d0
                    !x2=0.05d0*x2
                    !x3=(x3e-x3s)/2

                    !if(regen==0) x1 = 0.001d0

                    v1=0.366d0!1000.d0!-0.5d0+v1
                    v2=0.d0
                    v3=0.d0!1000.d0!-0.5d0+v3

                    u01=0.d0
                    u02=0.d0
                    u03=0.d0

                    u11=0.d0
                    u12=0.d0
                    u13=0.d0

                    om01=0.d0
                    om02=0.d0
                    om03=0.d0

                    om11=0.d0
                    om12=0.d0
                    om13=0.d0

                    call get_processor_number((/x1,x2,x3/), 1, proc_num)

                    bb%x1=x1
                    bb%x2=x2
                    bb%x3=x3

                    bb%v1=v1
                    bb%v2=v2
                    bb%v3=v3

                    bb%u01=0.d0
                    bb%u02=0.d0
                    bb%u03=0.d0

                    bb%u11=0.d0
                    bb%u12=0.d0
                    bb%u13=0.d0

                    bb%om01=0.d0
                    bb%om02=0.d0
                    bb%om03=0.d0

                    bb%om11=0.d0
                    bb%om12=0.d0
                    bb%om13=0.d0

                    bb%a1=0.d0
                    bb%a2=0.d0
                    bb%a3=0.d0

                    bb%b1=0.d0
                    bb%b2=0.d0
                    bb%b3=0.d0

                    bb%c1=0.d0
                    bb%c2=0.d0
                    bb%c3=0.d0

                    bb%d=0.d0

                    bb%proc=proc_num*1.d0

                    !                    if(regen==0) call export_bubble_data(bb)
                    call post_bubble(proc_num, bb)

                end do

            end do

        endif

        call send_bubble(bubl, bubble_cpt)

    end subroutine generate_bubbles
    subroutine generate_bubbles2(decomp, regen)
        implicit none
        integer, intent(in) :: decomp, regen
        integer             :: i, j, k, bg, nb_bubbles
        real*8              :: x1,x2,x3, v1,v2,v3, u01,u02,u03, u11,u12,u13, om01,om02,om03, om11,om12,om13, a1,a2,a3, b1,b2,b3, c1,c2,c3, d
        integer             :: proc_num
        type(bubble)        :: bb
        real*8              :: x1s, x2s, x3s, x1e, x2e, x3e
        integer, parameter  :: nby=100, nbz=10

        if (nrank==0) then


            do bg = 1, bublgen_nb

                x1s=Z(max(1, bubl_gen(bg)%n1))
                x2s=Y(max(1, bubl_gen(bg)%n2))
                x3s=X(max(1, bubl_gen(bg)%n3))

                x1e=Z(min(n1, bubl_gen(bg)%n1+bubl_gen(bg)%s1))
                x2e=Y(min(n2, bubl_gen(bg)%n2+bubl_gen(bg)%s2))
                x3e=X(min(n3, bubl_gen(bg)%n3+bubl_gen(bg)%s3))

                if(regen==0) nb_bubbles=nby*nbz
                if(regen==1)  nb_bubbles=bubl_gen(bg)%nb_regen


                do k = 1, nbz
                    do j = 1, nby
                        call random_number(x1)
                        call random_number(x2)
                        call random_number(x3)

                        !                    call random_number(v1)
                        !                    call random_number(v2)
                        !                    call random_number(v3)

                        x1=x1s!+x1*(x1e-x1s)
                        x2=x2s+(j*1.d0/nby)*(x2e-x2s)
                        x3=x3s+(k*1.d0/nbz)*(x3e-x3s)

                        v1=0.d0!1000.d0!-0.5d0+v1
                        v2=0.d0!1000.d0!-0.5d0+v2
                        v3=0.d0!1000.d0!-0.5d0+v3

                        u01=0.d0
                        u02=0.d0
                        u03=0.d0

                        u11=0.d0
                        u12=0.d0
                        u13=0.d0

                        om01=0.d0
                        om02=0.d0
                        om03=0.d0

                        om11=0.d0
                        om12=0.d0
                        om13=0.d0

                        call get_processor_number((/x1,x2,x3/), 1, proc_num)

                        bb%x1=x1
                        bb%x2=x2
                        bb%x3=x3

                        bb%v1=v1
                        bb%v2=v2
                        bb%v3=v3

                        bb%u01=0.d0
                        bb%u02=0.d0
                        bb%u03=0.d0

                        bb%u11=0.d0
                        bb%u12=0.d0
                        bb%u13=0.d0

                        bb%om01=0.d0
                        bb%om02=0.d0
                        bb%om03=0.d0

                        bb%om11=0.d0
                        bb%om12=0.d0
                        bb%om13=0.d0

                        bb%a1=0.d0
                        bb%a2=0.d0
                        bb%a3=0.d0

                        bb%b1=0.d0
                        bb%b2=0.d0
                        bb%b3=0.d0

                        bb%c1=0.d0
                        bb%c2=0.d0
                        bb%c3=0.d0

                        bb%d=0.d0

                        bb%proc=proc_num*1.d0

                        !                    if(regen==0) call export_bubble_data(bb)
                        call post_bubble(proc_num, bb)

                    end do
                end do

            end do

        endif

        call send_bubble(bubl, bubble_cpt)

    end subroutine generate_bubbles2
    subroutine generate_bubbles3(decomp, regen)
        implicit none
        integer, intent(in) :: decomp, regen
        integer             :: i, j, k, bg, nb_bubbles
        real*8              :: x1,x2,x3, v1,v2,v3, u01,u02,u03, u11,u12,u13, om01,om02,om03, om11,om12,om13, a1,a2,a3, b1,b2,b3, c1,c2,c3, d
        integer             :: proc_num
        type(bubble)        :: bb
        real*8              :: x1s, x2s, x3s, x1e, x2e, x3e

        integer, parameter  :: nby=100, nbz=10

        if (nrank==0) then


            do bg = 1, bublgen_nb

                x1s=Z(max(1, bubl_gen(bg)%n1))
                x2s=Y(max(1, bubl_gen(bg)%n2))
                x3s=X(max(1, bubl_gen(bg)%n3))

                x1e=Z(min(n1, bubl_gen(bg)%n1+bubl_gen(bg)%s1))
                x2e=Y(min(n2, bubl_gen(bg)%n2+bubl_gen(bg)%s2))
                x3e=X(min(n3, bubl_gen(bg)%n3+bubl_gen(bg)%s3))

                if(regen==0) nb_bubbles=nbz*nby
                if(regen==1)  nb_bubbles=bubl_gen(bg)%nb_regen


                do k = 1, nbz
                    do j = 1, nby

                        call random_number(x1)
                        call random_number(x2)
                        call random_number(x3)

                        !                    call random_number(v1)
                        !                    call random_number(v2)
                        !                    call random_number(v3)

                        x1=x1s
                        x2=x2s+(j*1.d0/nby)*(x2e-x2s)
                        x3=x3s+(k*1.d0/nbz)*(x3e-x3s)

                        v1=0.d0!1000.d0!-0.5d0+v1
                        v2=0.d0!1000.d0!-0.5d0+v2
                        v3=0.d0!1000.d0!-0.5d0+v3

                        u01=0.d0
                        u02=0.d0
                        u03=0.d0

                        u11=0.d0
                        u12=0.d0
                        u13=0.d0

                        om01=0.d0
                        om02=0.d0
                        om03=0.d0

                        om11=0.d0
                        om12=0.d0
                        om13=0.d0

                        call get_processor_number((/x1,x2,x3/), 1, proc_num)

                        bb%x1=x1
                        bb%x2=x2
                        bb%x3=x3

                        bb%v1=v1
                        bb%v2=v2
                        bb%v3=v3

                        bb%u01=0.d0
                        bb%u02=0.d0
                        bb%u03=0.d0

                        bb%u11=0.d0
                        bb%u12=0.d0
                        bb%u13=0.d0

                        bb%om01=0.d0
                        bb%om02=0.d0
                        bb%om03=0.d0

                        bb%om11=0.d0
                        bb%om12=0.d0
                        bb%om13=0.d0

                        bb%a1=0.d0
                        bb%a2=0.d0
                        bb%a3=0.d0

                        bb%b1=0.d0
                        bb%b2=0.d0
                        bb%b3=0.d0

                        bb%c1=0.d0
                        bb%c2=0.d0
                        bb%c3=0.d0

                        bb%d=0.d0

                        bb%proc=proc_num*1.d0

                        !                    if(regen==0) call export_bubble_data(bb)
                        call post_bubble(proc_num, bb)

                    end do

                end do

            end do

        endif

        call send_bubble(bubl, bubble_cpt)

    end subroutine generate_bubbles3
    subroutine generate_bubbles4(decomp, regen)
        implicit none
        integer, intent(in) :: decomp, regen
        integer             :: i, j, k, bg, nb_bubbles
        real*8              :: x1,x2,x3, v1,v2,v3, u01,u02,u03, u11,u12,u13, om01,om02,om03, om11,om12,om13, a1,a2,a3, b1,b2,b3, c1,c2,c3, d
        integer             :: proc_num
        type(bubble)        :: bb
        real*8              :: x1s, x2s, x3s, x1e, x2e, x3e

        integer, parameter  :: nbx=8000, nbz=50

        if (nrank==0) then


            do bg = 1, bublgen_nb

                x1s=Z(max(1, bubl_gen(bg)%n1))
                x2s=Y(max(1, bubl_gen(bg)%n2))
                x3s=X(max(1, bubl_gen(bg)%n3))

                x1e=Z(min(n1, bubl_gen(bg)%n1+bubl_gen(bg)%s1))
                x2e=Y(min(n2, bubl_gen(bg)%n2+bubl_gen(bg)%s2))
                x3e=X(min(n3, bubl_gen(bg)%n3+bubl_gen(bg)%s3))

                if(regen==0) nb_bubbles=nbz*nbx
                if(regen==1)  nb_bubbles=bubl_gen(bg)%nb_regen


                do k = 1, nbz
                    do i = 1, nbx

                        call random_number(x1)
                        call random_number(x2)
                        call random_number(x3)

                        !                    call random_number(v1)
                        !                    call random_number(v2)
                        !                    call random_number(v3)

                        x1=x1s+(i*1.d0/nbx)*(x1e-x1s)
                        x2=x2s
                        x3=x3s+(k*1.d0/nbz)*(x3e-x3s)

                        v1=0.d0!1000.d0!-0.5d0+v1
                        v2=0.d0!1000.d0!-0.5d0+v2
                        v3=0.d0!1000.d0!-0.5d0+v3

                        u01=0.d0
                        u02=0.d0
                        u03=0.d0

                        u11=0.d0
                        u12=0.d0
                        u13=0.d0

                        om01=0.d0
                        om02=0.d0
                        om03=0.d0

                        om11=0.d0
                        om12=0.d0
                        om13=0.d0

                        call get_processor_number((/x1,x2,x3/), 1, proc_num)

                        bb%x1=x1
                        bb%x2=x2
                        bb%x3=x3

                        bb%v1=v1
                        bb%v2=v2
                        bb%v3=v3

                        bb%u01=0.d0
                        bb%u02=0.d0
                        bb%u03=0.d0

                        bb%u11=0.d0
                        bb%u12=0.d0
                        bb%u13=0.d0

                        bb%om01=0.d0
                        bb%om02=0.d0
                        bb%om03=0.d0

                        bb%om11=0.d0
                        bb%om12=0.d0
                        bb%om13=0.d0

                        bb%a1=0.d0
                        bb%a2=0.d0
                        bb%a3=0.d0

                        bb%b1=0.d0
                        bb%b2=0.d0
                        bb%b3=0.d0

                        bb%c1=0.d0
                        bb%c2=0.d0
                        bb%c3=0.d0

                        bb%d=0.d0

                        bb%proc=proc_num*1.d0

                        !                    if(regen==0) call export_bubble_data(bb)
                        call post_bubble(proc_num, bb)

                    end do

                end do

            end do

        endif

        call send_bubble(bubl, bubble_cpt)

    end subroutine generate_bubbles4


    subroutine perform_fields_at_bubble(it, refresh_velocity)
        use decomp_2d
        use interpol
        use decomp2D_utils
        use subdomains_view
        use BUBBLE_fields_utils
        use physical_fields

        implicit none

        integer, intent(in)                                                         :: it, refresh_velocity

        real*8, dimension(3)                                                        :: pt

        real*8, dimension(0:n3+1)                                                   :: Xc_h
        real*8, dimension(0:n2+1)                                                   :: Yc_h
        real*8, dimension(0:n1)                                                     :: Zc_h


        integer                                                                     :: i, j, k, error
        logical                                                                     :: isout

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))   :: tau3D_y
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))   :: gradtau3D_y
        real*8, dimension(xsize(1), xsize(3))                                       :: gradtaub20_x, gradtaub21_x
        integer                                                                         ::  n1e, n2e, n3e

        Xc_h=2.d0   ! Verifying everything is setted correctly
        Yc_h=2.d0
        Zc_h=2.d0

        pt=0.d0

        Xc_h(1:n3-1)=Xc
        Xc_h(0)=-Xc(1)          ! PERIODIC
        Xc_h(n3)=Xc(n3-1)+dx3   ! PERIODIC

        Yc_h(1:n2-1)=Yc
        Yc_h(0)=Y(1)            ! DIRICHLET
        Yc_h(n2)=Y(n2)          ! DIRICHLET

        Zc_h(1:n1-1)=Zc
        Zc_h(0)=-Zc(1)          ! PERIODIC
        Zc_h(n1)=Zc(n1-1)+dx1   ! PERIODIC


        ! Calculer tau_3D <= du/dy
        ! warning : décomposition y nécessaire pour effectuer des dérivées en y
        n1e=ysize(1)!(min(n1m, yend(1))-ystart(1))+1
        n3e=ysize(3)!(min(n3m, yend(3))-ystart(3))+1
        call D1c_3Dy(q1_y, tau3D_y, ysize(1),n1e,n2,ysize(3),n3e, dx2, .true., NS_Q1_BC2)

        tau3D_y(:,1,:)=(q1_y(:,2,:)-q1_y(:,1,:))/(Yc(2)-Yc(1))
        tau3D_y(:,n2-1,:)=(q1_y(:,n2-1,:)-q1_y(:,n2-2,:))/(Yc(n2-1)-Yc(n2-2))


        do i = ystart(1), yend(1)
            do k = ystart(3), yend(3)
                gradtau3D_y(i,1:n2-1,k)=tau3D_y(i,1:n2-1,k)*Yc_to_YcTr_for_D2(:,1)
            end do
        end do

        call D2c_MULTACC_3Dy(q1_y, gradtau3D_y, ysize(1),n1e,n2,ysize(3),n3e, dx2, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,2))

        gradtau3D_y=gradtau3D_y/tau3D_y

        ! Changement de décomposition y=>x
        call transpose2Dy_y_to_x(gradtau3D_y(:,1,:), gradtaub20_x, xstart_glob(:,2))
        call transpose2Dy_y_to_x(gradtau3D_y(:,n2-1,:), gradtaub21_x, xstart_glob(:,2))
        call transpose_y_to_x(gradtau3D_y, gradtau3D_x)
        if (refresh_velocity==1) then

            call prepare_fields ! Preparing centered velocity fields and walls

            call treat_boundaries(q1c_x, q1b20_x, q1b21_x, q1h2) ! Working on 0, L1, L2 et L3 values
            call treat_boundaries(q2c_x, q2b20_x, q2b21_x, q2h2) ! Working on 0, L1, L2 et L3 values
            call treat_boundaries(q3c_x, q3b20_x, q3b21_x, q3h2) ! Working on 0, L1, L2 et L3 values

            call treat_boundaries(om1_x, om1b20_x, om1b21_x, om1h2)
            call treat_boundaries(om2_x, om2b20_x, om2b21_x, om2h2)
            call treat_boundaries(om3_x, om3b20_x, om3b21_x, om3h2)


            call treat_boundaries(om1_x, om1b20_x, om1b21_x, om1h2)

            call treat_boundaries(gradtau3D_x, gradtaub20_x, gradtaub21_x, gradtau3Dh2)

        endif

        ! DEBUG: pourquoi ne pas faire l'interpolation que pour refresh_velocity=1 ???
        do i = 1, bubble_cpt
            isout=.false.

            if ((bubl(i)%x1<=Z(1)).or.(bubl(i)%x1>=Z(n1))) isout=.true.
            if ((bubl(i)%x2<=Y(1)).or.(bubl(i)%x2>=Y(n2))) isout=.true.
            if ((bubl(i)%x3<=X(1)).or.(bubl(i)%x3>=X(n3))) isout=.true.

            if (.not. isout) then

                pt=(/bubl(i)%x1, bubl(i)%x2, bubl(i)%x3/)

                if (it==0) then

                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), q1h2, bubl(i)%u01, error)
                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), q2h2, bubl(i)%u02, error)
                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), q3h2, bubl(i)%u03, error)

                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), om1h2, bubl(i)%om01, error)
                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), om2h2, bubl(i)%om02, error)
                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), om3h2, bubl(i)%om03, error)

                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), gradtau3Dh2, bubl(i)%gradtau_tau0, error)

                !                    bubl(i)%u01=0.d0
                !                    bubl(i)%u02=0.d0
                !                    bubl(i)%u03=0.d0

                !                    bubl(i)%v1=bubl(i)%u01+1.d-5
                !                    bubl(i)%v2=bubl(i)%u02+1.d-5
                !                    bubl(i)%v3=bubl(i)%u03+1.d-5

                elseif (it==1) then

                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), q1h2, bubl(i)%u11, error)
                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), q2h2, bubl(i)%u12, error)
                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), q3h2, bubl(i)%u13, error)

                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), om1h2, bubl(i)%om11, error)
                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), om2h2, bubl(i)%om12, error)
                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), om3h2, bubl(i)%om13, error)

                !                    bubl(i)%u11=0.d0
                !                    bubl(i)%u12=0.d0
                !                    bubl(i)%u13=0.d0

                endif

            endif

        end do

    contains

        ! --------------------------------------------------------------------------------
        !                                Filling TABLES (INITIALIZATION)
        !---------------------------------------------------------------------------------

        subroutine prepare_fields()
            implicit none

            real*8, dimension(xsize(1), xsize(3))                                       :: q1_wall20_x
            real*8, dimension(xsize(1), xsize(3))                                       :: q1_wall21_x
            real*8, dimension(zstart(1):zend(1), zstart(3):zend(3))                     :: q3_wall20_z, q3_wall21_z
            real*8, dimension(zstart(1):zend(1), zstart(3):zend(3))                     :: q3b20_z, q3b21_z

            q1b20_x=0.d0
            q1b21_x=0.d0

            q3b20_z=0.d0
            q3b21_z=0.d0

            om1b20_x(:,:)=0.d0
            om2b20_x(:,:)=0.d0
            om3b20_x(:,:)=0.d0

            om1b21_x(:,:)=0.d0
            om2b21_x(:,:)=0.d0
            om3b21_x(:,:)=0.d0

            ! Transpose q1_wall2 in X1-decomp
            call transpose2Dy_y_to_x(q1_wall20, q1_wall20_x, xstart_glob(:,2))
            call transpose2Dy_y_to_x(q1_wall21, q1_wall21_x, xstart_glob(:,2))

            ! Prepare Wall - X1 direction
            do i = 1, n1-1
                q1b20_x(i,:)=0.5d0*(q1_wall20_x(i,:)+q1_wall20_x(i+1,:))
                q1b21_x(i,:)=0.5d0*(q1_wall21_x(i,:)+q1_wall21_x(i+1,:))
            end do

            ! Prepare Wall - X2 direction
            call transpose2Dy_y_to_x(q2_wall20, q2b20_x, xstart_glob(:,2))
            call transpose2Dy_y_to_x(q2_wall21, q2b21_x, xstart_glob(:,2))

            ! Transpose q3_wall2 in X3 decomp
            call transpose2Dy_y_to_z(q3_wall20, q3_wall20_z, zstart_glob(:,2))
            call transpose2Dy_y_to_z(q3_wall21, q3_wall21_z, zstart_glob(:,2))

            ! Prepare Wall - X3 direction

            do k = 1, n3-1
                q3b20_z(:,k)=0.5d0*(q3_wall20_z(:,k)+q3_wall20_z(:,k+1))
                q3b21_z(:,k)=0.5d0*(q3_wall21_z(:,k)+q3_wall21_z(:,k+1))
            end do

            call transpose2Dy_z_to_x(q3b20_z, q3b20_x, xstart_glob(:,2))
            call transpose2Dy_z_to_x(q3b21_z, q3b21_x, xstart_glob(:,2))

            call perform_velocity_at_center_x(q3_z, q2_y, q1_x, q3c_x, q2c_x, q1c_x)
            call perform_vorticity_x(q1c_x, q2c_x, q3c_x, om1_x, om2_x, om3_x)

            if(xstart(2)==1) then
                om1b20_x(:,:)=om1_x(:,1,:)
                om2b20_x(:,:)=om2_x(:,1,:)
                om3b20_x(:,:)=om3_x(:,1,:)
            endif

            if(xstart(2)==n2) then
                om1b21_x(:,:)=om1_x(:,n2-1,:)
                om2b21_x(:,:)=om2_x(:,n2-1,:)
                om3b21_x(:,:)=om3_x(:,n2-1,:)
            endif

        end subroutine prepare_fields



        ! --------------------------------------------------------------------------------
        !               Determinating the BOUNDARIES VALUES at 0, L1, L2 & L3
        !               -> Generates manual halo points to one direction X, Y or Z
        !                  to the intern workspace/area
        !---------------------------------------------------------------------------------


        subroutine treat_boundaries(qc_x, qb20_x,qb21_x, qc_x_h2)

            implicit none

            real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in)    :: qc_x
            real*8, dimension(xsize(1), xsize(3))                                                   :: qb20_x,qb21_x
            real*8, dimension(0:n1, 0:xsize(2)+1, 0:xsize(3)+1), intent(out)                        :: qc_x_h2

            real*8, dimension(:,:,:), allocatable                                                   :: qc_x_h
            real*8, dimension(:,:), allocatable                                                     :: qb20_x_h
            real*8, dimension(:,:), allocatable                                                     :: qb21_x_h

            integer                                                                                 :: end_proc
            integer                                                                                 :: start_proc
            integer                                                                                 :: mpi_err

            integer                                                                                 :: k1, k2, j1, j2

            qc_x_h2=100.d0

            ! halo values
            call update_halo(qc_x,qc_x_h,level=1) ! --> 0:n3+1, 0:n2+1, 1:n1

            ! Each centered halo value are copied to another table wich includes a wider range
            do i = 1, n1
                do j = 0, xsize(2)+1
                    do k = 0, xsize(3)+1
                        qc_x_h2(i, j, k)=qc_x_h(i, j, k) ! --> 0:n3+1, 0:n2+1, 0:n1
                    end do
                end do
            end do

            j1=0
            j2=xsize(2)+1
            if (xstart(2)==1) j1=1
            if (xend(2)==n2) j2=xsize(2)-1

            k1=0
            k2=xsize(3)+1
            if(xstart(3)==1) k1=1
            if(xend(3)==n3) k2=xsize(3)-1

            ! Set halo point for walls
            call halo_2Dx_2(qb20_x, qb20_x_h)
            call halo_2Dx_2(qb21_x, qb21_x_h)

            if (xstart(2)==1) then

                do i = 1, n1-1
                    do k = k1, k2
                        qc_x_h2(i, 0, k)=qb20_x_h(i,k)
                    end do
                end do

            endif

            ! --- Add of a ghost point at Y=L2
            if (xend(2)==n2) then

                do i = 1, n1-1
                    do k = k1, k2
                        qc_x_h2(i, xsize(2), k)=qb21_x_h(i,k)
                    end do
                end do

            endif

            ! Boundaries direction 3 ***************************************
            ! --- Add of a ghost point at X=0

            if(xstart(3)==1) then
                call get_processor_number(xstart(1),xstart(2),n3-1, 1, end_proc)
            endif


            if(xend(3)==n3) then
                call get_processor_number(xstart(1),xstart(2),1, 1, start_proc)
            endif

            if (xsize(3)/=n3) then

                if(xstart(3)==1) then
                    call MPI_SEND (qc_x_h2(1:n1-1,j1:j2,1),((j2-j1+1)*(n1-1)), MPI_DOUBLE_PRECISION, end_proc, 100, MPI_COMM_WORLD, mpi_err)
                    call MPI_RECV (qc_x_h2(1:n1-1,j1:j2,0),((j2-j1+1)*(n1-1)), MPI_DOUBLE_PRECISION ,end_proc,100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
                endif

                if(xend(3)==n3) then
                    call MPI_RECV (qc_x_h2(1:n1-1,j1:j2,xsize(3)),((j2-j1+1)*(n1-1)), MPI_DOUBLE_PRECISION ,start_proc,100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
                    call MPI_SEND (qc_x_h2(1:n1-1,j1:j2,xsize(3)-1),((j2-j1+1)*(n1-1)), MPI_DOUBLE_PRECISION, start_proc, 100, MPI_COMM_WORLD, mpi_err)
                endif

            else
                qc_x_h2(1:n1-1,j1:j2,0)=qc_x_h2(1:n1-1,j1:j2,xsize(3)-1)
                qc_x_h2(1:n1-1,j1:j2,xsize(3))=qc_x_h2(1:n1-1,j1:j2,1)
            end if

            ! Boundaries direction 1 ***************************************

            ! --- Add of a ghost point at Z=0

            do j = j1, j2

                do k = k1, k2
                    qc_x_h2(0, j, k)=qc_x_h2(n1-1, j, k)
                end do
            end do

            ! --- Add of a ghost point at Z=L1

            do j = j1, j2
                do k = k1, k2
                    qc_x_h2(n1, j, k)=qc_x_h2(1, j, k)
                end do
            end do

            ! --------------------------------------------------------------------------------
            !               Filling the boundaries conditions related to the edges
            !               -> Generates manual halo points on the edges of
            !                  the intern workspace/area, without the last
            !                  points at the end of the edges
            !---------------------------------------------------------------------------------

            ! Edge 3 - 2 ***************************************

            if (xsize(3)/=n3) then

                if((xstart(3)==1).and.(xstart(2)==1)) then
                    call MPI_SEND (qc_x_h2(:,0,1),(n1+1), MPI_DOUBLE_PRECISION, end_proc, 100, MPI_COMM_WORLD, mpi_err)
                    call MPI_RECV (qc_x_h2(:,0,0),(n1+1), MPI_DOUBLE_PRECISION ,end_proc,100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
                endif

                if((xend(3)==n3).and.(xstart(2)==1)) then
                    call MPI_RECV (qc_x_h2(:,0,xsize(3)),(n1+1), MPI_DOUBLE_PRECISION ,start_proc,100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
                    call MPI_SEND (qc_x_h2(:,0,xsize(3)-1),(n1+1), MPI_DOUBLE_PRECISION, start_proc, 100, MPI_COMM_WORLD, mpi_err)
                endif

                if ((xstart(3)==1).and.(xend(2)==n2)) then
                    call MPI_SEND (qc_x_h2(:,xsize(2),1),(n1+1), MPI_DOUBLE_PRECISION, end_proc, 100, MPI_COMM_WORLD, mpi_err)
                    call MPI_RECV (qc_x_h2(:,xsize(2),0),(n1+1), MPI_DOUBLE_PRECISION ,end_proc,100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
                endif

                if ((xend(3)==n3).and.(xend(2)==n2)) then
                    call MPI_RECV (qc_x_h2(:,xsize(2),xsize(3)),(n1+1), MPI_DOUBLE_PRECISION ,start_proc,100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
                    call MPI_SEND (qc_x_h2(:,xsize(2),xsize(3)-1),(n1+1), MPI_DOUBLE_PRECISION, start_proc, 100, MPI_COMM_WORLD, mpi_err)
                endif

            else

                qc_x_h2(:,0,0)=qc_x_h2(:,0,xsize(3)-1)
                qc_x_h2(:,0,xsize(3))=qc_x_h2(:,0,1)

                qc_x_h2(:,xsize(2),0)=qc_x_h2(:,xsize(2),xsize(3)-1)
                qc_x_h2(:,xsize(2),xsize(3))=qc_x_h2(:,xsize(2),1)

            end if

            ! Edge 3 - 1 ***************************************

            if (xstart(3)==1) then

                qc_x_h2(0, :, 0)=qc_x_h2(n1-1, :, 0)
                qc_x_h2(n1, :, 0)=qc_x_h2(1, :, 0)

            end if

            if (xend(3)==n3) then
                qc_x_h2(0, :, xsize(3))=qc_x_h2(n1-1, :, xsize(3))
                qc_x_h2(n1, :, xsize(3))=qc_x_h2(1, :, xsize(3))
            end if


            ! Edge 2 - 1 ***************************************

            if ((xstart(2)==1)) then

                do k = k1, k2
                    qc_x_h2(0, 0, :)=qc_x_h2(n1-1, 0, :)
                    qc_x_h2(n1, 0, :)=qc_x_h2(1, 0, :)
                end do

            endif

            if ((xend(2)==n2)) then

                do k = k1, k2
                    qc_x_h2(0, xsize(2), :)=qc_x_h2(n1-1, xsize(2), :)
                    qc_x_h2(n1, xsize(2), :)=qc_x_h2(1, xsize(2), :)
                end do

            endif

            ! --------------------------------------------------------------------------------
            !               Filling the last points boundaries conditions
            !               -> Generates manual halo points at the end of
            !                  the edges
            !---------------------------------------------------------------------------------


            ! Points on 1 - 3 ***************************************

            if ((xstart(1)==1).and.(xstart(2)==1).and.(xstart(3)==1)) then

                qc_x_h2(0, 0, 0)=qc_x_h2(n1-1, 0, 0)

            end if


            if ((xend(1)==n1).and.(xstart(2)==1).and.(xstart(3)==1)) then

                qc_x_h2(xsize(1), 0, 0)=qc_x_h2(1, 0, 0)

            end if

            if ((xstart(1)==1).and.(xend(2)==n2).and.(xstart(3)==1)) then
                qc_x_h2(0, xsize(2), 0)=qc_x_h2(n1-1, xsize(2), 0)
            end if


            if ((xend(1)==n1).and.(xend(2)==n2).and.(xstart(3)==1)) then
                qc_x_h2(xsize(1), xsize(2), 0)=qc_x_h2(1, xsize(2), 0)
            end if




            if ((xstart(1)==1).and.(xstart(2)==1).and.(xend(3)==n3)) then
                qc_x_h2(0, 0, xsize(3))=qc_x_h2(n1-1, 0, xsize(3))
            end if


            if ((xend(1)==n1).and.(xstart(2)==1).and.(xend(3)==n3)) then
                qc_x_h2(xsize(1), 0, xsize(3))=qc_x_h2(1, 0, xsize(3))
            end if

            if ((xstart(1)==1).and.(xend(2)==n2).and.(xend(3)==n3)) then
                qc_x_h2(0, xsize(2), xsize(3))=qc_x_h2(n1-1, xsize(2), xsize(3))
            end if


            if ((xend(1)==n1).and.(xend(2)==n2).and.(xend(3)==n3)) then
                qc_x_h2(xsize(1), xsize(2), xsize(3))=qc_x_h2(1, xsize(2), xsize(3))
            end if

        end subroutine treat_boundaries

    end subroutine perform_fields_at_bubble


    subroutine move_bubble(ntime,substep)

        use bubble_data
!        use VELOCITY_solver  !ATTENTION BUBBLECHANGE
        use time_schemes

        implicit none

        integer, intent(in)                 :: ntime
        integer                             :: ns, sub_cpt

        integer, intent(in)                 :: substep
        real*8, dimension(8)                :: t_tmp
        real*8                              :: dt2

        real*8, dimension(3)                :: u0, u1, v0, v1, om0, om1, alpha, alpha0
        real*8, dimension(3)                :: vec_u0, vec_u1, vec_v0, vec_v1
        real*8                              :: beta, gam, teta, arkhimedes, com1, com2, com4
        integer                             :: i, j, mpierr, k
        integer, save                       :: it=1
        real*8                              :: inter_ac

        real*8                              :: Re_particle, drag_coef
        real*8                              :: drag_force_1, drag_force_2, VM_force_1, VM_force_2, lift_force_1, lift_force_2, bouyant_force
        real*8                              :: alpha_bub, dalpha_dx2, f_alpha, beta_alpha, D_perp = 1.d0, u_stokes, gam_shear, v_dif_2, v_dif_2_old
        real*8                              :: x1_bub, x2_bub, xelec
        real*8                              :: Num_Bubbles
        real*8                              :: gradtau_tau0, gradtau_tau1, k_alpha

        u_stokes = -g_modified*(rad_modified**2)*ren/3 !Non-dimensionnal Stokes velocity

        !        open(16, file="BUBBLE.csv", position="append")

!        ATTENTION BUBBLECHANGE
!        do ns=1,nb_substep
!            ! First estimation of velocity with NO bubble effect
!            call update_velocity(ntime, ns)
!        enddo

        vec_u0=0.d0
        vec_u1=0.d0
        vec_v0=0.d0
        vec_v1=0.d0

        dt2=0.d0
        dt2=dt/substep

        com1=1.d0
        com2=1.d0
        com4=1.d0

        v_dif_2 = 0.d0
        v_dif_2_old = 0.d0

        open(1522, file="Number_of_Bubbles_Proc1", position="append")
        if(nrank==1) write(1522,*)'Number of bubbles in the Procesor 1:',bubble_cpt
        close(1522)

        do i = 1, bubble_cpt !We save the initial position of the bubbles before the iteration loop

            bubl_old(i)%x1 = bubl(i)%x1
            bubl_old(i)%x2 = bubl(i)%x2
            bubl_old(i)%x3 = bubl(i)%x3

        end do

        do k = 1,1 !We make n iterations in our coupled code

            do i = 1, bubble_cpt !We correct the initial positions of the bubble

                bubl(i)%x1 = bubl_old(i)%x1
                bubl(i)%x2 = bubl_old(i)%x2
                bubl(i)%x3 = bubl_old(i)%x3

            end do

            do sub_cpt = 1, sub_step !We move the bubble for all intermediate steps

                call perform_fields_at_bubble(1, sub_cpt) !We perform the new fields at the bubble

                do i = 1, bubble_cpt

                    !if ((bubl(i)%x1<=Z(1)).or.(bubl(i)%x1>=Z(n1))) then
                    !    bubl(i)%v1=-bubl(i)%v1
                    !end if

                    !if ((bubl(i)%x2<=Y(1)).or.(bubl(i)%x2>=Y(n2))) then
                    !    bubl(i)%v2=-bubl(i)%v2
                    !end if

                    !if ((bubl(i)%x3<=X(1)).or.(bubl(i)%x3>=X(n3))) then
                    !    bubl(i)%v3=-bubl(i)%v3
                    !end if

                    u0(1) = bubl(i)%u01
                    u1(1) = bubl(i)%u11
                    v0(1) = bubl(i)%v1

                    u0(2) = bubl(i)%u02
                    u1(2) = bubl(i)%u12
                    v0(2) = bubl(i)%v2

                    u0(3) = bubl(i)%u03
                    u1(3) = bubl(i)%u13
                    v0(3) = bubl(i)%v3

                    om0(1) = bubl(i)%om01
                    om1(1) = bubl(i)%om11

                    om0(2) = bubl(i)%om02
                    om1(2) = bubl(i)%om12

                    om0(3) = bubl(i)%om03
                    om1(3) = bubl(i)%om13

                    gradtau_tau0 = bubl(i)%gradtau_tau0
                    gradtau_tau1 = bubl(i)%gradtau_tau1

                    ! *****************Kernel of the solver*******************************

                    Re_particle = d*Uc*((u0(1)-v0(1))**2+(u0(2)-v0(2))**2+(u0(3)-v0(3))**2)**0.5d0/nu_viscosity
                    drag_coef = 16/Re_particle

                    beta = com1*(3.d0/8.d0)*(drag_coef/d_modified)*(divro)*dt2/(1.d0+com2*0.5d0*divro)
                    gam = com2*(3.d0/2.d0)*divro/(1.d0+0.5d0*divro)
                    teta = 0.5d0*Cl*divro*dt2/(1.d0+com2*0.5d0*divro)
                    arkhimedes = (1.d0-divro)*g_modified*dt2/(1.d0+com2*0.5d0*divro)

                    alpha0(1) = beta*abs(u0(1)-v0(1))*(u0(1)-v0(1))+v0(1)+gam*(u1(1)-u0(1))+com4*arkhimedes
                    alpha0(2) = beta*abs(u0(2)-v0(2))*(u0(2)-v0(2))+v0(2)+gam*(u1(2)-u0(2))
                    alpha0(3) = beta*abs(u0(3)-v0(3))*(u0(3)-v0(3))+v0(3)+gam*(u1(3)-u0(3))

                    call vect(u0, om0, vec_u0) !! Warning : Call vect was after the line
                    call vect(u1, om1, vec_u1) !! alpha0=alpha0+ teta*(vec_u0+vec_u1-vec_v0)
                    call vect(v0, om0, vec_v0)

                    alpha0=alpha0+ teta*(vec_u0+vec_u1-vec_v0) ! New position for this line

                    !! First estimation
                    v1=v0

                    !! CORE OF ITERATIVE METHOD
                    do j = 1, 10

                        call vect(v1, om1, vec_v1)

                        alpha = alpha0 - teta*vec_v1

                        call solve_CN(alpha(1), beta, u1(1), v1(1))
                        call solve_CN(alpha(2), beta, u1(2), v1(2))
                        call solve_CN(alpha(3), beta, u1(3), v1(3))

                    end do

                    ! ***************** End of the Kernel of the solver ********************+
                    ! **********Void Fraction mapping from COMSOL****************************
                    v_dif_2_old = v_dif_2
                    x1_bub = bubl(i)%x1/(0.5d0*L1)
                    x2_bub = bubl(i)%x2/(0.5d0*L1)
                    xelec=Z(max(1, bubl_gen(1)%n1))/(0.5d0*L1)
                    x1_bub=x1_bub-xelec
                    if (x1_bub>1.0d0)x1_bub=0.d0
                    if (x1_bub<0.0d0)x1_bub=0.d0
                    !*****Void mapping (given by Jonathan)***********************************
                    !alpha_bub = 0.043d0*(x1_bub**0.075d0)*EXP(-630.d0*(x1_bub**(-0.668d0))*(x2_bub**1.28d0))  ! Forced Conv. I=62, Re=27
                    !alpha_bub = 0.0288d0*(x1_bub**0.266d0)*EXP(-2266.d0*(x1_bub**(-0.434d0))*(x2_bub**1.43d0))  ! Forced Conv. I=62, Re=108
                    alpha_bub = 0.086d0*(x1_bub**0.127d0)*EXP(-836.d0*(x1_bub**(-0.583d0))*(x2_bub**1.29d0))  ! Forced Conv. I=125, Re=36
                    !alpha_bub = 0.0559d0*(x1_bub**0.277d0)*EXP(-2745.d0*(x1_bub**(-0.429d0))*(x2_bub**1.43d0)) ! Forced Conv. I=125, Re=144
                    !alpha_bub = 0.248d0*(x1_bub**0.381d0)*EXP(-1071.d0*(x1_bub**(-0.387d0))*(x2_bub**1.29d0)) ! Forced Conv. I=250, Re=50
                    !alpha_bub = 0.127d0*(x1_bub**0.344d0)*EXP(-3087.d0*(x1_bub**(-0.41d0))*(x2_bub**1.41d0)) ! Forced Conv. I=250, Re=200
                    !alpha_bub = 0.312d0*(x1_bub**0.439d0)*EXP(-14485*(x1_bub**(-0.555d0))*(x2_bub**1.73d0))  ! Forced Conv. I=500, Re=280
                    !*************************************************************************
                    !****** Derivative of alpha wuth respect to x2****************************
                    !dalpha_dx2 = -27.09d0*(x1_bub**(0.075-0.667))*EXP(-630.d0*(x1_bub**(-0.668d0))*(x2_bub**1.28d0))*(x2_bub**(0.28d0)/((0.5d0*L1)**1.28))
                    !dalpha_dx2 = -65.2608d0*(x1_bub**(0.266-0.434))*EXP(-2266.d0*(x1_bub**(-0.434d0))*(x2_bub**1.43d0))*(x2_bub**(0.43d0))/((0.5d0*L1)**1.43)
                    dalpha_dx2 = -71.896d0*(x1_bub**(0.127-0.583))*EXP(-836.d0*(x1_bub**(-0.583d0))*(x2_bub**1.29d0))*(x2_bub**(0.29d0))/((0.5d0*L1)**1.29)
                    !dalpha_dx2 = -153.4455d0*(x1_bub**(0.277-0.429))*EXP(-2745.d0*(x1_bub**(-0.429d0))*(x2_bub**1.43d0))*(x2_bub**(0.43d0))/((0.5d0*L1)**1.43)
                    !dalpha_dx2 = -265.608d0*(x1_bub**(0.381-0.387))*EXP(-1071.d0*(x1_bub**(-0.387d0))*(x2_bub**1.29d0))*(x2_bub**(0.29d0))/((0.5d0*L1)**1.29)
                    !dalpha_dx2 = -392.05d0*(x1_bub**(0.344-0.41))*EXP(-3087.d0*(x1_bub**(-0.41d0))*(x2_bub**1.41d0))*(x2_bub**(0.41d0))/((0.5d0*L1)**1.41)
                    !dalpha_dx2 = -4519.3d0*(x1_bub**(0.439-0.555))*EXP(-14485*(x1_bub**(-0.555d0))*(x2_bub**1.73d0))*(x2_bub**(0.73d0))/((0.5d0*L1)**1.73)
                    !*************************************************************************
                    gam_shear = abs(om0(3))
                    beta_alpha = (1/3)*(alpha_bub**2)*(1.d0+0.5d0*EXP(8.8d0*alpha_bub))
                    k_alpha = 0.6d0*alpha_bub**2 !*
                    f_alpha = (1-alpha_bub)**5
                    v_dif_2 = -((rad_modified**2)*gam_shear*beta_alpha+rad_modified*u_stokes*f_alpha*D_perp)*dalpha_dx2-(rad_modified**2)*gam_shear*k_alpha*gradtau_tau0    !/((1-alpha_bub)*alpha_bub)
                    !v1(2) = v_dif_2+v1(2)
                    ! ************************************************************************

                    bubl(i)%v1=v1(1)
                    bubl(i)%v2=v1(2)
                    bubl(i)%v3=v1(3)

                    ! Correct unexpected errors
                    !if (abs(bubl(i)%v1)>5.d0*Uc) then
                    !    bubl(i)%v1 = Uc+u_stokes
                    !end if
                    !if (abs(bubl(i)%v2)>Uc) then
                    !    bubl(i)%v2 = Uc/20.d0
                    !end if
                    !if (abs(bubl(i)%v3)>0.001d0) then
                    !    bubl(i)%v3 = 0.d0
                    !end if

                    ! We actualize the position to perform the action in the right place
                    bubl(i)%x1=bubl(i)%x1+(v0(1)+bubl(i)%v1)*dt2*0.5d0
                    bubl(i)%x2=bubl(i)%x2+(v0(2)+bubl(i)%v2)*dt2*0.5d0
                    bubl(i)%x3=bubl(i)%x3+(v0(3)+bubl(i)%v3)*dt2*0.5d0

                    ! We actualize the accelerration to perform the bubble action

                    !bubl(i)%a1 = bubl(i)%a1+(1/divro)*(v1(1)-v0(1))/dt-arkhimedes/substep
                    !************New way for performing acceleration of bubbles****************************************
                    bubl(i)%a1 = bubl(i)%a1+(3.d0/4.d0)*(drag_coef/d_modified)*(divro)*abs(bubl(i)%u11-bubl(i)%v1)*(bubl(i)%u11-bubl(i)%v1)/((1.d0-divro)*g_modified)/substep*(-dt)
                    bubl(i)%a2 = bubl(i)%a2+(3.d0/4.d0)*(drag_coef/d_modified)*(divro)*abs(bubl(i)%u12-bubl(i)%v2)*(bubl(i)%u12-bubl(i)%v2)/((1.d0-divro)*g_modified)/substep*(-dt)
                    bubl(i)%a3 = bubl(i)%a3+(1/divro)*(v1(3)-v0(3))/dt*(-dt)
                    !****************************************************************************************************
                    !*************Old way for performing the acceleration of bubbles*************************************
                    !bubl(i)%a1 = bubl(i)%a1+(1/divro)*(v1(1)-v0(1))/dt+arkhimedes/substep
                    !bubl(i)%a2 = bubl(i)%a2+(1/divro)*(v1(2)-v0(2))/dt
                    !bubl(i)%a3 = bubl(i)%a3+(1/divro)*(v1(3)-v0(3))/dt
                    !*****************************************************************************************************
                end do

            end do

            !call bubble_colllision

            !************Functions that includes the alpha map as an active force ******************
            !call include_active_alpha ! Comment to de-activate
            !***************************************************************************************

            call reassign_particle

            !if(mod(ntime, 10).eq.0) then
            call perform_bubble_action !We perform the bubble action for the new accleration and position => fb1,fb2,fb3

             !ATTENTION BUBBLECHANGE !!!!
!            call add_action(fb1, fb2, fb3) !We add the action to the body force in the NS solver
            !end if

!            ATTENTION BUBBLECHANGE !!!!
!            do ns=1,nb_substep
!
!                call update_multiphysics_velocity(ntime, ns) !We update the velocities for each node
!
!            end do

            !if(mod(ntime, 10).eq.0) then
            !    call substract_action(fb1, fb2, fb3) ! We clean the adde action once the velocity have been updated
            !end if

        end do

        open(1517, file="bubble_acceleration", position="append")
        open(1518, file="bubble_velocity", position="append")
        open(1519, file="bubble_position", position="append")
        open(1520, file="bubble_velocity_field", position="append")
        open(1521, file="bubble_vorticity_field", position="append")
        !do i = 1, bubble_cpt
        i=1
        if(nrank==1)    write(1517,*) bubl(i)%a1, bubl(i)%a2, bubl(i)%a3
        if(nrank==1)    write(1518,*) bubl(i)%v1, bubl(i)%v2, bubl(i)%v3
        if(nrank==1)    write(1519,*) x1_bub, x2_bub, bubl(1)%x3
        if(nrank==1)    write(1520,*) bubl(i)%u11, bubl(i)%u12, bubl(i)%u13
        if(nrank==1)    write(1521,*) bubl(i)%om11, bubl(i)%om12, bubl(i)%om13
        !end do
        close(1517)
        close(1518)
        close(1519)
        close(1520)
        close(1521)


        !**********Writting the relative forces in one bubble**************************************************************
        open(1516, file="Relative_forces", position="append")
        !do i = 1,bubble_cpt
        i=1
        Re_particle = d*Uc*((bubl(i)%u11-bubl(i)%v1)**2+(bubl(i)%u12-bubl(i)%v2)**2+(bubl(i)%u13-bubl(i)%v3)**2)**0.5d0/nu_viscosity
        drag_coef = 16/Re_particle
        bouyant_force = (1.d0-divro)*g
        drag_force_1 = (3.d0/4.d0)*(drag_coef/d)*(divro)*abs(bubl(i)%u11-bubl(i)%v1)*(bubl(i)%u11-bubl(i)%v1)*Uc**2/bouyant_force
        drag_force_2 = (3.d0/4.d0)*(drag_coef/d)*(divro)*abs(bubl(i)%u12-bubl(i)%v2)*(bubl(i)%u12-bubl(i)%v2)*Uc**2/bouyant_force
        VM_force_1 = (3.d0/2.d0)*divro*(bubl(i)%u11-bubl(i)%u01)*Uc**2/(h_height*dt2*bouyant_force)
        VM_force_2 = (3.d0/2.d0)*divro*(bubl(i)%u12-bubl(i)%u02)*Uc**2/(h_height*dt2*bouyant_force)
        lift_force_1 = Cl*divro*(vec_u1(1)-vec_v1(1))*Uc**2/(bouyant_force*h_height)
        lift_force_2 = Cl*divro*(vec_u1(2)-vec_v1(2))*Uc**2/(bouyant_force*h_height)
        bouyant_force = bouyant_force/bouyant_force
        if(nrank==1)   write(1516,*) drag_force_1, drag_force_2, VM_force_1, VM_force_2, lift_force_1, lift_force_2, bouyant_force
        !end do
        close(1516)
        !*****************************************************************************************************************

        do i = 1, bubble_cpt !We update the final fields at the bubble

            bubl(i)%v1 = v1(1)
            bubl(i)%v2 = v1(2)
            bubl(i)%v3 = v1(3)

            bubl(i)%u01 = u1(1)
            bubl(i)%u02 = u1(2)
            bubl(i)%u03 = u1(3)

            bubl(i)%om01 = om1(1)
            bubl(i)%om02 = om1(2)
            bubl(i)%om03 = om1(3)

            bubl(i)%d = g*(1.d0-divro)

            !*************Cleaning of the accelerations for the next time step*************************
            bubl(i)%a1 = 0.d0
            bubl(i)%a2 = 0.d0
            bubl(i)%a3 = 0.d0
            !******************************************************************************************
        end do


        call reassign_particle

    contains

        subroutine solve_CN(al, be, u, v)
            implicit none

            real*8, intent(in)  ::  be, al, u
            real*8, intent(out) ::  v

            real*8              ::  a, b, c, x1, x2

            integer, parameter  :: fid=125

            a=0.d0
            b=0.d0
            c=0.d0

            v=0.d0

            if ((al-u)<0.d0) then
                a = be
                b = - be*u*2.d0 - 1.d0
                c = al + be*u*u
                call solve_square(a, b, c, x1, x2)

                v = min(x1, x2)

            elseif ((al-u)>0.d0)then
                a = - be
                b = be*u*2.d0 - 1.d0
                c = al - be*u*u
                call solve_square(a, b, c, x1, x2)

                v = max(x1, x2)

            endif

        end subroutine solve_CN

        subroutine solve_square(a, b, c, x1, x2)
            implicit none

            real*8, intent(in)  ::  a, b, c
            real*8, intent(out) ::  x1, x2

            real*8              ::  delta

            x1=0.d0
            x2=0.d0
            delta=0.d0

            delta = b*b - 4.d0*a*c

            if (delta==0.d0) then
                x1 = (-b)/(2.d0*a)
                x2 = x1

            elseif (delta>0.d0) then
                x1 = (-b-dsqrt(delta))/(2.d0*a)
                x2 = (-b+dsqrt(delta))/(2.d0*a)
            endif

            if (a==0.d0) then
                x1 = -c/b
                x2 = x1
            endif

        end subroutine solve_square


        subroutine vect(A, B, AB)
            implicit none

            real*8, dimension(3), intent(in)    :: A, B
            real*8, dimension(3), intent(out)   :: AB

            AB=0.d0

            AB(1)=A(2)*B(3)-A(3)*B(2)
            AB(2)=A(3)*B(1)-A(1)*B(3)
            AB(3)=A(1)*B(2)-A(2)*B(1)

        end subroutine vect


    end subroutine move_bubble

    subroutine include_active_alpha

!        use VELOCITY_solver    !ATTENTION BUBBLECHANGE
        use bubble_data
        use boundaries
        use time_schemes

        implicit none

        real*8                              :: alpha_cell
        real*8                              :: x1_cell, x2_cell, xelec_ini, xelec_fin
        integer                             :: i, j, k
        integer                             :: n1s, n1e, n2s, n2e, n3s, n3e

        ! ************Determination of the procesor domain**************************************
        if (BC1==OPEN) then
            n1s=max(2, xstart(1))
            n1e=min(n1-2, xend(1))
        else
            n1s=xstart(1)
            n1e=min(n1-1, xend(1))
        endif

        n2s=xstart(2)
        n2e=min(n2-1, xend(2))

        n3s=xstart(3)
        n3e=min(n3-1, xend(3))

        xelec_ini=Z(max(1, bubl_gen(1)%n1))
        xelec_fin=Z(bubl_gen(1)%n1+INT(n1/2.d0))

        do k = n3s, n3e !Loop for x3
            do i = n1s, n1e !Loop for x1
                do j = n2s, n2e !Loop for x2
                    ! ***************************************************************************************
                    if (Zc(i) > xelec_ini .AND. Zc(i) <= xelec_fin) then
                        ! ***************** End of the Kernel of the solver ********************+
                        ! **********Void Fraction mapping from COMSOL****************************
                        x1_cell = (Zc(i)-xelec_ini)/(0.5d0*L1)
                        x2_cell = Y(j)/(0.5d0*L1)
                        !*****Void mapping (given by Jonathan)***********************************
                        !alpha_cell = 0.043d0*(x1_cell**0.075d0)*EXP(-630.d0*(x1_cell**(-0.668d0))*(x2_cell**1.28d0))
                        alpha_cell = 0.086d0*(x1_cell**0.127d0)*EXP(-836.d0*(x1_cell**(-0.583d0))*(x2_cell**1.29d0))
                        !alpha_cel = 0.248d0*(x1_cell**0.576d0)*EXP(-4768.d0*(x2_cell**1.54d0))
                        !**************Adition as an active velocity into de NS equation**********
                        fb1_alpha_x(i,j,k)=alpha_cell*g*dt
                        fb2_alpha_x(i,j,k)=0.d0
                        fb3_alpha_x(i,j,k)=0.d0
                        !**************************************************************************
                    endif
                end do
            end do
        end do

        !call transpose_y_to_x(fb1_alpha_y, fb1_alpha_x)
        !call transpose_y_to_z(fb1_alpha_y, fb1_alpha_z)
        !call transpose_y_to_x(fb2_alpha_y, fb2_alpha_x)
        !call transpose_y_to_z(fb2_alpha_y, fb2_alpha_z)

!       ATTENTION BUBBLECHANGE !!!!
!        call add_action(fb1_alpha_x, fb2_alpha_x, fb3_alpha_x)


    end subroutine include_active_alpha



    subroutine bubble_colllision

        use bubble_data

        implicit none

        integer                             :: i, j, num_collision, num_wall_collision
        real*8                              :: Ysup, Yinf, rad_bubble, distance, overlap, collision_distance, pos2

        num_collision = 0.d0
        num_wall_collision = 0.d0


        !Check bubble on bubble colllision
        rad_bubble = d/2.d0
        collision_distance = d*((1.d0+0.3d0)**(1.d0/3.d0))
        !open(1524, file="Collision", position="append")
        !if(nrank==1)   write(1524,*) 'Collision distance', collision_distance
        !if(nrank==1)   write(1524,*) '*******************************************************************'
        close(1524)
        do i = 1, bubble_cpt-1
            do j = i+1,bubble_cpt
                distance = sqrt((bubl(i)%x1-bubl(j)%x1)**2+(bubl(i)%x2-bubl(j)%x2)**2+(bubl(i)%x3-bubl(j)%x3)**2)
                if (distance < collision_distance) then
                    overlap = collision_distance-distance
                    num_collision = num_collision+1
                    call perform_collision(i,j,distance,overlap)
                    !open(1524, file="Collision", position="append")
                    !if(nrank==1)   write(1524,*) 'Final Distance', sqrt((bubl(i)%x1-bubl(j)%x1)**2+(bubl(i)%x2-bubl(j)%x2)**2+(bubl(i)%x3-bubl(j)%x3)**2)
                    !if(nrank==1)   write(1524,*) '*******************************************************************'
                    !close(1524)
                end if
            end do
        end do

        ! Check inferior wall collison
        Yinf = Y(1)
        do i = 1, bubble_cpt
            if(bubl(i)%x2<Yinf) then
                open(1525, file="Wall_Collison", position="append")
                if(nrank==1)   write(1525,*) 'Initial bubble position:', bubl(i)%x1, bubl(i)%x2, bubl(i)%x3
                if(nrank==1)   write(1525,*) 'Initial bubble velocity:', bubl(i)%v1, bubl(i)%v2, bubl(i)%v3
                call random_number(pos2)
                bubl(i)%x2 = pos2*0.01d0
                if (bubl(i)%v2>0) then
                    !bubl(i)%v2 = bubl(i)%v2*1.2d0
                    bubl(i)%v2 = abs(bubl(i)%v2)*1.2d0
                    bubl(i)%v3 = 0.d0
                else
                    !bubl(i)%v2 = -bubl(i)%v2*1.2d0
                    bubl(i)%v2 = abs(bubl(i)%v2)*1.2d0
                    bubl(i)%v3 = 0.d0
                end if
                if(nrank==1)   write(1525,*) 'Final bubble position:', bubl(i)%x1, bubl(i)%x2, bubl(i)%x3
                if(nrank==1)   write(1525,*) 'Final bubble velocity:', bubl(i)%v1, bubl(i)%v2, bubl(i)%v3
                if(nrank==1)   write(1525,*) Yinf
                if(nrank==1)   write(1525,*) '*********************************************************************'
                close(1525)
                num_wall_collision = num_wall_collision+1
            end if
        end do

        ! Check superior wall collision
        Ysup = Y(n2)
        do i = 1, bubble_cpt
            if(bubl(i)%x2+d/2>Ysup) then
                bubl(i)%x2 = Ysup
                bubl(i)%v2 = -bubl(i)%v2
                num_wall_collision = num_wall_collision+1
            end if
        end do


        open(1523, file="Number_of_collisions", position="append")
        if(nrank==1)   write(1523,*) 'Num of bub-bub colissions:', num_collision, 'Num of bub-wall colissions:', num_wall_collision, 'Num of bubbles:', bubble_cpt
        close(1523)

    end subroutine bubble_colllision

    subroutine perform_collision(i,j,distance,overlap)

        implicit none

        integer, intent(in)    :: i, j
        real*8, intent(in)     :: distance, overlap
        real*8, dimension(3)   :: dv, n, m, v_ns1, v_ns2
        real*8                 :: v_aux

        !open(1524, file="Collision", position="append")

        !if(nrank==1)   write(1524,*) 'Bubbles collision', i, '-', j
        !if(nrank==1)   write(1524,*) 'Initial velocity of bubble 1', bubl(i)%v1,bubl(i)%v2,bubl(i)%v3
        !if(nrank==1)   write(1524,*) 'Initial velocity of bubble 2', bubl(j)%v1,bubl(j)%v2,bubl(j)%v3
        !if(nrank==1)   write(1524,*) 'Initial position of bubble 1', bubl(i)%x1,bubl(i)%x2,bubl(i)%x3
        !if(nrank==1)   write(1524,*) 'Initial position of bubble 2', bubl(j)%x1,bubl(j)%x2,bubl(j)%x3
        !if(nrank==1)   write(1524,*) 'Distance:', distance, 'Overlap:', overlap

        ! Direction versor of the collision
        dv(1) = (bubl(i)%x1-bubl(j)%x1)/distance
        dv(2) = (bubl(i)%x2-bubl(j)%x2)/distance
        dv(3) = (bubl(i)%x3-bubl(j)%x3)/distance

        !if(nrank==1)   write(1524,*) 'Directoin of collision', dv(1),dv(2),dv(3)

        ! First perpendicular versor
        n(1) = 1/sqrt(1+(dv(1)/dv(3))**2)
        n(2) = 0
        n(3) = -dv(1)*n(1)/dv(3)

        !if(nrank==1)   write(1524,*) 'Normal 1', n(1),n(2),n(3)

        ! Second perpendicular vector
        m(1) = 1/sqrt(1+((dv(3)**2/dv(1)+dv(1))/dv(2))**2+(dv(3)/dv(1))**2)
        m(2) = -m(1)*(dv(3)**2/dv(1)+dv(1))/dv(2)
        m(3) = m(1)*dv(3)/dv(1)

        !if(nrank==1)   write(1524,*) 'Normal 2', m(1),m(2),m(3)

        ! Descomposition of the velocities of the bubbles that will collide
        v_ns1(1) = bubl(i)%v1*dv(1)+bubl(i)%v2*dv(2)+bubl(i)%v3*dv(3)
        v_ns1(2) = bubl(i)%v1*n(1)+bubl(i)%v2*n(2)+bubl(i)%v3*n(3)
        v_ns1(3) = bubl(i)%v1*m(1)+bubl(i)%v2*m(2)+bubl(i)%v3*m(3)

        !if(nrank==1)   write(1524,*) 'Velocity 1 in the  new system', v_ns1(1),v_ns1(2),v_ns1(3)

        v_ns2(1) = bubl(j)%v1*dv(1)+bubl(j)%v2*dv(2)+bubl(j)%v3*dv(3)
        v_ns2(2) = bubl(j)%v1*n(1)+bubl(j)%v2*n(2)+bubl(j)%v3*n(3)
        v_ns2(3) = bubl(j)%v1*m(1)+bubl(j)%v2*m(2)+bubl(j)%v3*m(3)

        !if(nrank==1)   write(1524,*) 'Velocity 2 in the  new system', v_ns2(1),v_ns2(2),v_ns2(3)

        ! Perform the collision (supposing equal mass on the bubbles)
        v_aux = v_ns1(1)
        v_ns1(1) = v_ns2(1)
        v_ns2(1) = v_aux

        !if(nrank==1)   write(1524,*) 'New velocities post collision', v_ns1(1),v_ns2(1)

        !Actualize the final velocities
        bubl(i)%v1 = v_ns1(1)*dv(1)+v_ns1(2)*n(1)+v_ns1(3)*m(1)
        bubl(i)%v2 = v_ns1(1)*dv(2)+v_ns1(2)*n(2)+v_ns1(3)*m(2)+abs(bubl(i)%v2)*0.2d0
        bubl(i)%v3 = v_ns1(1)*dv(3)+v_ns1(2)*n(3)+v_ns1(3)*m(3)

        !if(nrank==1)   write(1524,*) 'Final velocity of bubble 1', bubl(i)%v1,bubl(i)%v2,bubl(i)%v3

        bubl(j)%v1 = v_ns2(1)*dv(1)+v_ns2(2)*n(1)+v_ns2(3)*m(1)
        bubl(j)%v2 = v_ns2(1)*dv(2)+v_ns2(2)*n(2)+v_ns2(3)*m(2)+abs(bubl(j)%v2)*0.2d0
        bubl(j)%v3 = v_ns2(1)*dv(3)+v_ns2(2)*n(3)+v_ns2(3)*m(3)

        !if(nrank==1)   write(1524,*) 'Final velocity of bubble 2', bubl(j)%v1,bubl(j)%v2,bubl(j)%v3

        !Actualize the final positions
        bubl(i)%x1 = bubl(i)%x1+overlap*dv(1)/2.0
        bubl(i)%x2 = bubl(i)%x2+overlap*dv(2)/2.0
        bubl(i)%x3 = bubl(i)%x3+overlap*dv(3)/2.0

        !if(nrank==1)   write(1524,*) 'Final position of bubble 1:', bubl(i)%x1,bubl(i)%x2,bubl(i)%x3

        bubl(j)%x1 = bubl(j)%x1-overlap*dv(1)/2.0
        bubl(j)%x2 = bubl(j)%x2-overlap*dv(2)/2.0
        bubl(j)%x3 = bubl(j)%x3-overlap*dv(3)/2.0

        !if(nrank==1)   write(1524,*) 'Final position of bubble 2:', bubl(j)%x1,bubl(j)%x2,bubl(j)%x3

        !close(1524)

    end subroutine perform_collision


    subroutine reassign_particle() ! Checks the assignment of each bubble. A processor is able to send a bubble to another processor if it's not its property.
        implicit none

        real*8, dimension(3)    :: pt
        integer                 :: i, proc_num, mpi_err
        integer                 :: send, send_glob

        send=0

        do i = 1, bubble_cpt
            pt(1) = bubl(i)%x1
            pt(2) = bubl(i)%x2
            pt(3) = bubl(i)%x3
            call get_processor_number(pt, 1, proc_num)

            if ((proc_num/=nrank).and.(proc_num/=-1)) then
                send=1
                bubl(i)%proc = proc_num
                call post_bubble(proc_num, bubl(i))
                call remove_bubble(i)

            elseif (proc_num==-1) then
                call remove_bubble(i)

            ! Bubbles are eliminated before the end for stability
            elseif (pt(1) > (L1-L1/3.d0)) then
                call remove_bubble(i)
            end if

        end do

        call MPI_ALLREDUCE (send, send_glob, 1, MPI_INTEGER, MPI_SUM , MPI_COMM_WORLD , mpi_err)
        if (send_glob>0) call send_bubble(bubl, bubble_cpt)


        send=0
        bubble_tmp_cpt=0

        open(17, file="REASSIGN", position="append")

        do i = 1, bubble_cpt

            pt(1) = bubl(i)%x1
            pt(2) = bubl(i)%x2
            pt(3) = bubl(i)%x3

            if (xend(1)/=n1) then

                call get_processor_number((/Zc(xend(1)+1), pt(2), pt(3) /), 1, proc_num)

                if ((proc_num/=-1).and.(pt(1)>=Zc(xend(1))).and.(pt(1)<=Z(xend(1)+1))) then
                    send=1
                    call post_bubble(proc_num, bubl(i))
                endif

            end if

            if (xend(2)/=n2) then

                call get_processor_number((/pt(1), Yc(xend(2)+1), pt(3) /), 1, proc_num)
                !                write(17,*) 'X2: ', nrank, proc_num, pt(2), Yc(xend(2)), Y(xend(2)+1), Yc(xend(2)+1)!, xpos_end_glob(nrank, 2)

                if ((proc_num/=-1).and.(pt(2)>=Yc(xend(2))).and.(pt(2)<=Y(xend(2)+1))) then
                    write(17,*) 'X2: ', nrank, proc_num, pt(2), Yc(xend(2)), Y(xend(2)+1), Yc(xend(2)+1)!, xpos_end_glob(nrank, 2)
                    send=1
                    call post_bubble(proc_num, bubl(i))
                !                    write(17,*) send_list(proc_num, send_nb(proc_num), :)
                endif

            end if

            if (xend(3)/=n3) then

                call get_processor_number((/pt(1), pt(2), Xc(xend(3)+1) /), 1, proc_num)

                if ((proc_num/=-1).and.(pt(3)>=Xc(xend(3))).and.(pt(3)<=X(xend(3)+1))) then
                    write(17,*) 'X3: ', nrank, proc_num, pt(3), Xc(xend(3)), X(xend(3)+1), Xc(xend(3)+1)!, xpos_end_glob(nrank, 3)
                    send=1
                    call post_bubble(proc_num, bubl(i))
                !                    write(17,*) send_list(proc_num, send_nb(proc_num), :)
                endif

            end if

        end do

        close(17)

        !        write(*,*) '/////////////////////////////////////////// BUBL SL:', nrank, sum(abs(send_list)), bubble_tmp_cpt

        call MPI_ALLREDUCE (send, send_glob, 1, MPI_INTEGER, MPI_SUM , MPI_COMM_WORLD , mpi_err)
        if (send_glob>0) call send_bubble(bubl_tmp, bubble_tmp_cpt)

    !        write(*,*) '///////////////////////////////////////////buble_tmp_cpt:', nrank, bubble_tmp_cpt
    !        write(*,*) '////////////', xend(2), xend(3)

    end subroutine reassign_particle


    subroutine associate_volume_control(import, direction, volume, no_bound_err)
        implicit none

        type(bubble), intent(in)            :: import
        integer, intent(in)                 :: direction
        integer, dimension(3), intent(out)  :: volume
        logical, intent(out)                :: no_bound_err

        real*8, dimension(3)                :: pt
        integer                             :: proc_num
        real*8                              :: x0
        integer                             :: i, k
        integer                             :: i1, i2

        no_bound_err=.true.
        volume=0

        pt(1) = import%x1
        pt(2) = import%x2
        pt(3) = import%x3


        if (direction==1) then

            x0 = Zc(1)
            volume(1) = floor(((pt(1)-x0)/dx1))+2

            i2 = xend(2)
            if (xend(2)==n2) i2 = xend(2)-1
            do i = xstart(2), i2
                if ((pt(2)>=Y(i)).and.(pt(2)<=Y(i+1))) then
                    volume(2) = i
                end if
            end do

            volume(3) = floor(pt(3)/dx3) + 1

        end if

        ! ATTENTION, Conditions valables dans un fonctionnement avec parallélisme, mettre à jour en rajoutant une condition sur xstart ET xend
        if (direction==2) then

            volume(1) = floor(pt(1)/dx1) + 1

            i1=xstart(2)
            i2=xend(2)

            ! Volume(2): case bubble in the first (half) cell
            if (xstart(2)==1) then

                i1=2

                if ((pt(2)>=Y(1)).and.(pt(2)<=Yc(1))) then
                    volume(2) = 1
                endif
            endif

            ! Volume(2): case bubble in the last (half) cell
            if (xend(2)==n2) then

                i2=n2-1

                if ((pt(2)>=Yc(n2-1)).and.(pt(2)<=Y(n2))) then
                    volume(2) = n2
                endif
            endif

            ! Inner cases
            do i = i1, i2
                if ((pt(2)>=Yc(i-1)).and.(pt(2)<=Yc(i))) then
                    volume(2) = i
                end if
            end do

            volume(3) = floor(pt(3)/dx3) + 1

        end if


        if (direction==3) then

            volume(1) = floor(pt(1)/dx1) + 1

            i2 = xend(2)
            if (xend(2)==n2) i2 = xend(2)-1
            do i = xstart(2), i2
                if ((pt(2)>=Y(i)).and.(pt(2)<=Y(i+1))) then
                    volume(2) = i
                end if
            end do

            x0 = Xc(1)
            volume(3) = floor((pt(3)-x0)/dx3) + 2

        end if

        if (volume(1)>xend(1)) Then
            volume(1)=-1
            no_bound_err=.false.
        endif
        if (volume(2)>xend(2)) Then
            volume(2)=-1
            no_bound_err=.false.
        endif
        if (volume(3)>xend(3)) Then
            volume(3)=-1
            no_bound_err=.false.
        endif
        if (volume(1)<xstart(1)) Then
            volume(1)=-1
            no_bound_err=.false.
        endif
        if (volume(2)<xstart(2)) Then
            volume(2)=-1
            no_bound_err=.false.
        endif
        if (volume(3)<xstart(3)) Then
            volume(3)=-1
            no_bound_err=.false.
        endif

    end subroutine associate_volume_control

    subroutine perform_bubble_action()

        use bubble_data
        use mathematical_constants

        implicit none

        integer, dimension(3)   :: vc           ! Index of control volume
        integer                 :: i, j , k
        real*8                  :: vc_value, vbubble, alphad, global_arkhimedes, Conv_Fact       ! Volume of the indexed control volume
        logical                 :: no_bound_err

        global_arkhimedes = (1.d0-divro)*g*dt/(1.d0+0.5d0*divro)

        fb1=0.d0
        fb2=0.d0
        fb3=0.d0

        void_x=0.d0

        Conv_Fact = 0.d0

        vc_value= 0.d0
        alphad  = 0.d0
        vbubble = (1/6.d0)*pi*d*d*d


        do i = 1, bubble_cpt

            call associate_volume_control(bubl(i), 1, vc, no_bound_err)

            if (no_bound_err) then

                vc_value=dx1 * (Y(vc(2)+1)-Y(vc(2))) * dx3
                ! Limit cases
                if (vc(1)==1) vc_value=0.5d0*dx1* (Y(vc(2)+1)-Y(vc(2))) * dx3
                if (vc(1)==n1) vc_value=0.5d0*dx1* (Y(vc(2)+1)-Y(vc(2))) * dx3

                alphad = vbubble/vc_value
                fb1(vc(1), vc(2), vc(3))=fb1(vc(1), vc(2), vc(3)) - (bubl(i)%a1)*alphad*Conv_Fact
                void_x(vc(1), vc(2), vc(3)) = void_x(vc(1), vc(2), vc(3)) + alphad
            end if


            call associate_volume_control(bubl(i), 2, vc, no_bound_err)

            if (no_bound_err) then

                if (vc(2)==1)                   vc_value=dx1 * (Yc(1)-Y(1)) * dx3
                if (vc(2)==n2)                  vc_value=dx1 * (Y(n2)-Yc(n2-1)) * dx3
                if ((vc(2)/=1).and.(vc(2)/=n2)) vc_value=dx1 * (Yc(vc(2))-Yc(vc(2)-1)) * dx3

                alphad = vbubble/vc_value
                fb2(vc(1), vc(2), vc(3))=fb2(vc(1), vc(2), vc(3)) - (bubl(i)%a2)*alphad*Conv_Fact
            end if


            call associate_volume_control(bubl(i), 3, vc, no_bound_err)

            if (no_bound_err) then

                vc_value=dx1 * (Y(vc(2)+1)-Y(vc(2))) * dx3
                ! Limit cases
                if (vc(3)==1) vc_value=0.5d0*dx1* (Y(vc(2)+1)-Y(vc(2))) * dx3
                if (vc(3)==n3) vc_value=0.5d0*dx1* (Y(vc(2)+1)-Y(vc(2))) * dx3

                alphad = vbubble/vc_value
                fb3(vc(1), vc(2), vc(3))=fb3(vc(1), vc(2), vc(3)) !- (bubl(i)%a3)*alphad*Conv_Fact
            end if

            bubl(i)%a1 = 0.d0
            bubl(i)%a2 = 0.d0
            bubl(i)%a3 = 0.d0

        end do


        do i = 1, bubble_tmp_cpt

            call associate_volume_control(bubl_tmp(i), 1, vc, no_bound_err)

            if (no_bound_err) then

                vc_value=dx1 * (Y(vc(2)+1)-Y(vc(2))) * dx3
                ! Limit cases
                if (vc(1)==1) vc_value=0.5d0*dx1* (Y(vc(2)+1)-Y(vc(2))) * dx3
                if (vc(1)==n1) vc_value=0.5d0*dx1* (Y(vc(2)+1)-Y(vc(2))) * dx3

                alphad = vbubble/vc_value
                fb1(vc(1), vc(2), vc(3))=fb1(vc(1), vc(2), vc(3)) !- bubl_tmp(i)%a1*alphad*Conv_Fact
                void_x(vc(1), vc(2), vc(3)) = void_x(vc(1), vc(2), vc(3)) + alphad
            end if


            call associate_volume_control(bubl_tmp(i), 2, vc, no_bound_err)

            if (no_bound_err) then


                if (vc(2)==1)                   vc_value=dx1 * (Yc(1)-Y(1)) * dx3
                if (vc(2)==n2)                  vc_value=dx1 * (Y(n2)-Yc(n2-1)) * dx3
                if ((vc(2)/=1).and.(vc(2)/=n2)) vc_value=dx1 * (Yc(vc(2))-Yc(vc(2)-1)) * dx3

                alphad = vbubble/vc_value
                fb2(vc(1), vc(2), vc(3))=fb2(vc(1), vc(2), vc(3)) !- (bubl_tmp(i)%a2)*alphad*Conv_Fact
            end if


            call associate_volume_control(bubl_tmp(i), 3, vc, no_bound_err)

            if (no_bound_err) then

                vc_value=dx1 * (Y(vc(2)+1)-Y(vc(2))) * dx3
                ! Limit cases
                if (vc(3)==1) vc_value=0.5d0*dx1* (Y(vc(2)+1)-Y(vc(2))) * dx3
                if (vc(3)==n3) vc_value=0.5d0*dx1* (Y(vc(2)+1)-Y(vc(2))) * dx3

                alphad = vbubble/vc_value
                fb3(vc(1), vc(2), vc(3))=fb3(vc(1), vc(2), vc(3)) !- bubl_tmp(i)%a3*alphad*Conv_Fact
            end if


            !do k = 1, bubble_cpt ! Clean the acceleration vector
               !open(1768, file="bubble_acceleration", position="append")
               !write(1768,*) bubl(k)%a1, bubl(k)%a2, bubl(k)%a3
               !close(1768)
            !end do


            ! Clean the acceleration vector
            bubl(i)%a1 = 0.d0
            bubl(i)%a2 = 0.d0
            bubl(i)%a3 = 0.d0

        end do

        !open(1769, file="bubble_action_2", position="append")
        !if (nrank==1) write(1769,*) sum(abs(fb1)), sum(abs(fb2)), sum(abs(fb3)), no_bound_err
        !close(1769)

    end subroutine perform_bubble_action

    subroutine export_bubble_data(position, bub)

        implicit none

        type(bubble), intent(in)    :: bub
        integer, intent(in)         :: position

        open(23, file="bubble_data.csv", position="append")
        write(23,*) position, nrank, bub%x1, ',', bub%x2, ',', bub%x3, ',', bub%v1, ',', bub%v2, ',', bub%v3, ',', bub%u01, ',', bub%u02, ',', bub%u03, ',', bub%u11, ',', bub%u12, ',',  bub%u13, ',', bub%om01, ',', bub%om02, ',', bub%om03, ',', bub%om11, ',', bub%om12, ',', bub%om13, ',', bub%a1, ',', bub%a2, ',', bub%a3, ',', bub%d, ',', bub%proc
        close(23)

    end subroutine export_bubble_data

end module Bubble_generator

module Bubble_generator_test

    use decomp_2d
    use mesh
    use DNS_settings
    !    use physical_fields
    use subdomains_view
    use bubble_data
    use bubble_parallel_data
    use interpol

    implicit none

    real*8,dimension (:,:,:), allocatable       :: q1c_x, q2c_x, q3c_x
    real*8,dimension (:,:,:), allocatable       :: om1_x, om2_x, om3_x

    !    real*8,dimension (:,:), allocatable        :: q1b20_x, q2b20_x, q3b20_x
    !    real*8,dimension (:,:), allocatable        :: q1b21_x, q2b21_x, q3b21_x

    real*8, dimension(:,:), allocatable         :: q1b20_x, q2b20_x, q3b20_x
    real*8, dimension(:,:), allocatable         :: q1b21_x, q2b21_x, q3b21_x

    real*8, dimension(:,:), allocatable         :: om1b20_x, om2b20_x, om3b20_x
    real*8, dimension(:,:), allocatable         :: om1b21_x, om2b21_x, om3b21_x

    real*8,dimension (:,:,:), allocatable       :: q1h2, q2h2, q3h2
    real*8,dimension (:,:,:), allocatable       :: om1h2, om2h2, om3h2

contains

    subroutine init_bubbles()
        implicit none

        allocate(q3c_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        allocate(q2c_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        allocate(q1c_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        allocate(om1_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        allocate(om2_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        allocate(om3_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        !
        allocate(q1b20_x(xsize(1), xsize(3)))
        allocate(q2b20_x(xsize(1), xsize(3)))
        allocate(q3b20_x(xsize(1), xsize(3)))
        allocate(q1b21_x(xsize(1), xsize(3)))
        allocate(q2b21_x(xsize(1), xsize(3)))
        allocate(q3b21_x(xsize(1), xsize(3)))

        allocate(om1b20_x(xsize(1), xsize(3)))
        allocate(om2b20_x(xsize(1), xsize(3)))
        allocate(om3b20_x(xsize(1), xsize(3)))

        allocate(om1b21_x(xsize(1), xsize(3)))
        allocate(om2b21_x(xsize(1), xsize(3)))
        allocate(om3b21_x(xsize(1), xsize(3)))

        allocate(q1h2(0:n1, 0:xsize(2)+1, 0:xsize(3)+1))
        allocate(q2h2(0:n1, 0:xsize(2)+1, 0:xsize(3)+1))
        allocate(q3h2(0:n1, 0:xsize(2)+1, 0:xsize(3)+1))
        allocate(om1h2(0:n1, 0:xsize(2)+1, 0:xsize(3)+1))
        allocate(om2h2(0:n1, 0:xsize(2)+1, 0:xsize(3)+1))
        allocate(om3h2(0:n1, 0:xsize(2)+1, 0:xsize(3)+1))

        allocate(fb1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        allocate(fb2(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        allocate(fb3(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))

    end subroutine init_bubbles

    ! BUBBLES GENERATION---------------------------------------------------------
    subroutine generate_bubbles(decomp, regen)
        implicit none
        integer, intent(in) :: decomp, regen
        integer             :: i, j, k, bg, nb_bubbles
        real*8              :: x1,x2,x3, v1,v2,v3, u01,u02,u03, u11,u12,u13, om01,om02,om03, om11,om12,om13
        integer             :: proc_num
        type(bubble)        :: bb
        real*8              :: x1s, x2s, x3s, x1e, x2e, x3e

        if (nrank==0) then


            do bg = 1, bublgen_nb

                x1s=Z(max(1, bubl_gen(bg)%n1))
                x2s=Y(max(1, bubl_gen(bg)%n2))
                x3s=X(max(1, bubl_gen(bg)%n3))

                x1e=Z(min(n1, bubl_gen(bg)%n1+bubl_gen(bg)%s1))
                x2e=Y(min(n2, bubl_gen(bg)%n2+bubl_gen(bg)%s2))
                x3e=X(min(n3, bubl_gen(bg)%n3+bubl_gen(bg)%s3))

                if(regen==0) nb_bubbles=bubl_gen(bg)%nb_start
                if(regen==1)  nb_bubbles=bubl_gen(bg)%nb_regen


                do i = 1, nb_bubbles
                    call random_number(x1)
                    call random_number(x2)
                    call random_number(x3)

                    call random_number(v1)
                    call random_number(v2)
                    call random_number(v3)

                    x1=x1s+x1*(x1e-x1s)
                    x2=x2s+x2*(x2e-x2s)
                    x3=x3s+x3*(x3e-x3s)

                    v1=-0.5d0+v1
                    v2=-0.5d0+v2
                    !                    v2=-10.d0*v2
                    v3=-0.5d0+v3
                    !                    v3=10.d0*v3

                    u01=0.d0
                    u02=0.d0
                    u03=0.d0

                    u11=0.d0
                    u12=0.d0
                    u13=0.d0

                    om01=0.d0
                    om02=0.d0
                    om03=0.d0

                    om11=0.d0
                    om12=0.d0
                    om13=0.d0

                    call get_processor_number((/x1,x2,x3/), 1, proc_num)

                    bb%x1=x1
                    bb%x2=x2
                    bb%x3=x3

                    bb%v1=v1
                    bb%v2=v2
                    bb%v3=v3

                    bb%u01=0.d0
                    bb%u02=0.d0
                    bb%u03=0.d0

                    bb%u11=0.d0
                    bb%u12=0.d0
                    bb%u13=0.d0

                    bb%om01=0.d0
                    bb%om02=0.d0
                    bb%om03=0.d0

                    bb%om11=0.d0
                    bb%om12=0.d0
                    bb%om13=0.d0

                    bb%a1=0.d0
                    bb%a2=0.d0
                    bb%a3=0.d0

                    bb%b1=0.d0
                    bb%b2=0.d0
                    bb%b3=0.d0

                    bb%c1=0.d0
                    bb%c2=0.d0
                    bb%c3=0.d0

                    bb%d=0.d0

                    bb%proc=proc_num*1.d0

                    call post_bubble(proc_num, bb)

                end do

            end do

        endif

        call send_bubble(bubl, bubble_cpt)

    end subroutine generate_bubbles


    subroutine perform_fields_at_bubble(it, refresh_velocity)

        use interpol
        use decomp2D_utils
        use subdomains_view

        use BUBBLE_fields_utils
        !        use physical_fields

        implicit none

        real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))    :: q3_z                 !! Used to test move_bubble subroutine with init_fields subroutine
        real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))    :: q2_y                 !! Used to test move_bubble subroutine with init_fields subroutine
        real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))    :: q1_x                 !! Used to test move_bubble subroutine with init_fields subroutine
        real*8, dimension(ystart(1):yend(1), ystart(3):yend(3))                     :: q1_wall20, q2_wall20, q3_wall20 !! Used to test move_bubble subroutine with init_fields subroutine
        real*8, dimension(ystart(1):yend(1), ystart(3):yend(3))                     :: q1_wall21, q2_wall21, q3_wall21 !! Used to test move_bubble subroutine with init_fields subroutine

        integer, intent(in)                                                         :: it, refresh_velocity

        real*8, dimension(3)                                                        :: pt

        real*8, dimension(0:n3+1)                                                   :: Xc_h
        real*8, dimension(0:n2+1)                                                   :: Yc_h
        real*8, dimension(0:n1)                                                     :: Zc_h



        integer                                                                     :: i, j, k, error
        logical                                                                     :: isout

        real*8, dimension(3)                                                        :: value_th, value, delta

        character(200)                                                              :: file_path
        character*10                                                                :: file_nrank
        integer                                                                     :: funit

        integer                                                                     :: n1s, n2s, n3s
        integer                                                                     :: n1e, n2e, n3e

        funit=150+nrank
        write(file_nrank, "(i10)")nrank
        file_path='interpol_test_'//trim(adjustl(file_nrank))

        n1s=xstart(1)
        n2s=xstart(2)
        n3s=xstart(3)

        n1e=xend(1)
        n2e=xend(2)
        n3e=xend(3)

        Xc_h=2.d0   ! Verifying everything is setted correctly
        Yc_h=2.d0
        Zc_h=2.d0

        Xc_h(1:n3-1)=Xc
        Xc_h(0)=-Xc(1)          ! PERIODIC
        Xc_h(n3)=Xc(n3-1)+dx3   ! PERIODIC

        Yc_h(1:n2-1)=Yc
        Yc_h(0)=Y(1)            ! DIRICHLET
        Yc_h(n2)=Y(n2)          ! DIRICHLET

        Zc_h(1:n1-1)=Zc
        Zc_h(0)=-Zc(1)          ! PERIODIC
        Zc_h(n1)=Zc(n1-1)+dx1   ! PERIODIC

        open(funit, file=trim(file_path))
        write(funit, *) Z(n1s), Y(n2s), X(n3s)
        write(funit, *) Z(n1e), Y(n2e), X(n3e)

        if (refresh_velocity==1) then

            call init_fields
            call prepare_fields ! Preparing centered velocity fields and walls

            call treat_boundaries(q1c_x, q1b20_x, q1b21_x, q1h2) ! Working on 0, L1, L2 et L3 values
            call treat_boundaries(q2c_x, q2b20_x, q2b21_x, q2h2) ! Working on 0, L1, L2 et L3 values
            call treat_boundaries(q3c_x, q3b20_x, q3b21_x, q3h2) ! Working on 0, L1, L2 et L3 values

            write(funit, *) size(q1c_x,1), size(q1c_x,2), size(q1c_x,3)

            write(funit, *)
            write(funit, *) '0  : ', Zc_h(0), Yc_h(0), Xc_h(0)
            write(funit, *) '1  : ', Zc_h(1), Yc_h(1), Xc_h(1)
            write(funit, *) 'n-1: ', Zc_h(n1-1), Yc_h(n2-1), Xc_h(n3-1)
            write(funit, *) 'n  : ', Zc_h(n1), Yc_h(n2), Xc_h(n3)
            write(funit, *)

        endif

        do i = 1, bubble_cpt
            isout=.false.

            if ((bubl(i)%x1<=Z(1)).or.(bubl(i)%x1>=Z(n1))) isout=.true.
            if ((bubl(i)%x2<=Y(1)).or.(bubl(i)%x2>=Y(n2))) isout=.true.
            if ((bubl(i)%x3<=X(1)).or.(bubl(i)%x3>=X(n3))) isout=.true.

            if (.not. isout) then

                pt=(/bubl(i)%x1, bubl(i)%x2, bubl(i)%x3/)

                if (it==0) then

                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), q1h2, bubl(i)%u01, error)
                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), q2h2, bubl(i)%u02, error)
                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), q3h2, bubl(i)%u03, error)

                    value_th(1)=dcos(pt(1))*dcos(pt(2))*dcos(pt(3))
                    delta(1)=abs(bubl(i)%u01-value_th(1))/abs(value_th(1))

                    value_th(2)=dsin(pt(1))*dcos(pt(2))*dcos(pt(3))
                    delta(2)=abs(bubl(i)%u02-value_th(2))/abs(value_th(2))

                    value_th(3)=dcos(pt(1))*dsin(pt(2))*dcos(pt(3))
                    delta(3)=abs(bubl(i)%u03-value_th(3))/abs(value_th(3))

                    if ((pt(1)<Zc(1)).and.(pt(3)>Xc(n3-1))) then

                        write(funit, *) "ARRETE31 numero 3"
                        write(funit, *) pt, delta, error, sum(abs(q1h2(0,:,xsize(3))))
                        write(funit, *)
                        call interpolate_point2(pt, Zc_h(0:n1), Yc_h(n2s-1:n2e+1), Xc_h(n3s-1:n3e+1), q1h2, bubl(i)%u03, error, funit)
                        call interpolate_point2(pt, Zc_h(0:n1), Yc_h(n2s-1:n2e+1), Xc_h(n3s-1:n3e+1), q2h2, bubl(i)%u03, error, funit)
                        call interpolate_point2(pt, Zc_h(0:n1), Yc_h(n2s-1:n2e+1), Xc_h(n3s-1:n3e+1), q3h2, bubl(i)%u03, error, funit)
                    endif

                    if(maxval(abs(delta))>0.01d0) then
                        if (abs(pt(2)-1.57d0)>0.01d0) then

                            write(funit, *) "WARNING"
                            write(funit, *) pt, delta, error
                            write(funit, *)

                            if(nrank==0) call interpolate_point2(pt, Zc_h(0:n1), Yc_h(n2s-1:n2e+1), Xc_h(n3s-1:n3e+1), q3h2, bubl(i)%u03, error, funit)

                        endif
                    else
                        write(funit, *) pt, delta, error
                    endif

                elseif (it==1) then

                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), q1h2, bubl(i)%u11, error)
                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), q2h2, bubl(i)%u12, error)
                    call interpolate_point(pt, Zc_h(0:n1), Yc_h(xstart(2)-1:xend(2)+1), Xc_h(xstart(3)-1:xend(3)+1), q3h2, bubl(i)%u13, error)

                    value_th(1)=dcos(pt(1))*dcos(pt(2))*dcos(pt(3))
                    delta(1)=abs(bubl(i)%u11-value_th(1))/abs(value_th(1))

                    value_th(2)=dsin(pt(1))*dcos(pt(2))*dcos(pt(3))
                    delta(2)=abs(bubl(i)%u12-value_th(2))/abs(value_th(2))

                    value_th(3)=dcos(pt(1))*dsin(pt(2))*dcos(pt(3))
                    delta(3)=abs(bubl(i)%u13-value_th(3))/abs(value_th(3))

                    if ((pt(1)<Zc(1)).and.(pt(3)>Xc(n3-1))) then

                        write(funit, *) "ARRETE31 numero 3"
                        write(funit, *) pt, delta, error, sum(abs(q1h2(0,:,xsize(3))))
                        write(funit, *)
                        call interpolate_point2(pt, Zc_h(0:n1), Yc_h(n2s-1:n2e+1), Xc_h(n3s-1:n3e+1), q1h2, bubl(i)%u13, error, funit)
                        call interpolate_point2(pt, Zc_h(0:n1), Yc_h(n2s-1:n2e+1), Xc_h(n3s-1:n3e+1), q2h2, bubl(i)%u13, error, funit)
                        call interpolate_point2(pt, Zc_h(0:n1), Yc_h(n2s-1:n2e+1), Xc_h(n3s-1:n3e+1), q3h2, bubl(i)%u13, error, funit)
                    endif

                    if(maxval(abs(delta))>0.01d0) then
                        if (abs(pt(2)-1.57d0)>0.01d0) then

                            write(funit, *) "WARNING"
                            write(funit, *) pt, delta, error
                            write(funit, *)

                            if(nrank==0) call interpolate_point2(pt, Zc_h(0:n1), Yc_h(n2s-1:n2e+1), Xc_h(n3s-1:n3e+1), q3h2, bubl(i)%u13, error, funit)

                        endif
                    else
                        write(funit, *) pt, delta, error
                    endif

                endif

                if ((pt(1)<Zc(1)).and.(pt(3)<Xc(1))) then

                    write(funit, *) "ARRETE31 numero 1"
                    write(funit, *) pt, delta, error
                    write(funit, *)
                endif

                if ((pt(1)>Zc(n1-1)).and.(pt(3)<Xc(1))) then

                    write(funit, *) "ARRETE31 numero 2"
                    write(funit, *) pt, delta, error
                    write(funit, *)
                endif

                if ((pt(1)>Zc(n1-1)).and.(pt(3)>Xc(n3-1))) then

                    write(funit, *) "ARRETE31 numero 4"
                    write(funit, *) pt, delta, error
                    write(funit, *)
                endif

            endif

        end do

        close(funit)

    contains

        ! --------------------------------------------------------------------------------
        !                                Filling TABLES (INITIALIZATION)
        !---------------------------------------------------------------------------------
        subroutine init_fields()

            implicit none

            ! *****************************************************************
            ! Filling of boundaries conditions - Velocity   *******************
            ! *****************************************************************

            ! q1
            do i = 1, n1
                do j = xstart(2), min(xend(2), n2-1)
                    do k = xstart(3), min(xend(3), n3-1)
                        q1_x(i, j, k)=dcos(Z(i))*dcos(Yc_h(j))*dcos(Xc_h(k))
                    end do
                end do
            end do

            ! q2
            do i = ystart(1), min(yend(1), n1-1)
                do j = 1, n2
                    do k = ystart(3), min(yend(3), n3-1)
                        q2_y(i, j, k)=dsin(Zc_h(i))*dcos(Y(j))*dcos(Xc_h(k))
                    end do
                end do
            end do

            ! q3
            do i = zstart(1), min(zend(1), n1-1)
                do j = zstart(2), min(zend(2), n2-1)
                    do k = 1, n3
                        q3_z(i, j, k)=dcos(Zc_h(i))*dsin(Yc_h(j))*dcos(X(k))
                    end do
                end do
            end do

            ! *****************************************************************
            ! Filling of boundaries conditions - Walls      *******************
            ! *****************************************************************

            ! q2_wall ---------------------------------------------------------
            do i = ystart(1), min(yend(1),n1-1)
                do k = ystart(3), min(yend(3),n3-1)

                    q2_wall20(i,k)=dsin(Zc_h(i))*dcos(Yc_h(0))*dcos(Xc_h(k))
                    q2_wall21(i,k)=dsin(Zc_h(i))*dcos(Yc_h(n2))*dcos(Xc_h(k))

                end do
            end do


            ! q1_wall ---------------------------------------------------------
            q1_wall20=0.d0
            q1_wall21=0.d0
            do i = ystart(1), yend(1)
                do k = ystart(3), min(yend(3),n3-1)

                    q1_wall20(i,k)=dcos(Z(i))*dcos(Yc_h(0))*dcos(Xc_h(k))
                    q1_wall21(i,k)=dcos(Z(i))*dcos(Yc_h(n2))*dcos(Xc_h(k))

                end do
            end do

            ! q3_wall ---------------------------------------------------------
            q3_wall20=0.d0
            q3_wall21=0.d0
            do k = ystart(3), yend(3)
                do i = ystart(1), min(yend(1),n1-1)
                    q3_wall20(i,k)=dcos(Zc_h(i))*dsin(Yc_h(0))*dcos(X(k))
                    q3_wall21(i,k)=dcos(Zc_h(i))*dsin(Yc_h(n2))*dcos(X(k))
                enddo
            end do

        end subroutine init_fields

        subroutine prepare_fields()
            implicit none

            real*8, dimension(xsize(1), xsize(3))                                       :: q1_wall20_x
            real*8, dimension(xsize(1), xsize(3))                                       :: q1_wall21_x
            real*8, dimension(zstart(1):zend(1), zstart(3):zend(3))                     :: q3_wall20_z, q3_wall21_z
            real*8, dimension(zstart(1):zend(1), zstart(3):zend(3))                     :: q3b20_z, q3b21_z

            !                q1_wall20_x=0.d0
            !                q1_wall21_x=0.d0

            q1b20_x=0.d0
            q1b21_x=0.d0

            !                q2b20_x=0.d0
            !                q2b21_x=0.d0
            !
            !                q3_wall20_z=0.d0
            !                q3_wall21_z=0.d0

            q3b20_z=0.d0
            q3b21_z=0.d0

            !                q3b20_x=0.d0
            !                q3b21_x=0.d0
            !
            !                q3c_x=0.d0
            !                q2c_x=0.d0
            !                q1c_x=0.d0
            !
            !                om1_x=0.d0
            !                om2_x=0.d0
            !                om3_x=0.d0

            om1b20_x(:,:)=0.d0
            om2b20_x(:,:)=0.d0
            om3b20_x(:,:)=0.d0

            om1b21_x(:,:)=0.d0
            om2b21_x(:,:)=0.d0
            om3b21_x(:,:)=0.d0

            ! Transpose q1_wall2 in X1-decomp
            call transpose2Dy_y_to_x(q1_wall20, q1_wall20_x, xstart_glob(:,2))
            call transpose2Dy_y_to_x(q1_wall21, q1_wall21_x, xstart_glob(:,2))

            ! Prepare Wall - X1 direction
            do i = 1, n1-1
                q1b20_x(i,:)=0.5d0*(q1_wall20_x(i,:)+q1_wall20_x(i+1,:))
                q1b21_x(i,:)=0.5d0*(q1_wall21_x(i,:)+q1_wall21_x(i+1,:))
            end do

            ! Prepare Wall - X2 direction
            call transpose2Dy_y_to_x(q2_wall20, q2b20_x, xstart_glob(:,2))
            call transpose2Dy_y_to_x(q2_wall21, q2b21_x, xstart_glob(:,2))

            ! Transpose q3_wall2 in X3 decomp
            call transpose2Dy_y_to_z(q3_wall20, q3_wall20_z, zstart_glob(:,2))
            call transpose2Dy_y_to_z(q3_wall21, q3_wall21_z, zstart_glob(:,2))

            ! Prepare Wall - X3 direction
            q3b20_z=0.d0
            q3b21_z=0.d0

            do k = 1, n3-1
                q3b20_z(:,k)=0.5d0*(q3_wall20_z(:,k)+q3_wall20_z(:,k+1))
                q3b21_z(:,k)=0.5d0*(q3_wall21_z(:,k)+q3_wall21_z(:,k+1))
            end do

            call transpose2Dy_z_to_x(q3b20_z, q3b20_x, xstart_glob(:,2))
            call transpose2Dy_z_to_x(q3b21_z, q3b21_x, xstart_glob(:,2))

            call perform_velocity_at_center_x(q3_z, q2_y, q1_x, q3c_x, q2c_x, q1c_x)
            call perform_vorticity_x(q1c_x, q2c_x, q3c_x, om1_x, om2_x, om3_x)

            if(xstart(2)==1) then
                om1b20_x(:,:)=om1_x(:,1,:)
                om2b20_x(:,:)=om2_x(:,1,:)
                om3b20_x(:,:)=om3_x(:,1,:)
            endif

            if(xstart(2)==n2) then
                om1b21_x(:,:)=om1_x(:,n2-1,:)
                om2b21_x(:,:)=om2_x(:,n2-1,:)
                om3b21_x(:,:)=om3_x(:,n2-1,:)
            endif

        end subroutine prepare_fields



        ! --------------------------------------------------------------------------------
        !               Determinating the BOUNDARIES VALUES at 0, L1, L2 & L3
        !               -> Generates manual halo points to one direction X, Y or Z
        !                  to the intern workspace/area
        !---------------------------------------------------------------------------------


        subroutine treat_boundaries(qc_x, qb20_x,qb21_x, qc_x_h2)

            implicit none
            integer :: end_proc
            integer :: start_proc
            integer :: mpi_err
            real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))    :: qc_x                 !! Used to test move_bubble subroutine with init_fields subroutine
            real*8, dimension(xsize(1), xsize(3))                                       :: qb20_x,qb21_x

            real*8, dimension(:,:,:), allocatable                                       :: qc_x_h
            real*8, dimension(0:n1, 0:xsize(2)+1, 0:xsize(3)+1)                         :: qc_x_h2

            real*8, dimension(:,:), allocatable                                         :: qb20_x_h
            real*8, dimension(:,:), allocatable                                         :: qb21_x_h

            integer                                                                     :: k1, k2, j1, j2

            qc_x_h2=100.d0

            ! halo values
            call update_halo(qc_x,qc_x_h,level=1) ! --> 0:n3+1, 0:n2+1, 1:n1

            ! Each centered halo value are copied to another table wich includes a wider range
            do i = 1, n1
                do j = 0, xsize(2)+1
                    do k = 0, xsize(3)+1
                        qc_x_h2(i, j, k)=qc_x_h(i, j, k) ! --> 0:n3+1, 0:n2+1, 0:n1
                    end do
                end do
            end do

            j1=0
            j2=xsize(2)+1
            if (xstart(2)==1) j1=1
            if (xend(2)==n2) j2=xsize(2)-1

            k1=0
            k2=xsize(3)+1
            if(xstart(3)==1) k1=1
            if(xend(3)==n3) k2=xsize(3)-1

            ! Set halo point for walls
            call halo_2Dx_2(qb20_x, qb20_x_h)

            call halo_2Dx_2(qb21_x, qb21_x_h)

            if (xstart(2)==1) then

                do i = 1, n1-1
                    do k = k1, k2
                        qc_x_h2(i, 0, k)=qb20_x_h(i,k)
                    end do
                end do

            endif

            ! --- Add of a ghost point at Y=L2
            if (xend(2)==n2) then

                do i = 1, n1-1
                    do k = k1, k2
                        qc_x_h2(i, xsize(2), k)=qb21_x_h(i,k)
                    end do
                end do

            endif

            ! Boundaries direction 3 ***************************************
            ! --- Add of a ghost point at X=0

            if(xstart(3)==1) then
                call get_processor_number(xstart(1),xstart(2),n3-1, 1, end_proc)
            endif
            write(*,*)nrank, '////////endproc,', end_proc


            if(xend(3)==n3) then
                call get_processor_number(xstart(1),xstart(2),1, 1, start_proc)
            endif
            write(*,*)nrank, '////////start_proc,', start_proc

            if (xsize(3)/=n3) then

                if(xstart(3)==1) then
                    call MPI_SEND (qc_x_h2(1:n1-1,j1:j2,1),((j2-j1+1)*(n1-1)), MPI_DOUBLE_PRECISION, end_proc, 100, MPI_COMM_WORLD, mpi_err)
                    call MPI_RECV (qc_x_h2(1:n1-1,j1:j2,0),((j2-j1+1)*(n1-1)), MPI_DOUBLE_PRECISION ,end_proc,100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
                endif

                if(xend(3)==n3) then
                    call MPI_RECV (qc_x_h2(1:n1-1,j1:j2,xsize(3)),((j2-j1+1)*(n1-1)), MPI_DOUBLE_PRECISION ,start_proc,100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
                    call MPI_SEND (qc_x_h2(1:n1-1,j1:j2,xsize(3)-1),((j2-j1+1)*(n1-1)), MPI_DOUBLE_PRECISION, start_proc, 100, MPI_COMM_WORLD, mpi_err)
                endif

            else
                qc_x_h2(1:n1-1,j1:j2,0)=qc_x_h2(1:n1-1,j1:j2,xsize(3)-1)
                qc_x_h2(1:n1-1,j1:j2,xsize(3))=qc_x_h2(1:n1-1,j1:j2,1)
            end if



            ! Boundaries direction 1 ***************************************

            ! --- Add of a ghost point at Z=0

            do j = j1, j2

                do k = k1, k2
                    qc_x_h2(0, j, k)=qc_x_h2(n1-1, j, k)
                end do
            end do

            ! --- Add of a ghost point at Z=L1

            do j = j1, j2
                do k = k1, k2
                    qc_x_h2(n1, j, k)=qc_x_h2(1, j, k)
                end do
            end do


            ! --------------------------------------------------------------------------------
            !               Filling the boundaries conditions related to the edges
            !               -> Generates manual halo points on the edges of
            !                  the intern workspace/area, without the last
            !                  points at the end of the edges
            !---------------------------------------------------------------------------------

            ! Edge 3 - 2 ***************************************
            if (xsize(3)/=n3) then

                if((xstart(3)==1).and.(xstart(2)==1)) then
                    call MPI_SEND (qc_x_h2(:,0,1),(n1+1), MPI_DOUBLE_PRECISION, end_proc, 100, MPI_COMM_WORLD, mpi_err)
                    call MPI_RECV (qc_x_h2(:,0,0),(n1+1), MPI_DOUBLE_PRECISION ,end_proc,100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
                endif

                if((xend(3)==n3).and.(xstart(2)==1)) then
                    call MPI_RECV (qc_x_h2(:,0,xsize(3)),(n1+1), MPI_DOUBLE_PRECISION ,start_proc,100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
                    call MPI_SEND (qc_x_h2(:,0,xsize(3)-1),(n1+1), MPI_DOUBLE_PRECISION, start_proc, 100, MPI_COMM_WORLD, mpi_err)
                endif

                if ((xstart(3)==1).and.(xend(2)==n2)) then
                    call MPI_SEND (qc_x_h2(:,xsize(2),1),(n1+1), MPI_DOUBLE_PRECISION, end_proc, 100, MPI_COMM_WORLD, mpi_err)
                    call MPI_RECV (qc_x_h2(:,xsize(2),0),(n1+1), MPI_DOUBLE_PRECISION ,end_proc,100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
                endif

                if ((xend(3)==n3).and.(xend(2)==n2)) then
                    call MPI_RECV (qc_x_h2(:,xsize(2),xsize(3)),(n1+1), MPI_DOUBLE_PRECISION ,start_proc,100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
                    call MPI_SEND (qc_x_h2(:,xsize(2),xsize(3)-1),(n1+1), MPI_DOUBLE_PRECISION, start_proc, 100, MPI_COMM_WORLD, mpi_err)
                endif

            else

                qc_x_h2(:,0,0)=qc_x_h2(:,0,xsize(3)-1)
                qc_x_h2(:,0,xsize(3))=qc_x_h2(:,0,1)

                qc_x_h2(:,xsize(2),0)=qc_x_h2(:,xsize(2),xsize(3)-1)
                qc_x_h2(:,xsize(2),xsize(3))=qc_x_h2(:,xsize(2),1)

            end if


            ! Edge 3 - 1 ***************************************

            if (xstart(3)==1) then

                qc_x_h2(0, :, 0)=qc_x_h2(n1-1, :, 0)
                qc_x_h2(n1, :, 0)=qc_x_h2(1, :, 0)

            end if

            if (xend(3)==n3) then
                qc_x_h2(0, :, xsize(3))=qc_x_h2(n1-1, :, xsize(3))
                qc_x_h2(n1, :, xsize(3))=qc_x_h2(1, :, xsize(3))
            end if


            ! Edge 2 - 1 ***************************************

            if ((xstart(2)==1)) then

                do k = k1, k2
                    qc_x_h2(0, 0, :)=qc_x_h2(n1-1, 0, :)
                    qc_x_h2(n1, 0, :)=qc_x_h2(1, 0, :)
                end do

            endif

            if ((xend(2)==n2)) then

                do k = k1, k2
                    qc_x_h2(0, xsize(2), :)=qc_x_h2(n1-1, xsize(2), :)
                    qc_x_h2(n1, xsize(2), :)=qc_x_h2(1, xsize(2), :)
                end do

            endif

            ! --------------------------------------------------------------------------------
            !               Filling the last points boundaries conditions
            !               -> Generates manual halo points at the end of
            !                  the edges
            !---------------------------------------------------------------------------------


            ! Points on 1 - 3 ***************************************

            if ((xstart(1)==1).and.(xstart(2)==1).and.(xstart(3)==1)) then

                qc_x_h2(0, 0, 0)=qc_x_h2(n1-1, 0, 0)

            end if


            if ((xend(1)==n1).and.(xstart(2)==1).and.(xstart(3)==1)) then

                qc_x_h2(xsize(1), 0, 0)=qc_x_h2(1, 0, 0)

            end if

            if ((xstart(1)==1).and.(xend(2)==n2).and.(xstart(3)==1)) then
                qc_x_h2(0, xsize(2), 0)=qc_x_h2(n1-1, xsize(2), 0)
            end if


            if ((xend(1)==n1).and.(xend(2)==n2).and.(xstart(3)==1)) then
                qc_x_h2(xsize(1), xsize(2), 0)=qc_x_h2(1, xsize(2), 0)
            end if




            if ((xstart(1)==1).and.(xstart(2)==1).and.(xend(3)==n3)) then
                qc_x_h2(0, 0, xsize(3))=qc_x_h2(n1-1, 0, xsize(3))
            end if


            if ((xend(1)==n1).and.(xstart(2)==1).and.(xend(3)==n3)) then
                qc_x_h2(xsize(1), 0, xsize(3))=qc_x_h2(1, 0, xsize(3))
            end if

            if ((xstart(1)==1).and.(xend(2)==n2).and.(xend(3)==n3)) then
                qc_x_h2(0, xsize(2), xsize(3))=qc_x_h2(n1-1, xsize(2), xsize(3))
            end if


            if ((xend(1)==n1).and.(xend(2)==n2).and.(xend(3)==n3)) then
                qc_x_h2(xsize(1), xsize(2), xsize(3))=qc_x_h2(1, xsize(2), xsize(3))
            end if

        end subroutine treat_boundaries

    end subroutine perform_fields_at_bubble


    subroutine move_bubble(sub_step, dt_tbl)
        implicit none

        integer, intent(in)                 :: sub_step
        real*8, dimension(8), intent(out)   :: dt_tbl
        real*8, dimension(8)                :: t_tmp
        real*8                              :: dt2

        real*8, parameter               :: Cd=160.d0, d=1.71d-4, divro=833.d0, g=-13380d0, Cl=0.5d0
        real*8, dimension(3)            :: u0, u1, v0, v1, om0, om1, alpha, alpha0, v1save
        real*8, dimension(3)            :: vec_u0, vec_u1, vec_v0, vec_v1
        real*8                          :: beta, gam, com1, com2, com4
        integer                         :: i, j, mpierr
        integer, save                   :: it=1



        dt2=0.d0
        dt2=dt/sub_step

        com1=1.d0
        com2=1.d0
        com4=1.d0

        open(15, file="BUBBLE", position="append")
        open(16, file="BUBBLE.csv", position="append")

        do i = 1, bubble_cpt

            if ((bubl(i)%x1<=Z(1)).or.(bubl(i)%x1>=Z(n1))) then
                bubl(i)%v1=-bubl(i)%v1
            end if

            if ((bubl(i)%x2<=Y(1)).or.(bubl(i)%x2>=Y(n2))) then
                bubl(i)%v2=-bubl(i)%v2
            end if

            if ((bubl(i)%x3<=X(1)).or.(bubl(i)%x3>=X(n3))) then
                bubl(i)%v3=-bubl(i)%v3
            end if

            v0(1) = bubl(i)%v1
            v0(2) = bubl(i)%v2
            v0(3) = bubl(i)%v3

            bubl(i)%x1=bubl(i)%x1+v0(1)*dt2
            bubl(i)%x2=bubl(i)%x2+v0(2)*dt2
            bubl(i)%x3=bubl(i)%x3+v0(3)*dt2

            write(15,*)v0(2)
            write(16,*)dt2, ',', v0(2)

        end do

        call reassign_particle
        call perform_bubble_action


        close(15)
        close(16)

    contains

        subroutine solve_CN(al, be, u, v)
            implicit none

            real*8, intent(in)  ::  be, al, u
            real*8, intent(out) ::  v

            real*8              ::  a, b, c, x1, x2

            integer, parameter  :: fid=125

            v=0.d0

            if ((al-u)<0.d0) then
                a = be
                b = - be*u*2.d0 - 1.d0
                c = al + be*u*u
                call solve_square(a, b, c, x1, x2)

                v = min(x1, x2)

            elseif ((al-u)>0.d0)then
                a = - be
                b = be*u*2.d0 - 1.d0
                c = al - be*u*u
                call solve_square(a, b, c, x1, x2)

                v = max(x1, x2)

            endif

        end subroutine solve_CN

        subroutine solve_square(a, b, c, x1, x2)
            implicit none

            real*8, intent(in)  ::  a, b, c
            real*8, intent(out) ::  x1, x2

            real*8              ::  delta

            delta = b*b - 4.d0*a*c

            if (delta==0.d0) then
                x1 = (-b)/(2.d0*a)
                x2 = x1

            elseif (delta>0.d0) then
                x1 = (-b-dsqrt(delta))/(2.d0*a)
                x2 = (-b+dsqrt(delta))/(2.d0*a)
            endif

            if (a==0.d0) then
                x1 = -c/b
                x2 = x1
            endif

        end subroutine solve_square


        subroutine vect(A, B, AB)
            implicit none

            real*8, dimension(3), intent(in)    :: A, B
            real*8, dimension(3), intent(out)   :: AB

            AB(1)=A(2)*B(3)-A(3)*B(2)
            AB(2)=A(3)*B(1)-A(1)*B(3)
            AB(3)=A(1)*B(2)-A(2)*B(1)

        end subroutine vect


    end subroutine move_bubble


    subroutine reassign_particle() ! Checks the assignment of each bubble. A processor is able to send a bubble to another processor if it's not its property.
        implicit none

        real*8, dimension(3)    :: pt
        integer                 :: i, proc_num, mpi_err
        integer                 :: send, send_glob

        send=0

        do i = 1, bubble_cpt
            pt(1) = bubl(i)%x1
            pt(2) = bubl(i)%x2
            pt(3) = bubl(i)%x3
            call get_processor_number(pt, 1, proc_num)

            if ((proc_num/=nrank).and.(proc_num/=-1)) then
                send=1
                bubl(i)%proc = proc_num
                call post_bubble(proc_num, bubl(i))
                call remove_bubble(i)

            elseif (proc_num==-1) then
                call remove_bubble(i)
            end if

        end do

        call MPI_ALLREDUCE (send, send_glob, 1, MPI_INTEGER, MPI_SUM , MPI_COMM_WORLD , mpi_err)
        if (send_glob>0) call send_bubble(bubl, bubble_cpt)


        send=0
        bubble_tmp_cpt=0

        open(17, file="REASSIGN", position="append")

        do i = 1, bubble_cpt

            pt(1) = bubl(i)%x1
            pt(2) = bubl(i)%x2
            pt(3) = bubl(i)%x3

            if (xend(1)/=n1) then

                call get_processor_number((/Zc(xend(1)+1), pt(2), pt(3) /), 1, proc_num)

                if ((proc_num/=-1).and.(pt(1)>=Zc(xend(1))).and.(pt(1)<=Z(xend(1)+1))) then
                    send=1
                    call post_bubble(proc_num, bubl(i))
                endif

            end if

            if (xend(2)/=n2) then

                call get_processor_number((/pt(1), Yc(xend(2)+1), pt(3) /), 1, proc_num)
                !                write(17,*) 'X2: ', nrank, proc_num, pt(2), Yc(xend(2)), Y(xend(2)+1), Yc(xend(2)+1)!, xpos_end_glob(nrank, 2)

                if ((proc_num/=-1).and.(pt(2)>=Yc(xend(2))).and.(pt(2)<=Y(xend(2)+1))) then
                    write(17,*) 'X2: ', nrank, proc_num, pt(2), Yc(xend(2)), Y(xend(2)+1), Yc(xend(2)+1)!, xpos_end_glob(nrank, 2)
                    send=1
                    call post_bubble(proc_num, bubl(i))
                !                    write(17,*) send_list(proc_num, send_nb(proc_num), :)
                endif

            end if

            if (xend(3)/=n3) then

                call get_processor_number((/pt(1), pt(2), Xc(xend(3)+1) /), 1, proc_num)

                if ((proc_num/=-1).and.(pt(3)>=Xc(xend(3))).and.(pt(3)<=X(xend(3)+1))) then
                    write(17,*) 'X3: ', nrank, proc_num, pt(3), Xc(xend(3)), X(xend(3)+1), Xc(xend(3)+1)!, xpos_end_glob(nrank, 3)
                    send=1
                    call post_bubble(proc_num, bubl(i))
                !                    write(17,*) send_list(proc_num, send_nb(proc_num), :)
                endif

            end if

        end do

        close(17)

        !        write(*,*) '/////////////////////////////////////////// BUBL SL:', nrank, sum(abs(send_list)), bubble_tmp_cpt

        call MPI_ALLREDUCE (send, send_glob, 1, MPI_INTEGER, MPI_SUM , MPI_COMM_WORLD , mpi_err)
        if (send_glob>0) call send_bubble(bubl_tmp, bubble_tmp_cpt)

    !        write(*,*) '///////////////////////////////////////////buble_tmp_cpt:', nrank, bubble_tmp_cpt
    !        write(*,*) '////////////', xend(2), xend(3)

    end subroutine reassign_particle


    subroutine associate_volume_control(import, direction, volume, no_bound_err)
        implicit none

        type(bubble), intent(in)            :: import
        integer, intent(in)                 :: direction
        integer, dimension(3), intent(out)  :: volume
        logical, intent(out)                :: no_bound_err

        real*8, dimension(3)                :: pt
        integer                             :: proc_num
        real*8                              :: x0
        integer                             :: i
        integer                             :: i1, i2

        no_bound_err=.true.
        volume=0

        pt(1) = import%x1
        pt(2) = import%x2
        pt(3) = import%x3


        if (direction==1) then

            x0 = Zc(1)
            volume(1) = floor(((pt(1)-x0)/dx1))+2

            i2 = xend(2)
            if (xend(2)==n2) i2 = xend(2)-1
            do i = xstart(2), i2
                if ((pt(2)>=Y(i)).and.(pt(2)<=Y(i+1))) then
                    volume(2) = i
                end if
            end do

            volume(3) = floor(pt(3)/dx3) + 1

        end if


        if (direction==2) then

            volume(1) = floor(pt(1)/dx1) + 1

            if (xstart(2)==1) then
                if ((pt(2)>=Y(xstart(2))).and.(pt(2)<=Yc(xstart(2)))) then
                    volume(2) = 1
                else

                    do i = xstart(2)+1, xend(2)
                        if ((pt(2)>=Yc(i-1)).and.(pt(2)<=Yc(i))) then
                            volume(2) = i
                        end if
                    end do

                end if

            else
                i2 = xend(2)
                if (xend(2)==n2) i2 = xend(2)-1
                do i = xstart(2), i2
                    if ((pt(2)>=Yc(i-1)).and.(pt(2)<=Yc(i))) then
                        volume(2) = i
                    end if
                end do
            endif

            volume(3) = floor(pt(3)/dx3) + 1

        end if


        if (direction==3) then

            volume(1) = floor(pt(1)/dx1) + 1

            i2 = xend(2)
            if (xend(2)==n2) i2 = xend(2)-1
            do i = xstart(2), i2
                if ((pt(2)>=Y(i)).and.(pt(2)<=Y(i+1))) then
                    volume(2) = i
                end if
            end do

            x0 = Xc(1)
            volume(3) = floor((pt(3)-x0)/dx3) + 2

        end if

        if (volume(1)>xend(1)) Then
            volume(1)=-1
            no_bound_err=.false.
        endif
        if (volume(2)>xend(2)) Then
            volume(2)=-1
            no_bound_err=.false.
        endif
        if (volume(3)>xend(3)) Then
            volume(3)=-1
            no_bound_err=.false.
        endif
        if (volume(1)<xstart(1)) Then
            volume(1)=-1
            no_bound_err=.false.
        endif
        if (volume(2)<xstart(2)) Then
            volume(2)=-1
            no_bound_err=.false.
        endif
        if (volume(3)<xstart(3)) Then
            volume(3)=-1
            no_bound_err=.false.
        endif

    end subroutine associate_volume_control

    subroutine perform_bubble_action()
        implicit none

        integer, dimension(3)   :: volume
        integer                 :: i, j , k
        logical                 :: no_bound_err

        fb1=0.d0
        fb2=0.d0
        fb3=0.d0

        do i = 1, bubble_cpt
            call associate_volume_control(bubl(i), 1, volume, no_bound_err)
            if (no_bound_err) fb1(volume(1), volume(2), volume(3))=fb1(volume(1), volume(2), volume(3))+1

            call associate_volume_control(bubl(i), 2, volume, no_bound_err)
            if (no_bound_err) fb2(volume(1), volume(2), volume(3))=fb2(volume(1), volume(2), volume(3))+1

            call associate_volume_control(bubl(i), 3, volume, no_bound_err)
            if (no_bound_err) fb3(volume(1), volume(2), volume(3))=fb3(volume(1), volume(2), volume(3))+1
        end do

        do i = 1, bubble_tmp_cpt
            call associate_volume_control(bubl_tmp(i), 1, volume, no_bound_err)
            if (no_bound_err) fb1(volume(1), volume(2), volume(3))=fb1(volume(1), volume(2), volume(3))+1

            call associate_volume_control(bubl_tmp(i), 2, volume, no_bound_err)
            if (no_bound_err) fb2(volume(1), volume(2), volume(3))=fb2(volume(1), volume(2), volume(3))+1

            call associate_volume_control(bubl_tmp(i), 3, volume, no_bound_err)
            if (no_bound_err) fb3(volume(1), volume(2), volume(3))=fb3(volume(1), volume(2), volume(3))+1
        end do

    end subroutine perform_bubble_action

    subroutine test_associate(direction)
        implicit none

        integer, intent(in)     :: direction

        type(bubble)            :: bubble_test
        integer, dimension(3)   :: volume_test
        integer                 :: ii, ji, ki
        real*8                  :: i, j, k
        logical                 :: no_bound_err


        i=0.d0
        j=0.d0
        k=0.d0

        ii=0
        ji=0
        ki=0

        bubble_test%v1=0.d0   ! TEST
        bubble_test%v2=0.d0   ! TEST
        bubble_test%v3=0.d0   ! TEST

        bubble_test%u01=0.d0   ! TEST
        bubble_test%u02=0.d0   ! TEST
        bubble_test%u03=0.d0   ! TEST

        bubble_test%u11=0.d0   ! TEST
        bubble_test%u12=0.d0   ! TEST
        bubble_test%u13=0.d0   ! TEST

        bubble_test%om01=0.d0   ! TEST
        bubble_test%om02=0.d0   ! TEST
        bubble_test%om03=0.d0   ! TEST

        bubble_test%om11=0.d0   ! TEST
        bubble_test%om12=0.d0   ! TEST
        bubble_test%om13=0.d0   ! TEST

        bubble_test%a1=0.d0   ! TEST
        bubble_test%a2=0.d0   ! TEST
        bubble_test%a3=0.d0   ! TEST

        bubble_test%b1=0.d0   ! TEST
        bubble_test%b2=0.d0   ! TEST
        bubble_test%b3=0.d0   ! TEST

        bubble_test%c1=0.d0   ! TEST
        bubble_test%c2=0.d0   ! TEST
        bubble_test%c3=0.d0   ! TEST

        bubble_test%d=0.d0   ! TEST

        bubble_test%proc=nrank


        open(29, file="VOLUME", position="append")
        write(29,*) 'nrank:', nrank

        if (direction==1) then

            call random_number(i)
            ii = floor(xstart(1) + i*(xend(1)-xstart(1)))
            if (ii==1) i = (Z(ii)+Zc(ii))*0.5d0
            if (ii>1) i = (Zc(ii-1)+Zc(ii))*0.5d0

            call random_number(j)
            ji = floor(xstart(2) + j*(xend(2)-xstart(2)))
            j = (Y(ji)+Y(ji+1))*0.5d0

            call random_number(k)
            ki = floor(xstart(3) + k*(xend(3)-xstart(3)))
            k = (X(ki)+X(ki+1))*0.5d0

            bubble_test%x1=i
            bubble_test%x2=j
            bubble_test%x3=k

            call associate_volume_control(bubble_test, 1, volume_test, no_bound_err)

        end if

        if (direction==2) then

            call random_number(i)
            ii = floor(xstart(1) + i*(xend(1)-xstart(1)))
            i = (Z(ii)+Z(ii+1))*0.5d0

            call random_number(j)
            ji = floor(xstart(2) + j*(xend(2)-xstart(2)))
            if (ji==1) j = (Y(ji)+Yc(ji))*0.5d0
            if (ji>1) j = (Yc(ji-1)+Yc(ji))*0.5d0

            call random_number(k)
            ki = floor(xstart(3) + k*(xend(3)-xstart(3)))
            k = (X(ki)+X(ki+1))*0.5d0

            bubble_test%x1=i
            bubble_test%x2=j
            bubble_test%x3=k

            call associate_volume_control(bubble_test, 2, volume_test, no_bound_err)

        end if

        if (direction==3) then

            call random_number(i)
            ii = floor(xstart(1) + i*(xend(1)-xstart(1)))
            i = (Z(ii)+Z(ii+1))*0.5d0


            call random_number(j)
            ji = floor(xstart(2) + j*(xend(2)-xstart(2)))
            j = (Y(ji)+Y(ji+1))*0.5d0

            call random_number(k)
            ki = floor(xstart(3) + k*(xend(3)-xstart(3)))
            if (ki==1) k = (X(ki)+Xc(ki))*0.5d0
            if (ki>1) k = (Xc(ki-1)+Xc(ki))*0.5d0

            bubble_test%x1=i
            bubble_test%x2=j
            bubble_test%x3=k

            call associate_volume_control(bubble_test, 3, volume_test, no_bound_err)

        end if

        if ((volume_test(1)==ii).and.(volume_test(2)==ji).and.(volume_test(3)==ki)) then
            write(29,*) 'Test dans la direction:', direction
            write(29,*) 'OK'
        else
            write(29,*) '--------------------------------------------------------------------------'
            write(29,*) 'Test dans la direction:', direction
            write(29,*) '!!!ERREUR!!!'
            write(29,*) 'coordonnées:', i, j, k
            write(29,*) 'Z(1), Zc(1):', Z(1), Zc(1)
            write(29,*) 'dx1, dx3:', dx1, dx3
            write(29,*) 'bornes 1:', xstart(1), xend(1)
            write(29,*) 'bornes 2:', xstart(2), xend(2)
            write(29,*) 'bornes 3:', xstart(3), xend(3)
            write(29,*) ii, ji, ki
            write(29,*) volume_test
            write(29,*) '--------------------------------------------------------------------------'
        endif

        write(29,*) '***************************************************************************************'
        write(29,*) '***************************************************************************************'
        close(29)

    end subroutine test_associate

end module Bubble_generator_test
