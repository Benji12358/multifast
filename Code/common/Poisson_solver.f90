
module FFTW3
    use, intrinsic :: iso_c_binding
         include 'fftw3.f03'
end module


module Fourier_in_xz_plane_seq
    implicit none

contains


    subroutine real_to_spectral(Fxz, Fxz_hat)

        use FFTW3
        use mesh

        implicit none
        integer :: i, j, k
        real*8, dimension(n1, n2, n3)      :: Fxz
        complex*16, dimension(n1, n2, n3m/2+1)  :: Fxz_hat
        integer*8 plan

        real*8, dimension(n3m)           :: rin_z
        complex*16, dimension(n3m/2+1)   :: cout_z

        complex*16, dimension(n1m)       :: cin_x
        complex*16, dimension(n1m)       :: cout_x


        call dfftw_plan_dft_r2c_1d(plan, n3m, rin_z(:), cout_z(:), FFTW_ESTIMATE)

        do j = 1, n2m
            do i = 1, n1m

                rin_z(:)=Fxz(i,j,1:n3m)
                call dfftw_execute_dft_r2c(plan, rin_z, cout_z)
                Fxz_hat(i,j,1:n3m/2+1)=cout_z(:)

            end do
        end do

        call dfftw_destroy_plan(plan)



        call dfftw_plan_dft_1d(plan, n1m, cin_x(:), cout_x(:), FFTW_FORWARD, FFTW_ESTIMATE)

        do j = 1, n2m
            do k = 1, n3m/2+1

                cin_x(:)=Fxz_hat(1:n1m,j,k)
                call dfftw_execute_dft(plan, cin_x(:), cout_x(:))
                Fxz_hat(1:n1m,j,k)=cout_x(:)

            end do
        end do

        call dfftw_destroy_plan(plan)

    end subroutine real_to_spectral


    subroutine spectral_to_real(Fxz_hat, Fxz)

        use FFTW3
        use mesh

        implicit none
        integer :: i, j, k
        real*8, dimension(n1, n2, n3)    :: Fxz
        complex*16, dimension(n1, n2, n3m/2+1)  :: Fxz_hat
        integer*8 plan

        complex*16, dimension(n1m)       :: cin_x
        complex*16, dimension(n1m)       :: cout_x

        real*8, dimension(n3m)           :: rout_z
        complex*16, dimension(n3m/2+1)   :: cin_z


        call dfftw_plan_dft_1d(plan, n1m, cin_x, cout_x, FFTW_BACKWARD, FFTW_ESTIMATE)

        do j = 1, n2m
            do k = 1, n3m/2+1

                cin_x=Fxz_hat(1:n1m,j,k)
                call dfftw_execute_dft(plan, cin_x, cout_x)
                Fxz_hat(1:n1m,j,k)=cout_x/n1m

            end do
        end do

        call dfftw_destroy_plan(plan)



        call dfftw_plan_dft_c2r_1d(plan, n3m, cin_z, rout_z, FFTW_ESTIMATE)

        do j = 1, n2m
            do i = 1, n1m

                cin_z=Fxz_hat(i,j,1:n3m/2+1)
                call dfftw_execute_dft_c2r(plan, cin_z, rout_z)
                Fxz(i,j,1:n3m)=rout_z/n3m

            end do
        end do

        call dfftw_destroy_plan(plan)


    end subroutine spectral_to_real

end module Fourier_in_xz_plane_seq


module Fourier_in_xz_plane
    implicit none

contains


    subroutine real_to_spectral(Fxz_z, Fxz_hat_x, sp_decomp)

        use FFTW3
        use mesh
        use decomp_2d

        implicit none
        integer :: i, j, k
        integer*8 plan

        real*8, dimension(n3m)           :: rin_z
        complex*16, dimension(n3m/2+1)   :: cout_z

        complex*16, dimension(n1m)       :: cin_x
        complex*16, dimension(n1m)       :: cout_x


        type(decomp_info) :: sp_decomp
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))      :: Fxz_z

        complex*16, dimension(      &
        sp_decomp%xst(1):sp_decomp%xen(1), &
        sp_decomp%xst(2):sp_decomp%xen(2), &
        sp_decomp%xst(3):sp_decomp%xen(3))   :: Fxz_hat_x

        complex*16, dimension(      &
        sp_decomp%yst(1):sp_decomp%yen(1), &
        sp_decomp%yst(2):sp_decomp%yen(2), &
        sp_decomp%yst(3):sp_decomp%yen(3))   :: Fxz_hat_y

        complex*16, dimension(      &
        sp_decomp%zst(1):sp_decomp%zen(1), &
        sp_decomp%zst(2):sp_decomp%zen(2), &
        sp_decomp%zst(3):sp_decomp%zen(3))   :: Fxz_hat_z

        call dfftw_plan_dft_r2c_1d(plan, n3m, rin_z(:), cout_z(:), FFTW_ESTIMATE)

        !do j = 1, n2m
        do j = sp_decomp%zst(2), min(n2m, sp_decomp%zen(2))

            !do i=1,n1m
            do i=sp_decomp%zst(1), min(n1m, sp_decomp%zen(1))

                rin_z(:)=Fxz_z(i,j,1:n3m)
                call dfftw_execute_dft_r2c(plan, rin_z, cout_z)
                Fxz_hat_z(i,j,1:n3m/2+1)=cout_z(:)

            end do
        end do

        call dfftw_destroy_plan(plan)

        call transpose_z_to_y(Fxz_hat_z, Fxz_hat_y, sp_decomp)
        call transpose_y_to_x(Fxz_hat_y, Fxz_hat_x, sp_decomp)



        call dfftw_plan_dft_1d(plan, n1m, cin_x(:), cout_x(:), FFTW_FORWARD, FFTW_ESTIMATE)

        do j = sp_decomp%xst(2), min(n2m, sp_decomp%xen(2))
            do k = sp_decomp%xst(3), sp_decomp%xen(3)

                cin_x(:)=Fxz_hat_x(1:n1m,j,k)
                call dfftw_execute_dft(plan, cin_x(:), cout_x(:))
                Fxz_hat_x(1:n1m,j,k)=cout_x(:)

            end do
        end do

        call dfftw_destroy_plan(plan)

    end subroutine real_to_spectral


    subroutine spectral_to_real(Fxz_hat_x, Fxz_z, sp_decomp)

        use FFTW3
        use mesh
        use decomp_2d

        implicit none
        integer :: i, j, k
        integer*8 plan

        type(decomp_info) :: sp_decomp

        complex*16, dimension(n1m)       :: cin_x
        complex*16, dimension(n1m)       :: cout_x

        real*8, dimension(n3m)           :: rout_z
        complex*16, dimension(n3m/2+1)   :: cin_z
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))      :: Fxz_z

        complex*16, dimension(      &
        sp_decomp%xst(1):sp_decomp%xen(1), &
        sp_decomp%xst(2):sp_decomp%xen(2), &
        sp_decomp%xst(3):sp_decomp%xen(3))   :: Fxz_hat_x

        complex*16, dimension(      &
        sp_decomp%yst(1):sp_decomp%yen(1), &
        sp_decomp%yst(2):sp_decomp%yen(2), &
        sp_decomp%yst(3):sp_decomp%yen(3))   :: Fxz_hat_y

        complex*16, dimension(      &
        sp_decomp%zst(1):sp_decomp%zen(1), &
        sp_decomp%zst(2):sp_decomp%zen(2), &
        sp_decomp%zst(3):sp_decomp%zen(3))   :: Fxz_hat_z





        call dfftw_plan_dft_1d(plan, n1m, cin_x, cout_x, FFTW_BACKWARD, FFTW_ESTIMATE)

        do j = sp_decomp%xst(2), min(n2m, sp_decomp%xen(2))
            do k = sp_decomp%xst(3), sp_decomp%xen(3)

                cin_x=Fxz_hat_x(1:n1m,j,k)
                call dfftw_execute_dft(plan, cin_x, cout_x)
                Fxz_hat_x(1:n1m,j,k)=cout_x/n1m

            end do
        end do

        call dfftw_destroy_plan(plan)

        call transpose_x_to_y(Fxz_hat_x, Fxz_hat_y, sp_decomp)
        call transpose_y_to_z(Fxz_hat_y, Fxz_hat_z, sp_decomp)



        call dfftw_plan_dft_c2r_1d(plan, n3m, cin_z, rout_z, FFTW_ESTIMATE)

        do j = sp_decomp%zst(2), min(n2m, sp_decomp%zen(2))
            do i=sp_decomp%zst(1), min(n1m, sp_decomp%zen(1))

                cin_z=Fxz_hat_z(i,j,1:n3m/2+1)
                call dfftw_execute_dft_c2r(plan, cin_z, rout_z)
                Fxz_z(i,j,1:n3m)=rout_z/n3m

            end do
        end do

        call dfftw_destroy_plan(plan)


    end subroutine spectral_to_real

end module Fourier_in_xz_plane


module poisson_020_solver

    use mesh
    use decomp_2d

    implicit none

    ! integer, parameter   :: kl=5,ku=5
    ! integer, parameter   :: kl=1,ku=1
    integer   :: kl,ku
    real*8,dimension(:,:), allocatable :: lapY_matrix
    real*8,dimension(:, :, :, :), allocatable::lap_matrix_ik
    real*8, dimension(:), allocatable   :: k2eq1, k2eq3
    real*8, dimension(:,:,:), allocatable   :: k2eq1_x, k2eq3_z, k2eq1_y, k2eq3_y, k2eq1_z
    real*8, dimension(:,:,:), allocatable   :: tmp_k2eq1_x, tmp_k2eq3_z, tmp_k2eq1_y, tmp_k2eq1_z


    type(decomp_info) :: sp_decomp



    public :: solve_Poisson, init
    private

contains

    subroutine prepare_lapY_matrix

        use d2c_from_d1s
        use mesh
        use irregular_derivative_coefficients
        use numerical_methods_settings, only: schemes_configuration
        use schemes_loader

        use workspace_view, only: lapY_file

        implicit none

        integer i
        real*8,dimension (n2,-1:1) :: d1s_coeffs, geom_coeffs, d1_u, d2_u

        real*8 geomA_coeffs(-1:1), geomB_coeff, A(-1:1), B(-1:1)
        logical :: lapY_file_exists

        ! Definition of first derivatives schemes used for scheme derivation
        ! (d2c=d_B(d_A).
        ! A gives the operator d_A coefficients. Same for B with d_B.

        A(1)=1.d0/dx2
        A(-1)=-A(1)

        B=A

        !-----------------

        lapY_matrix=0.d0

        ! Matrix coefficients calculation
        lapY_matrix(1,-1)=0.d0
        lapY_matrix(1,0)=-Yc_to_YcTr_for_D1(1)*Y_to_YTr_for_D1(2)/dx2**2
        lapY_matrix(1,1)= Yc_to_YcTr_for_D1(1)*Y_to_YTr_for_D1(2)/dx2**2

        lapY_matrix(n2m,-1)=Yc_to_YcTr_for_D1(n2m) * Y_to_YTr_for_D1(n2m)/dx2**2
        lapY_matrix(n2m,0)=-Yc_to_YcTr_for_D1(n2m) * Y_to_YTr_for_D1(n2m)/dx2**2
        lapY_matrix(n2m,1)=0.d0


        inquire( file=lapY_file, exist=lapY_file_exists)

        if (lapY_file_exists) then
            call fill_LapY_from_file(n2m-2, lapY_matrix(2:n2m-1,-kl:ku))

        else

            if (schemes_configuration.eq.O2_SCHEMES) then

                call perform_lapY_O2(n2m-2, dx2, lapY_matrix(2:n2m-1,-kl:ku),       &
                Yc_to_YcTr_for_D1(2:n2m-1), Y_to_YTr_for_D1(2:n2m))

            else

                call perform_lapY(n2m-2, dx2, lapY_matrix(2:n2m-1,-kl:ku),       &
                Yc_to_YcTr_for_D1(2:n2m-1), Y_to_YTr_for_D1(2:n2m))


            endif

        endif

        !write(*,*)'Cmatrix'
        !do i = 2, n2m
            !write(*,*)'i=',i
            !write(*,*) lapY_matrix(i,-kl:ku)
        !end do

        return



    contains


        subroutine perform_lapY(n_max, h, C, GC_at_f, GC_at_df)

            use d2c_from_d1s

            implicit none
            integer :: i, n_max, j

            real*8 h, geomB_coeff
            integer A_size(-3:3)
            real *8 GC_at_f(n_max)
            real *8 GC_at_df(n_max+1)
            real*8 A(-3:3,-3:3), B(-3:3), C(n_max, -5:5), geomA_coeffs(-3:3)

            integer, parameter  :: MATRIX_FILE=98

            geomB_coeff=1.d0

            open(MATRIX_FILE, file="matrix_coeffs")


            ! j=1-------------------------------------------------------
            geomA_coeffs=0.d0

            C(1,:)=0.d0
            A=0.d0
            B=0.d0
            A_size=0

            A(-1,1) =  1.d0 / h
            A(-1,-1)=-A(-1,1)
            A_size(-1)=1

            A(1,1) =  9.d0 / (8.d0*h)
            A(1,2) =  -1.d0 / (24.d0*h)
            A(1,-1)=-A(1,1)
            A(1,-2)=-A(1,2)
            A_size(1)=2

            B=A(-1,:)

            geomA_coeffs(-1)=GC_at_df(1)
            geomA_coeffs(1)=GC_at_df(2)
            geomB_coeff=GC_at_f(1)


            call set_d2c_coeffs3(A(-1:1, -2:2), B(-1:1), C(1,-1:2),     &
            A_size(-1:1), 2, 1, geomA_coeffs(-1:1), geomB_coeff)

            ! j=2-------------------------------------------------------
            geomA_coeffs=0.d0

            C(2,:)=0.d0
            A=0.d0
            B=0.d0
            A_size=0

            A(-2,1) =  1.d0 / h
            A(-2,-1)=-A(-2,1)
            A_size(-2)=1

            A(-1,1) =  9.d0 / (8.d0*h)
            A(-1,2) =  -1.d0 / (24.d0*h)
            A(-1,-1)=-A(-1,1)
            A(-1,-2)=-A(-1,2)
            A_size(-1)=2


            A(1:2,1) =  1.189101951200031d0 / h
            A(1:2,2) =  -0.073717642266682d0 / h
            A(1:2,3) =  0.006410195120003d0  / h

            A(1:2,-1)=-A(1:2,1)
            A(1:2,-2)=-A(1:2,2)
            A(1:2,-3)=-A(1:2,3)

            A_size(1:2)=3

            B=A(-1,:)

            B=A(-1,:)

            geomA_coeffs(-2)=GC_at_df(1)
            geomA_coeffs(-1)=GC_at_df(2)
            geomA_coeffs(1)=GC_at_df(3)
            geomA_coeffs(2)=GC_at_df(4)

            geomB_coeff=GC_at_f(2)

            call set_d2c_coeffs3(A(-2:2, -3:3), B(-2:2), C(2,-2:4),     &
            A_size(-2:2), 3, 2, geomA_coeffs(-2:2), geomB_coeff)

            ! j=3-------------------------------------------------------
            geomA_coeffs=0.d0

            C(3,:)=0.d0
            A=0.d0
            B=0.d0
            A_size=0

            A(-3,1) =  1.d0 / h
            A(-3,-1)=-A(-3,1)
            A_size(-3)=1

            A(-2,1) =  9.d0 / (8.d0*h)
            A(-2,2) =  -1.d0 / (24.d0*h)
            A(-2,-1)=-A(-2,1)
            A(-2,-2)=-A(-2,2)
            A_size(-2)=2

            A(-1:3,1) =  1.189101951200031d0 / h
            A(-1:3,2) =  -0.073717642266682d0 / h
            A(-1:3,3) =  0.006410195120003d0  / h

            A(-1:3,-1)=-A(-1:3,1)
            A(-1:3,-2)=-A(-1:3,2)
            A(-1:3,-3)=-A(-1:3,3)

            A_size(-1:3)=3

            B=A(-1,:)

            geomA_coeffs(-3)=GC_at_df(1)
            geomA_coeffs(-2)=GC_at_df(2)
            geomA_coeffs(-1)=GC_at_df(3)
            geomA_coeffs(1)=GC_at_df(4)
            geomA_coeffs(2)=GC_at_df(5)
            geomA_coeffs(3)=GC_at_df(6)

            geomB_coeff=GC_at_f(3)

            call set_d2c_coeffs3(A(-3:3, -3:3), B(-3:3), C(3,-3:5),     &
            A_size(-3:3), 3, 3, geomA_coeffs(-3:3), geomB_coeff)

            ! j=4-------------------------------------------------------
            geomA_coeffs=0.d0

            C(4,:)=0.d0
            A=0.d0
            B=0.d0
            A_size=0

            A(-3,1) =  9.d0 / (8.d0*h)
            A(-3,2) =  -1.d0 / (24.d0*h)
            A(-3,-1)=-A(-3,1)
            A(-3,-2)=-A(-3,2)
            A_size(-3)=2

            A(-2:3,1) =  1.189101951200031d0 / h
            A(-2:3,2) =  -0.073717642266682d0 / h
            A(-2:3,3) =  0.006410195120003d0  / h

            A(-2:3,-1)=-A(-2:3,1)
            A(-2:3,-2)=-A(-2:3,2)
            A(-2:3,-3)=-A(-2:3,3)

            A_size(-2:3)=3

            B=A(-1,:)

            geomA_coeffs(-3)=GC_at_df(2)
            geomA_coeffs(-2)=GC_at_df(3)
            geomA_coeffs(-1)=GC_at_df(4)
            geomA_coeffs(1)=GC_at_df(5)
            geomA_coeffs(2)=GC_at_df(6)
            geomA_coeffs(3)=GC_at_df(7)

            geomB_coeff=GC_at_f(4)

            call set_d2c_coeffs3(A(-3:3, -3:3), B(-3:3), C(4,-4:5), A_size(-3:3), 3, 3,     &
            geomA_coeffs(-3:3), geomB_coeff)

            if(nrank==0) write(MATRIX_FILE,"(10E25.15)")(C(4,j), j=-4,5)


            ! j=5_to_n_max-4-------------------------------------------------------
            geomA_coeffs=0.d0

            A=0.d0
            B=0.d0
            A_size=0

            A(-3:3,1) =  1.189101951200031d0 / h
            A(-3:3,2) =  -0.073717642266682d0 / h
            A(-3:3,3) =  0.006410195120003d0  / h

            A(-3:3,-1)=-A(-3:3,1)
            A(-3:3,-2)=-A(-3:3,2)
            A(-3:3,-3)=-A(-3:3,3)

            A_size(-3:3)=3

            B=A(-1,:)

            do i = 5, n_max-4

                C(i,:)=0.d0

                geomA_coeffs(-3)=GC_at_df(i-2)
                geomA_coeffs(-2)=GC_at_df(i-1)
                geomA_coeffs(-1)=GC_at_df(i)
                geomA_coeffs(1)=GC_at_df(i+1)
                geomA_coeffs(2)=GC_at_df(i+2)
                geomA_coeffs(3)=GC_at_df(i+3)

                geomB_coeff=GC_at_f(i)


                call set_d2c_coeffs3(A(-3:3, -3:3), B(-3:3), C(i,-5:5),     &
                A_size(-3:3), 3, 3, geomA_coeffs(-3:3), geomB_coeff)

                if(nrank==0) write(MATRIX_FILE,"(11E25.15)")(C(i,j), j=-5,5)

            end do



            ! j=nmax----------------------------------------------------
            geomA_coeffs=0.d0

            C(n_max,:)=0.d0
            A=0.d0
            B=0.d0
            A_size=0

            A(1,1) =  1.d0 / h
            A(1,-1)=-A(1,1)
            A_size(1)=1

            A(-1,1) =  9.d0 / (8.d0*h)
            A(-1,2) =  -1.d0 / (24.d0*h)
            A(-1,-1)=-A(-1,1)
            A(-1,-2)=-A(-1,2)
            A_size(-1)=2

            B=A(1,:)

            geomA_coeffs(-1)=GC_at_df(n_max)
            geomA_coeffs(1)=GC_at_df(n_max+1)

            geomB_coeff=GC_at_f(n_max)

            call set_d2c_coeffs3(A(-1:1, -2:2), B(-1:1), C(n_max,-2:1), &
            A_size(-1:1), 2, 1, geomA_coeffs(-1:1), geomB_coeff)

            if(nrank==0) write(MATRIX_FILE,"(4E25.15)")(C(n_max,j), j=-2,1)


            ! j=nmax-1 -------------------------------------------------
            geomA_coeffs=0.d0

            C(n_max-1,:)=0.d0
            A=0.d0
            B=0.d0
            A_size=0

            A(2,1) =  1.d0 / h
            A(2,-1)=-A(2,1)
            A_size(2)=1

            A(1,1) =  9.d0 / (8.d0*h)
            A(1,2) =  -1.d0 / (24.d0*h)
            A(1,-1)=-A(1,1)
            A(1,-2)=-A(1,2)
            A_size(1)=2


            A(-2:-1,1) =  1.189101951200031d0 / h
            A(-2:-1,2) =  -0.073717642266682d0 / h
            A(-2:-1,3) =  0.006410195120003d0  / h

            A(-2:-1,-1)=-A(-2:-1,1)
            A(-2:-1,-2)=-A(-2:-1,2)
            A(-2:-1,-3)=-A(-2:-1,3)

            A_size(-2:-1)=3

            B=A(1,:)

            geomA_coeffs(-2)=GC_at_df(n_max-2)
            geomA_coeffs(-1)=GC_at_df(n_max-1)
            geomA_coeffs(1)=GC_at_df(n_max)
            geomA_coeffs(2)=GC_at_df(n_max+1)

            geomB_coeff=GC_at_f(n_max-1)

            call set_d2c_coeffs3(A(-2:2, -3:3), B(-2:2), C(n_max-1,-4:2), &
            A_size(-2:2), 3, 2, geomA_coeffs(-2:2), geomB_coeff)

            if(nrank==0) write(MATRIX_FILE,"(7E25.15)")(C(n_max-1,j), j=-4,2)


            ! j=nmax-2--------------------------------------------------
            geomA_coeffs=0.d0

            C(n_max-2,:)=0.d0
            A=0.d0
            B=0.d0
            A_size=0

            A(3,1) =  1.d0 / h
            A(3,-1)=-A(3,1)
            A_size(3)=1

            A(2,1) =  9.d0 / (8.d0*h)
            A(2,2) =  -1.d0 / (24.d0*h)
            A(2,-1)=-A(2,1)
            A(2,-2)=-A(2,2)
            A_size(2)=2

            A(-3:1,1) =  1.189101951200031d0 / h
            A(-3:1,2) =  -0.073717642266682d0 / h
            A(-3:1,3) =  0.006410195120003d0  / h

            A(-3:1,-1)=-A(-3:1,1)
            A(-3:1,-2)=-A(-3:1,2)
            A(-3:1,-3)=-A(-3:1,3)

            A_size(-3:1)=3

            B=A(1,:)

            geomA_coeffs(-3)=GC_at_df(n_max-4)
            geomA_coeffs(-2)=GC_at_df(n_max-3)
            geomA_coeffs(-1)=GC_at_df(n_max-2)
            geomA_coeffs(1)=GC_at_df(n_max-1)
            geomA_coeffs(2)=GC_at_df(n_max)
            geomA_coeffs(3)=GC_at_df(n_max+1)

            geomB_coeff=GC_at_f(n_max-2)

            call set_d2c_coeffs3(A(-3:3, -3:3), B(-3:3), C(n_max-2,-5:3), &
            A_size(-3:3), 3, 3, geomA_coeffs(-3:3), geomB_coeff)

            if(nrank==0) write(MATRIX_FILE,"(9E25.15)")(C(n_max-2,j), j=-5,3)

            ! j=nmax-3--------------------------------------------------
            geomA_coeffs=0.d0

            C(n_max-3,:)=0.d0
            A=0.d0
            B=0.d0
            A_size=0

            A(3,1) =  9.d0 / (8.d0*h)
            A(3,2) =  -1.d0 / (24.d0*h)
            A(3,-1)=-A(3,1)
            A(3,-2)=-A(3,2)
            A_size(3)=2

            A(-3:2,1) =  1.189101951200031d0 / h
            A(-3:2,2) =  -0.073717642266682d0 / h
            A(-3:2,3) =  0.006410195120003d0  / h

            A(-3:2,-1)=-A(-3:2,1)
            A(-3:2,-2)=-A(-3:2,2)
            A(-3:2,-3)=-A(-3:2,3)

            A_size(-3:2)=3

            B=A(1,:)

            geomA_coeffs(-3)=GC_at_df(n_max-5)
            geomA_coeffs(-2)=GC_at_df(n_max-4)
            geomA_coeffs(-1)=GC_at_df(n_max-3)
            geomA_coeffs(1)=GC_at_df(n_max-2)
            geomA_coeffs(2)=GC_at_df(n_max-1)
            geomA_coeffs(3)=GC_at_df(n_max)

            geomB_coeff=GC_at_f(n_max-3)

            call set_d2c_coeffs3(A(-3:3, -3:3), B(-3:3), C(n_max-3,-5:4), &
            A_size(-3:3), 3, 3, geomA_coeffs(-3:3), geomB_coeff)

            if(nrank==0) write(MATRIX_FILE,"(10E25.15)")(C(n_max-3,j), j=-5,4)

            close(MATRIX_FILE)

            return

        end subroutine

        subroutine perform_lapY_O2(n_max, h, C, GC_at_f, GC_at_df)

            use d2c_from_d1s

            implicit none
            integer :: i, n_max, j

            real*8 h, geomB_coeff
            integer A_size, B_size
            real *8 GC_at_f(n_max)
            real *8 GC_at_df(n_max+1)
            real*8 A(-1:1), B(-1:1), C(n_max, -1:1), geomA_coeffs(-1:1)

            integer, parameter  :: MATRIX_FILE=98

            geomB_coeff=1.d0

            open(MATRIX_FILE, file="matrix_coeffs")

            ! j=1 to n_max-------------------------------------------------------
            geomA_coeffs=0.d0

            A=0.d0
            B=0.d0

            A(1) =  1.d0 / h
            A(-1)=-A(1)

            B=A

            A_size = 1
            B_size = 1

            do i = 1, n_max

                C(i,:) = 0.d0
                geomA_coeffs(-1)=GC_at_df(i)
                geomA_coeffs(1)=GC_at_df(i+1)
                geomB_coeff=GC_at_f(i)

                geomB_coeff=GC_at_f(i)

                ! In contrary of general routine using DRP, here, we use set_d2c_coeffs and not set_d2c_coeffs3
                call set_d2c_coeffs(A(-1:1), B(-1:1), C(i,-1:1),     &
                    A_size, B_size, geomA_coeffs(-1:1), geomB_coeff)

                if(nrank==0) write(MATRIX_FILE,"(11E25.15)")(C(i,j), j=-1,1)

            end do
            ! ---------------------------------------------------------------------

            close(MATRIX_FILE)

            return

        end subroutine


        subroutine fill_LapY_from_file(n_max, C)

            use d2c_from_d1s

            implicit none
            integer :: i, n_max, j

            real*8 h, geomB_coeff
            integer A_size(-3:3)
            real *8 GC_at_f(n_max)
            real *8 GC_at_df(n_max+1)
            real*8 A(-3:3,-3:3), B(-3:3), C(n_max, -5:5), geomA_coeffs(-3:3)

            integer, parameter  :: MATRIX_FILE=97, MATRIX_FILE2=98

            geomB_coeff=1.d0

            open(MATRIX_FILE, file=lapY_file)
            !open(MATRIX_FILE2, file="matrix_coeffs_test")


            ! j=1-------------------------------------------------------
            read(MATRIX_FILE,*)(C(1, j), j=-1,2)
            !write(MATRIX_FILE2,"(4E25.15)")(C(1, j), j=-1,2)

            ! j=2-------------------------------------------------------
            read(MATRIX_FILE,*)(C(2,j), j=-2,4)
            !write(MATRIX_FILE2,"(7E25.15)")(C(2,j), j=-2,4)

            ! j=3-------------------------------------------------------
            read(MATRIX_FILE,*)(C(3,j), j=-3,5)
            !write(MATRIX_FILE2,"(9E25.15)")(C(3,j), j=-3,5)

            ! j=4-------------------------------------------------------
            read(MATRIX_FILE,*)(C(4,j), j=-4,5)
            !write(MATRIX_FILE2,"(10E25.15)")(C(4,j), j=-4,5)


            ! j=5-------------------------------------------------------

            do i = 5, n_max-4
                read(MATRIX_FILE,*)(C(i,j), j=-5,5)
                !write(MATRIX_FILE2,"(11E25.15)")(C(i,j), j=-5,5)
            end do

            ! j=nmax----------------------------------------------------
            read(MATRIX_FILE,*)(C(n_max,j), j=-2,1)
            !write(MATRIX_FILE2,"(4E25.15)")(C(n_max,j), j=-2,1)


            ! j=nmax-1 -------------------------------------------------
            read(MATRIX_FILE,*)(C(n_max-1,j), j=-4,2)
            !write(MATRIX_FILE2,"(7E25.15)")(C(n_max-1,j), j=-4,2)


            ! j=nmax-2--------------------------------------------------
            read(MATRIX_FILE,*)(C(n_max-2,j), j=-5,3)
            !write(MATRIX_FILE2,"(9E25.15)")(C(n_max-2,j), j=-5,3)

            ! j=nmax-3--------------------------------------------------
            read(MATRIX_FILE,*)(C(n_max-3,j), j=-5,4)
            !write(MATRIX_FILE2,"(10E25.15)")(C(n_max-3,j), j=-5,4)

            close(MATRIX_FILE)
            !close(MATRIX_FILE2)

            return

        end subroutine

    end subroutine

    subroutine init
        use IBM_settings
        use numerical_methods_settings
        use schemes_loader

        implicit none

        integer n3mh,k,i,j

        if (schemes_configuration.eq.O2_SCHEMES) then
            kl = 1
            ku = 1
        else
            kl = 5
            ku = 5
        endif

        n3mh=n3m/2+1
        call decomp_info_init(n1, n2, n3mh, sp_decomp)

        allocate(lapY_matrix(n2,-kl:ku))
        allocate(lap_matrix_ik(sp_decomp%yst(1):sp_decomp%yen(1),sp_decomp%yst(3):sp_decomp%yen(3),n2-1,-kl:ku))

        allocate(tmp_k2eq1_y(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
        allocate(tmp_k2eq1_x(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
        allocate(tmp_k2eq3_z(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
        allocate(tmp_k2eq1_z(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))


        allocate(k2eq1_z(sp_decomp%zst(1):sp_decomp%zen(1),sp_decomp%zst(2):sp_decomp%zen(2),sp_decomp%zst(3):sp_decomp%zen(3)))
        allocate(k2eq3_z(sp_decomp%zst(1):sp_decomp%zen(1),sp_decomp%zst(2):sp_decomp%zen(2),sp_decomp%zst(3):sp_decomp%zen(3)))
        allocate(k2eq1_y(sp_decomp%yst(1):sp_decomp%yen(1),sp_decomp%yst(2):sp_decomp%yen(2),sp_decomp%yst(3):sp_decomp%yen(3)))
        allocate(k2eq3_y(sp_decomp%yst(1):sp_decomp%yen(1),sp_decomp%yst(2):sp_decomp%yen(2),sp_decomp%yst(3):sp_decomp%yen(3)))
        allocate(k2eq1_x(sp_decomp%xst(1):sp_decomp%xen(1),sp_decomp%xst(2):sp_decomp%xen(2),sp_decomp%xst(3):sp_decomp%xen(3)))
        
        allocate(k2eq1(n1))
        allocate(k2eq3(n3))

        call perform_K2eq_for_ik

        do j=sp_decomp%xst(2),sp_decomp%xen(2)
            do k=sp_decomp%xst(3),sp_decomp%xen(3)
                k2eq1_x(:,j,k) = k2eq1(sp_decomp%xst(1):sp_decomp%xen(1))
            enddo
        enddo

        do j=sp_decomp%zst(2),sp_decomp%zen(2)
            do i=sp_decomp%zst(1),sp_decomp%zen(1)
                k2eq3_z(i,j,:) = k2eq3(sp_decomp%zst(3):sp_decomp%zen(3))
            enddo
        enddo

        call transpose_x_to_y(k2eq1_x,k2eq1_y,sp_decomp)
        call transpose_z_to_y(k2eq3_z,k2eq3_y,sp_decomp)


        call prepare_lapY_matrix


        ! Zero wave number i, k ----------------------------------
        if (sp_decomp%yst(1)==1) then

            if (sp_decomp%yst(3)==1) then

                lap_matrix_ik(1,1,1,-kl:-1)=0.d0
                lap_matrix_ik(1,1,1,0)=1.d0
                lap_matrix_ik(1,1,1,1:ku)=0.d0

                do j=2,n2m
                    lap_matrix_ik(1,1,j,-kl:-1)=lapY_matrix(j,-kl:-1)
                    ! lap_matrix_ik(1,1,j,0)=lapY_matrix(j,0)-k2eq1(1)-k2eq3(1)
                    lap_matrix_ik(1,1,j,0)=lapY_matrix(j,0)-k2eq1_y(1,j,1)-k2eq3_y(1,j,1)
                    lap_matrix_ik(1,1,j,1:ku)=lapY_matrix(j,1:ku)
                enddo

            end if

        end if

        ! Non zero wave number i ----------------------------------
        if (sp_decomp%yst(3)==1) then

            do j=1,n2m
                do i=max(2, sp_decomp%yst(1)), min(n1m, sp_decomp%yen(1))

                    lap_matrix_ik(i,1,j,-kl:-1)=lapY_matrix(j,-kl:-1)
                    ! lap_matrix_ik(i,1,j,0)=lapY_matrix(j,0)-k2eq1(i)-k2eq3(1)
                    lap_matrix_ik(i,1,j,0)=lapY_matrix(j,0)-k2eq1_y(i,j,1)-k2eq3_y(i,j,1)
                    lap_matrix_ik(i,1,j,1:ku)=lapY_matrix(j,1:ku)

                enddo
            enddo

        end if

        ! Non zero wave number k ----------------------------------
        do k=max(2, sp_decomp%yst(3)), min(n3mh, sp_decomp%yen(3))
            do j=1,n2m
                do i=sp_decomp%yst(1), min(n1m, sp_decomp%yen(1))

                    lap_matrix_ik(i,k,j,-kl:-1)=lapY_matrix(j,-kl:-1)
                    ! lap_matrix_ik(i,k,j,0)=lapY_matrix(j,0)-k2eq1(i)-k2eq3(k)
                    lap_matrix_ik(i,k,j,0)=lapY_matrix(j,0)-k2eq1_y(i,j,k)-k2eq3_y(i,j,k)
                    lap_matrix_ik(i,k,j,1:ku)=lapY_matrix(j,1:ku)

                enddo
            enddo
        enddo

        return

    contains

        subroutine perform_K2eq_for_ik

            use mathematical_constants
            use d2c_from_d1s
            use schemes_interface

            implicit none

            include 'fftw3.f'
            integer k,i

            do k=1,n3m
                k2eq3(k)=K2eq_of_D1s_twice(D1s, n3, k, L3/PI, .true.)
            enddo

            do i=1,n1m
                k2eq1(i)=K2eq_of_D1s_twice(D1s, n1, i, L1/PI, .true.)
            enddo

            return

        end subroutine

    end subroutine



    subroutine solve_Poisson(qcap_z, dph_z)

        !$ use OMP_LIB
        use Fourier_in_xz_plane
        use mesh, only : n1,n2,n3, n1m,n2m,n3m

        use decomp_2d

        implicit none

        include 'fftw3.f'

        real*8, dimension(:,:,:), intent(in)        :: qcap_z
        real*8, dimension(:,:,:), intent(out)       :: dph_z

        complex*16, dimension(      &
        sp_decomp%xst(1):sp_decomp%xen(1), &
        sp_decomp%xst(2):sp_decomp%xen(2), &
        sp_decomp%xst(3):sp_decomp%xen(3))   :: xz_modes_x

        complex*16, dimension(  &
        sp_decomp%yst(1):sp_decomp%yen(1), &
        sp_decomp%yst(2):sp_decomp%yen(2), &
        sp_decomp%yst(3):sp_decomp%yen(3))   :: xz_modes_y


        call real_to_spectral(qcap_z, xz_modes_x, sp_decomp)

        call transpose_x_to_y(xz_modes_x, xz_modes_y, sp_decomp)

        if ((sp_decomp%yst(1)==1).and.(sp_decomp%yst(3)==1)) then
            xz_modes_y(1,1,1)=dcmplx(0.d0,0.d0)
        endif

        call solve_poisson_in_fourier_space(xz_modes_y)


        call transpose_y_to_x(xz_modes_y, xz_modes_x, sp_decomp)
        call spectral_to_real(xz_modes_x, dph_z, sp_decomp)


        return

    contains



        subroutine solve_poisson_in_fourier_space(qk_y)

            use Lapack_wrappers

            implicit none

            integer n3mh,k,i,j



            real*8, dimension(  &
            sp_decomp%yst(1):sp_decomp%yen(1), &
            sp_decomp%yst(2):sp_decomp%yen(2), &
            sp_decomp%yst(3):sp_decomp%yen(3))   :: re_qk_y, im_qk_y

            complex*16, dimension(  &
            sp_decomp%yst(1):sp_decomp%yen(1), &
            sp_decomp%yst(2):sp_decomp%yen(2), &
            sp_decomp%yst(3):sp_decomp%yen(3))   :: qk_y

            n3mh=n3m/2+1

            do j=1,n2m
                !do k=1,n3mh
                do k=sp_decomp%yst(3), sp_decomp%yen(3)
                    !do i=1,n1m
                    do i=sp_decomp%yst(1), min(n1m, sp_decomp%yen(1))
                        re_qk_y(i,j,k)=dreal(qk_y(i,j,k))
                        im_qk_y(i,j,k)=dimag(qk_y(i,j,k))
                    enddo
                enddo
            enddo

            !********* Zero wave number k=1 **********
            !********* Zero wave number i=1 **********
            if ((sp_decomp%yst(1)==1).and.(sp_decomp%yst(3)==1)) then
                call Lapack_band_solver(lap_matrix_ik(1,1,:,:), kl, ku, n2m, re_qk_y(1,1:n2m,1))
                call Lapack_band_solver(lap_matrix_ik(1,1,:,:), kl, ku, n2m, im_qk_y(1,1:n2m,1))
            end if

            !********* Non zero wave number i ********
            if (sp_decomp%yst(3)==1) then

                !do i=2,n1m
                do i=max(2, sp_decomp%yst(1)),min(n1m, sp_decomp%yen(1))

                    call Lapack_band_solver(lap_matrix_ik(i,1,:,:), kl, ku, n2m, re_qk_y(i,1:n2m,1))
                    call Lapack_band_solver(lap_matrix_ik(i,1,:,:), kl, ku, n2m, im_qk_y(i,1:n2m,1))

                enddo

            end if

            !*********** Non zero wave number k ********
            !do k=2,n3mh
            do k=max(2, sp_decomp%yst(3)), sp_decomp%yen(3)
                !do i=1,n1m
                do i=sp_decomp%yst(1), min(n1m, sp_decomp%yen(1))

                    call Lapack_band_solver(lap_matrix_ik(i,k,:,:), kl, ku, n2m, re_qk_y(i,1:n2m,k))
                    call Lapack_band_solver(lap_matrix_ik(i,k,:,:), kl, ku, n2m, im_qk_y(i,1:n2m,k))
                enddo
            enddo

            do j=1,n2m
                !do k=1,n3mh
                do k=sp_decomp%yst(3), sp_decomp%yen(3)
                    !do i=1,n1m
                    do i=sp_decomp%yst(1), min(n1m, sp_decomp%yen(1))
                        qk_y(i,j,k)=dcmplx(re_qk_y(i,j,k),im_qk_y(i,j,k))
                    enddo
                enddo
            enddo

            return
        end subroutine


    end subroutine

end module poisson_020_solver
