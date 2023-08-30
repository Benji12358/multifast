module VELOCITY_solver

    use schemes3D_interface
    use time_schemes
    use physical_fields
    use boundaries
    use mesh

    use irregular_derivative_coefficients
    use DNS_settings
    use poisson_020_solver, poisson_020_init=>init, ORL_solve_Poisson=>solve_Poisson
    use poisson_interface, LAMB_solve_Poisson=>solve_Poisson
    use lamballais_transfo

    use decomp_2d
    use mpi

    use MHD_data, only:MHD_state
    use IBM, only:IBM_activated, interpol_type

    implicit none

    ! Only update_velocity is callable from the outside
    public :: update_velocity, init, save_state, update_multiphysics_velocity, add_action_x, add_action_z, substract_action
    public :: Fext1, Fext2, Fext3, mean_grad_P1, mean_grad_P3, mean_grad_P1_viscous
    public :: finalize

    private

    procedure(LAMB_solve_Poisson), pointer  :: solve_Poisson

    logical     :: previousRHS_are_available=.false.

    ! Mean pressure gradient, used for both cases, with IBM and without IBM
    real*8              :: mean_grad_P1, mean_grad_P3, mean_grad_P1_viscous

    ! Navier-stokes equation coefficients
    real*8      :: diff21_coef, diff22_coef
    real*8      :: diff1_coef, diff2_coef , diff3_coef
    real*8      :: conv1_coef, conv2_coef, conv3_coef

    ! Navier-stokes equation coefficients for second order box IBM
    real*8      :: diff21_coef_ibm, diff22_coef_ibm
    real*8      :: diff1_coef_ibm, diff2_coef_ibm , diff3_coef_ibm
    real*8      :: conv1_coef_ibm, conv2_coef_ibm, conv3_coef_ibm

    real*8, dimension(:,:,:), allocatable       :: p13_z, p13_y, p13_x
    real*8, dimension(:,:,:), allocatable       :: p12_y, p12_x, p12_z
    real*8, dimension(:,:,:), allocatable       :: p23_z, p23_y, p23_x

    real*8, dimension(:,:,:), allocatable     :: RHS1_x, RHS1_y, RHS1_z
    real*8, dimension(:,:,:), allocatable     :: RHS2_x, RHS2_y, RHS2_z
    real*8, dimension(:,:,:), allocatable     :: RHS3_x, RHS3_y, RHS3_z

    real*8, dimension(:,:,:), allocatable     :: gradP1_x
    real*8, dimension(:,:,:), allocatable     :: gradP2_x, gradP2_y
    real*8, dimension(:,:,:), allocatable     :: gradP3_x, gradP3_y, gradP3_z


    real*8, dimension(:,:,:), allocatable     :: previousRHS1_x
    real*8, dimension(:,:,:), allocatable     :: previousRHS2_x
    real*8, dimension(:,:,:), allocatable     :: previousRHS3_x

    real*8, dimension(:,:,:), allocatable     :: NSTERMS1, NSTERMS2, NSTERMS3
    real*8, dimension(:,:,:), allocatable     :: Fext1, Fext2, Fext3

    real*8, dimension(:,:,:), allocatable     :: q1_save, q2_save, q3_save
    real*8, dimension(:,:,:), allocatable     :: q1_save_x, q2_save_x, q3_save_x
    real*8, dimension(:,:,:), allocatable     :: q1_save_y, q2_save_y, q3_save_y
    real*8, dimension(:,:,:), allocatable     :: q1_save_z, q2_save_z, q3_save_z

    real*8, dimension(:,:), allocatable     :: outflow_conv1o, outflow_conv2o, outflow_conv3o
    real*8, dimension(:,:), allocatable     :: outflow_conv1c, outflow_conv2c, outflow_conv3c
    real*8, dimension(:,:), allocatable     :: outflow_diff1o, outflow_diff2o, outflow_diff3o
    real*8, dimension(:,:), allocatable     :: outflow_diff1c, outflow_diff2c, outflow_diff3c

    real*8, dimension(:,:,:), allocatable     :: RHS1_x_ibm, RHS2_x_ibm, RHS3_x_ibm
    real*8, dimension(:,:,:), allocatable     :: gradP1_x_ibm, gradP2_x_ibm, gradP3_x_ibm

    real*8, dimension(:,:,:), allocatable     :: NSTERMS1_ibm, NSTERMS2_ibm, NSTERMS3_ibm
    real*8, dimension(:,:,:), allocatable     :: Fext1_ibm, Fext2_ibm, Fext3_ibm

    real*8, dimension(:,:,:), allocatable     :: RHS1_y_ibm, RHS2_y_ibm, RHS3_y_ibm
    real*8, dimension(:,:,:), allocatable     :: gradP2_y_ibm, gradP3_y_ibm

    real*8, dimension(:,:,:), allocatable     :: RHS1_z_ibm, RHS2_z_ibm, RHS3_z_ibm
    real*8, dimension(:,:,:), allocatable     :: gradP3_z_ibm

    real*8, dimension(:,:,:), allocatable     :: previousRHS1_x_ibm, previousRHS2_x_ibm, previousRHS3_x_ibm

    real*8, dimension(:,:,:), allocatable     :: p13_x_q1_ibm, p13_x_q3_ibm, p12_x_q1_ibm, p12_x_q2_ibm
    real*8, dimension(:,:,:), allocatable     :: p23_y_q2_ibm, p23_y_q3_ibm, p12_y_q1_ibm, p12_y_q2_ibm
    real*8, dimension(:,:,:), allocatable     :: p23_z_q2_ibm, p23_z_q3_ibm, p13_z_q1_ibm, p13_z_q3_ibm

    real*8, dimension(:,:,:), allocatable     :: q1_save2_ibm, q2_save2_ibm, q3_save2_ibm

    real*8  :: al,gam,rom
    logical :: first_time_advancement=.true.

    integer :: NNN, NN
    integer :: debugproc

contains

    subroutine init()

        use run_ctxt_data
        use numerical_methods_settings
        use poisson_generic_solver
        use VELOCITY_workspace_view, only:recovery_RHS_dir
        use FRINGE_data
        use FRINGE_dao, only:fringe_infos

        use IBM_data
        use IBM

        implicit none

        if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) call activate_O2_for_pressure_correction

        ! Poisson solver initialization
        if (use_generic_poisson) then
            ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) then
            !     call poisson_init(n1_ibm,n2_ibm,n3_ibm, BC1, BC2, BC3, L1, L2, L3, istret, alpha, beta, 1)
            ! else
                call poisson_init(n1,n2,n3, BC1, BC2, BC3, L1, L2, L3, istret, alpha, beta, 1)
            ! endif
            call generic_poisson_infos
            solve_Poisson=>LAMB_solve_Poisson
        else
            call poisson_020_init
            solve_Poisson=>ORL_solve_Poisson
        end if
        ! call generic_poisson_infos

        if (use_fringe) call fringe_infos

        ! Allocate AND initialize the array used by the VELOCITY_solver (RHS, previousRHS etc.)
        call allocate_data

        if ((run_ctxt==CONTINUE_FROM_PREVIOUS_RUN).or.(run_ctxt==RECOVERY_A_RUN)) then

            inquire( file=trim(recovery_RHS_dir)//"/RHS1.h5", exist=previousRHS_are_available)
            if (previousRHS_are_available) then
                call read_previousRHS
            end if
        endif


    contains

        subroutine allocate_data()

            implicit none

            ! RHS arrays allocations
            allocate(RHS1_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            RHS1_x=0.d0
            allocate(RHS2_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            RHS2_x=0.d0
            allocate(RHS3_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            RHS3_x=0.d0
            allocate(gradP1_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            gradP1_x=0.d0
            allocate(gradP2_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            gradP2_x=0.d0
            allocate(gradP3_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            gradP3_x=0.d0

            allocate(NSTERMS1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            NSTERMS1=0.d0
            allocate(NSTERMS2(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            NSTERMS2=0.d0
            allocate(NSTERMS3(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            NSTERMS3=0.d0
            allocate(Fext1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            Fext1=0.d0
            allocate(Fext2(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            Fext2=0.d0
            allocate(Fext3(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            Fext3=0.d0

            allocate(RHS1_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            RHS1_y=0.d0
            allocate(RHS2_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            RHS2_y=0.d0
            allocate(RHS3_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            RHS3_y=0.d0
            allocate(gradP2_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            gradP2_y=0.d0
            allocate(gradP3_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            gradP3_y=0.d0

            allocate(RHS1_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
            RHS1_z=0.d0
            allocate(RHS2_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
            RHS2_z=0.d0
            allocate(RHS3_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
            RHS3_z=0.d0
            allocate(gradP3_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
            gradP3_z=0.d0


            allocate(previousRHS1_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            previousRHS1_x=0.d0
            allocate(previousRHS2_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            previousRHS2_x=0.d0
            allocate(previousRHS3_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            previousRHS3_x=0.d0

            allocate(q1_save_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            q1_save_x=0.d0
            allocate(q2_save_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            q2_save_x=0.d0
            allocate(q3_save_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            q3_save_x=0.d0

            allocate(q1_save_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            q1_save_y=0.d0
            allocate(q2_save_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            q2_save_y=0.d0
            allocate(q3_save_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            q3_save_y=0.d0

            allocate(q1_save_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
            q1_save_z=0.d0
            allocate(q2_save_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
            q2_save_z=0.d0
            allocate(q3_save_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
            q3_save_z=0.d0

            if (streamwise==1) then
                allocate(q1_save(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
                q1_save=0.d0
                allocate(q2_save(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
                q2_save=0.d0
                allocate(q3_save(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
                q3_save=0.d0
            endif

            if (streamwise==3) then
                allocate(q1_save(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
                q1_save=0.d0
                allocate(q2_save(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
                q2_save=0.d0
                allocate(q3_save(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
                q3_save=0.d0
            endif


            ! Cross quantities array allocation
            allocate(p13_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
            p13_z=0.d0
            allocate(p23_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
            p23_z=0.d0
            allocate(p12_z(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
            p12_z=0.d0

            allocate(p23_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            p23_y=0.d0
            allocate(p13_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            p13_y=0.d0
            allocate(p12_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            p12_y=0.d0

            allocate(p13_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            p13_x=0.d0
            allocate(p12_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            p12_x=0.d0
            allocate(p23_x(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            p23_x=0.d0

            allocate(outflow_conv1o(xstart(2):xend(2), xstart(3):xend(3)))
            outflow_conv1o=0.d0
            allocate(outflow_conv2o(xstart(2):xend(2), xstart(3):xend(3)))
            outflow_conv2o=0.d0
            allocate(outflow_conv3o(xstart(2):xend(2), xstart(3):xend(3)))
            outflow_conv3o=0.d0

            allocate(outflow_conv1c(xstart(2):xend(2), xstart(3):xend(3)))
            outflow_conv1c=0.d0
            allocate(outflow_conv2c(xstart(2):xend(2), xstart(3):xend(3)))
            outflow_conv2c=0.d0
            allocate(outflow_conv3c(xstart(2):xend(2), xstart(3):xend(3)))
            outflow_conv3c=0.d0

            allocate(outflow_diff1o(xstart(2):xend(2), xstart(3):xend(3)))
            outflow_diff1o=0.d0
            allocate(outflow_diff2o(xstart(2):xend(2), xstart(3):xend(3)))
            outflow_diff2o=0.d0
            allocate(outflow_diff3o(xstart(2):xend(2), xstart(3):xend(3)))
            outflow_diff3o=0.d0

            allocate(outflow_diff1c(xstart(2):xend(2), xstart(3):xend(3)))
            outflow_diff1c=0.d0
            allocate(outflow_diff2c(xstart(2):xend(2), xstart(3):xend(3)))
            outflow_diff2c=0.d0
            allocate(outflow_diff3c(xstart(2):xend(2), xstart(3):xend(3)))
            outflow_diff3c=0.d0

            if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) then

                ! q1_location

                allocate(q1_save2_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                q1_save2_ibm=0.d0
                allocate(q2_save2_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                q2_save2_ibm=0.d0
                allocate(q3_save2_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                q3_save2_ibm=0.d0

                allocate(RHS1_x_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                RHS1_x_ibm=0.d0
                allocate(gradP1_x_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                gradP1_x_ibm=0.d0
                allocate(NSTERMS1_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                NSTERMS1_ibm=0.d0
                allocate(previousRHS1_x_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                previousRHS1_x_ibm=0.d0
                allocate(Fext1_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                Fext1_ibm=0.d0
                allocate(Fext2_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                Fext2_ibm=0.d0
                allocate(Fext3_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                Fext3_ibm=0.d0

                allocate(RHS1_y_ibm(ystart_fine(1):yend_fine(1), ystart_fine(2):yend_fine(2), ystart_fine(3):yend_fine(3)))
                RHS1_y_ibm=0.d0

                allocate(RHS1_z_ibm(zstart_fine(1):zend_fine(1), zstart_fine(2):zend_fine(2), zstart_fine(3):zend_fine(3)))
                RHS1_z_ibm=0.d0


                ! q2_location

                allocate(RHS2_x_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                RHS2_x_ibm=0.d0
                allocate(gradP2_x_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                gradP2_x_ibm=0.d0
                allocate(NSTERMS2_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                NSTERMS2_ibm=0.d0
                allocate(previousRHS2_x_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                previousRHS2_x_ibm=0.d0

                allocate(RHS2_y_ibm(ystart_fine(1):yend_fine(1), ystart_fine(2):yend_fine(2), ystart_fine(3):yend_fine(3)))
                RHS2_y_ibm=0.d0
                allocate(gradP2_y_ibm(ystart_fine(1):yend_fine(1), ystart_fine(2):yend_fine(2), ystart_fine(3):yend_fine(3)))
                gradP2_y_ibm=0.d0

                allocate(RHS2_z_ibm(zstart_fine(1):zend_fine(1), zstart_fine(2):zend_fine(2), zstart_fine(3):zend_fine(3)))
                RHS2_z_ibm=0.d0


                ! q3_location

                allocate(RHS3_x_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                RHS3_x_ibm=0.d0
                allocate(gradP3_x_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                gradP3_x_ibm=0.d0
                allocate(NSTERMS3_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                NSTERMS3_ibm=0.d0
                allocate(previousRHS3_x_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                previousRHS3_x_ibm=0.d0

                allocate(RHS3_y_ibm(ystart_fine(1):yend_fine(1), ystart_fine(2):yend_fine(2), ystart_fine(3):yend_fine(3)))
                RHS3_y_ibm=0.d0
                allocate(gradP3_y_ibm(ystart_fine(1):yend_fine(1), ystart_fine(2):yend_fine(2), ystart_fine(3):yend_fine(3)))
                gradP3_y_ibm=0.d0

                allocate(RHS3_z_ibm(zstart_fine(1):zend_fine(1), zstart_fine(2):zend_fine(2), zstart_fine(3):zend_fine(3)))
                RHS3_z_ibm=0.d0
                allocate(gradP3_z_ibm(zstart_fine(1):zend_fine(1), zstart_fine(2):zend_fine(2), zstart_fine(3):zend_fine(3)))
                gradP3_z_ibm=0.d0
                
                
                ! for velocity products
                allocate(p13_x_q1_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                p13_x_q1_ibm = 0.d0
                allocate(p13_x_q3_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                p13_x_q3_ibm = 0.d0

                allocate(p12_x_q1_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                p12_x_q1_ibm = 0.d0
                allocate(p12_x_q2_ibm(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)))
                p12_x_q2_ibm = 0.d0


                allocate(p23_y_q2_ibm(ystart_fine(1):yend_fine(1), ystart_fine(2):yend_fine(2), ystart_fine(3):yend_fine(3)))
                p23_y_q2_ibm = 0.d0
                allocate(p23_y_q3_ibm(ystart_fine(1):yend_fine(1), ystart_fine(2):yend_fine(2), ystart_fine(3):yend_fine(3)))
                p23_y_q3_ibm = 0.d0

                allocate(p12_y_q1_ibm(ystart_fine(1):yend_fine(1), ystart_fine(2):yend_fine(2), ystart_fine(3):yend_fine(3)))
                p12_y_q1_ibm = 0.d0
                allocate(p12_y_q2_ibm(ystart_fine(1):yend_fine(1), ystart_fine(2):yend_fine(2), ystart_fine(3):yend_fine(3)))
                p12_y_q2_ibm = 0.d0


                allocate(p23_z_q2_ibm(zstart_fine(1):zend_fine(1), zstart_fine(2):zend_fine(2), zstart_fine(3):zend_fine(3)))
                p23_z_q2_ibm = 0.d0
                allocate(p23_z_q3_ibm(zstart_fine(1):zend_fine(1), zstart_fine(2):zend_fine(2), zstart_fine(3):zend_fine(3)))
                p23_z_q3_ibm = 0.d0

                allocate(p13_z_q1_ibm(zstart_fine(1):zend_fine(1), zstart_fine(2):zend_fine(2), zstart_fine(3):zend_fine(3)))
                p13_z_q1_ibm = 0.d0
                allocate(p13_z_q3_ibm(zstart_fine(1):zend_fine(1), zstart_fine(2):zend_fine(2), zstart_fine(3):zend_fine(3)))
                p13_z_q3_ibm = 0.d0

            endif


        end subroutine allocate_data

        subroutine read_previousRHS()

            use HDF5_IO

            implicit none
            character(200)    :: file_path

            file_path=trim(recovery_RHS_dir)//"/RHS1"
            call hdf_read_3Dfield(file_path, previousRHS1_x, "RHS1", nx_global,ny_global,nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))

            file_path=trim(recovery_RHS_dir)//"/RHS2"
            call hdf_read_3Dfield(file_path, previousRHS2_x, "RHS2", nx_global,ny_global,nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))

            file_path=trim(recovery_RHS_dir)//"/RHS3"
            call hdf_read_3Dfield(file_path, previousRHS3_x, "RHS3", nx_global,ny_global,nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))


        end subroutine read_previousRHS

    end subroutine init


    ! This subroutine transpose all physical quantities. Theses quantities must been known for all sub-domains
    ! orientations.
    ! At entry:
    !   - U, Pr     must be in      Z-configuration
    !   - V         must be in      Y-configuration
    !   - W         must be in      X-configuration


    subroutine update_velocity(ntime, ns)

        use numerical_methods_settings
        use VELOCITY_bc_controller
        use VELOCITY_operations
        use IBM_settings
        use inflow_settings

        use VELOCITY_inout_flow_old, only:    &
        OPEN_OLD_get_inflow=>get_inflow

        use VELOCITY_inout_flow, only:    &
        OPEN_get_inflow=>get_inflow

        use FRINGE_data
        use FRINGE_solver, only:    &
        FRINGE_set_inflow=>set_inflow
    
        use IBM_data
        use IBM
        use mpi
        
        use time_schemes, only:nb_substep

        implicit none

        integer, intent(in) :: ns, ntime
        integer             :: j, k, n1s, s
        integer :: mpi_err

        real*8 :: t_bef, t_aft

        debugproc=-1

        !***********************************************************
        !******************  CORE OF THE DNS ***********************
        !***********************************************************

        ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) call deactivate_O2_for_intermediate_velocity

        !   Time integration implicit viscous
        !   if nb_substep=1 Adams-Bashfort if nb_substep=3 Runge-Kutta
        call set_time_coeffs
        call set_NS_coeffs

        ! ****************  RHS COMPUTATION ************************
        call perform_velocity_products
        t_bef = MPI_WTIME()
        call perform_RHS

        ! write(*,*) 'aft perform_RHS'

        t_aft = MPI_WTIME()
        ! write(*,*) 'perform_RHS', t_aft-t_bef

        ! Define the velocity field at the inlet/outlet for an open boundary problem
        ! MUST BE AFTER the RHS calculation: the inflow and outflow obtained are for the next time step (n+1)
        ! and must not been applied to the current step velocity field
        ! MUST BE BEFORE calculation of u* : the outflow condition is performed from velocity field u at the current
        ! time step
        if (BC1==OPEN) then
            if (.not. inout_newversion) call OPEN_OLD_get_inflow(ntime, ns)
            if (inout_newversion) call OPEN_get_inflow
!            call SCALAR_OPEN_get_inflow(ntime, ns)
            call perform_outflow_velocity()
        endif

        ! ! If MHD used, external forces have to be added for before the computation of the pressure gradient
        if (MHD_state.ne.0) then
            RHS1_x = RHS1_x + Fext1
            RHS2_x = RHS2_x + Fext2
            RHS3_x = RHS3_x + Fext3
        endif

        ! The mean pressure gradient is performed so as to fullfill the flow rate conservation
        if (MHD_state.eq.0) then
            if (streamwise==3)  then
                call perform_mean_gradP3(mean_grad_P3)
                mean_grad_P1=0.d0
            elseif (streamwise==1)  then
                call perform_mean_gradP1(mean_grad_P1)
                mean_grad_P3=0.d0
            endif
        else
            if (streamwise==3) then
                call perform_mean_gradP3_MHD(dP_streamwise)
                call perform_mean_gradP1_MHD(dP_spanwise)
                mean_grad_P1 = dP_spanwise
                mean_grad_P3 = dP_streamwise
            elseif(streamwise==1) then
                call perform_mean_gradP3_MHD(dP_spanwise)
                call perform_mean_gradP1_MHD(dP_streamwise)
                mean_grad_P3 = dP_spanwise
                mean_grad_P1 = dP_streamwise
            endif
        endif

        ! In case the IBM is activated, the pressure gradient is imposed
        if (IBM_activated) then
            mean_grad_P1_viscous = mean_grad_P1
            call correct_gradP(mean_grad_P1)

            ! if fringe is used for sake of periodicity, we treat this region as "constant flowrate"
            ! Thus, we need to compute the viscous terms and the volume integral of Fext_fringe_x in this region and add it to dp_dx
            if (use_fringe) then
                dP_streamwise = 0.d0
                call perform_mean_gradP1_fringe(dP_streamwise)
                if (nrank==0) write(*,*) 'contribution of Fringe to dPdx', dP_streamwise
                mean_grad_P1 = mean_grad_P1 + dP_streamwise
            endif
        endif

        if (nrank==0) write(*,*) 'mean_grad_P1', mean_grad_P1

        RHS1_x =  RHS1_x - mean_grad_P1
        RHS3_x =  RHS3_x - mean_grad_P3

        ! ******************  UPDATE VELOCITY **********************
        ! RHS=(Non linear terms + Diffusive terms)
        ! PreviousRHS = RHS(n-1)
        ! Obtention of û*, v*, w*

        ! If MHD not used, external forces are added after computation of pressure gradient
        if (MHD_state.eq.0) then
            RHS1_x = RHS1_x + Fext1
            RHS2_x = RHS2_x + Fext2
            RHS3_x = RHS3_x + Fext3
        endif

        ! Store velocity at the current time step for support of an external force.
        ! The arrays "q" will be erased by the current estimation of velocity at n+1 (in the coupling process)
        ! So its value is saved in the q*_save arrays
        q1_save_x=q1_x
        q2_save_x=q2_x
        q3_save_x=q3_x

        call perform_NSterms

        ! Note : here the time advancement is restricted to point where the velocity field is not constrained
        ! by the inflow/outflow condition (
        call perform_intermediate_velocity
        
        ! Save the current RHS for the next iteration ...
        previousRHS1_x=RHS1_x
        previousRHS2_x=RHS2_x
        previousRHS3_x=RHS3_x

        if (streamwise==1) then
            q1_save=q1_x
            q2_save=q2_x
            q3_save=q3_x
        endif

        if (streamwise==3) then
            call velocity_transpose_x_to_z(q1_x, q2_x, q3_x)
            q1_save=q1_z
            q2_save=q2_z
            q3_save=q3_z
        endif

        call final_velocity(ntime, ns)

        return

    contains

        subroutine set_time_coeffs()

            if (first_time_advancement) then

                if (previousRHS_are_available) then

                    gam=ga(ns)
                    rom=ro(ns)

                else
                    gam=1.d0
                    rom=0.d0

                end if

            else
                gam=ga(ns)
                rom=ro(ns)

            end if

            first_time_advancement=.false.

            al=gam+rom

        end subroutine set_time_coeffs

        subroutine set_NS_coeffs()

            diff21_coef=dx2*ren

            diff1_coef=dx1*sqrt(ren)
            diff22_coef=dx2*sqrt(ren)
            diff2_coef=dx2*sqrt(ren)
            diff3_coef=dx3*sqrt(ren)

            conv1_coef=-dx1
            conv2_coef=-dx2
            conv3_coef=-dx3

        end subroutine set_NS_coeffs

        ! Flow rate conservation : mean of grad(P) in flow direction is performed so as to ensure conservation
        ! grad(P)=sum(laplacian(u)) in space
        subroutine perform_mean_gradP1(mean_grad_P1)

            use boundary_scheme
            use IBM
            use DRP_IBM
            use FRINGE_data, only: use_fringe, n_interest_region_start, n_fringe_end

            implicit none

            real*8, intent(out) :: mean_grad_P1

            real*8,dimension (n1)                           :: don1
            real*8,dimension (n2)                           :: don2
            real*8,dimension (n3)                           :: don3
            real*8,dimension (xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))           :: tmp1, tmp2_x, tmp3_x
            real*8,dimension (ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))           :: tmp2, tmp3_y
            real*8,dimension (zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))           :: tmp3

            real*8, dimension (:,:,:), allocatable :: tmp1_ibm, tmp2_ibm, tmp3_ibm

            real*8 s1tot, s1tot_glob, total_volume
            integer k,i,j, mpi_err, n1e, n2e, n3e, n1s, n2s, n3s

            s1tot=0.d0
            s1tot_glob=0.d0
            tmp1=0.d0
            tmp2=0.d0
            tmp3=0.d0

            ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) then
            !     allocate(tmp1_ibm(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            !     allocate(tmp2_ibm(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            !     allocate(tmp3_ibm(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

            !     tmp1_ibm=0.d0
            !     tmp2_ibm=0.d0
            !     tmp3_ibm=0.d0

            ! endif

            ! Z ------------------------------------------------------------------------

            n1e=(min(n1m, zend(1))-zstart(1))+1
            n2e=(min(n2m, zend(2))-zstart(2))+1

            call D2c_3Dz(q1_z, tmp3, zsize(1),n1e,zsize(2),n2e,n3, dx3, .true., NS_Q1_BC3)

            ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) then

            !     if (object_in_current_proc_ccm_z) then

            !         n1s=zstart_ibm_ccm(1)
            !         n2s=zstart_ibm_ccm(2)
            !         n3s=zstart_ibm_ccm(3)

            !         n1e=zend_ibm_ccm(1)
            !         n2e=zend_ibm_ccm(2)
            !         n3e=zend_ibm_ccm(3)

            !         call D2c_IBM_O2_3Dz(q1_z(n1s:n1e,n2s:n2e,n3s:n3e), tmp3_ibm(n1s:n1e,n2s:n2e,n3s:n3e),zsize_ibm_ccm(1),zsize_ibm_ccm(2),zsize_ibm_ccm(3),dx3)

            !     endif

            !     if (object_in_current_proc_cpm_z) then

            !         do i=zstart_ibm_cpm(1), zend_ibm_cpm(1)
            !             do j=zstart_ibm_cpm(2), zend_ibm_cpm(2)
            !                 do k=zstart_ibm_cpm(3), zend_ibm_cpm(3)

            !                     tmp3(i,j,k) = tmp3_ibm(i,j,k)

            !                 enddo
            !             enddo
            !         enddo

            !     endif

            ! endif

            ! Y ------------------------------------------------------------------------

            n1e=(min(n1m, yend(1))-ystart(1))+1
            n3e=(min(n3m, yend(3))-ystart(3))+1

            call D1c_MULT_3Dy(q1_y, tmp2, ysize(1),n1e,n2,ysize(3),n3e, dx2, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,1))
            call D2c_MULTACC_3Dy(q1_y, tmp2, ysize(1),n1e,n2,ysize(3),n3e, dx2, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,2))

            ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) then

            !     if (object_in_current_proc_ccm_y) then

            !         n1s=ystart_ibm_ccm(1)
            !         n2s=ystart_ibm_ccm(2)
            !         n3s=ystart_ibm_ccm(3)

            !         n1e=yend_ibm_ccm(1)
            !         n2e=yend_ibm_ccm(2)
            !         n3e=yend_ibm_ccm(3)

            !         call D1c_IBM_O2_MULTACC_3Dy(q1_y(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_ibm(n1s:n1e,n2s:n2e,n3s:n3e),ysize_ibm_ccm(1),ysize_ibm_ccm(2),ysize_ibm_ccm(3),dx2,Yc_to_YcTr_for_D2(:,1))
            !         call D2c_IBM_O2_MULTACC_3Dy(q1_y(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_ibm(n1s:n1e,n2s:n2e,n3s:n3e),ysize_ibm_ccm(1),ysize_ibm_ccm(2),ysize_ibm_ccm(3),dx2,Yc_to_YcTr_for_D2(:,2))

            !     endif

            !     if (object_in_current_proc_cpm_y) then

            !         do i=ystart_ibm_cpm(1), yend_ibm_cpm(1)
            !             do j=ystart_ibm_cpm(2), yend_ibm_cpm(2)
            !                 do k=ystart_ibm_cpm(3), yend_ibm_cpm(3)

            !                     tmp2(i,j,k) = tmp2_ibm(i,j,k)

            !                 enddo
            !             enddo
            !         enddo

            !     endif

            ! endif

            do k=ystart(3), min(n3m, yend(3))   !do k=1,n3m
                do i=ystart(1), min(n1m, yend(1))   !do i=1,n1m

                    tmp2(i,1,k)= ( q1_y(i,2,k)*a3_d + q1_y(i,1,k)*a2_d + q1_wall20(i,k)*a1_d )/dx2**2
                    tmp2(i,n2m,k)= ( q1_y(i,n2m-1,k)*a1_u + q1_y(i,n2m,k) *a2_u + q1_wall21(i,k)*a3_u)/dx2**2

                enddo
            enddo

            ! X ------------------------------------------------------------------------

            n2e=(min(n2m, xend(2))-xstart(2))+1
            n3e=(min(n3m, xend(3))-xstart(3))+1

            call D2c_3Dx(q1_x, tmp1, n1,xsize(2),n2e,xsize(3),n3e, dx1, .false., NS_Q1_BC1)

            ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) then

            !     if (object_in_current_proc_ccm_x) then

            !         n1s=xstart_ibm_ccm(1)
            !         n2s=xstart_ibm_ccm(2)
            !         n3s=xstart_ibm_ccm(3)

            !         n1e=xend_ibm_ccm(1)
            !         n2e=xend_ibm_ccm(2)
            !         n3e=xend_ibm_ccm(3)

            !         call D2c_IBM_O2_3Dx(q1_x(n1s:n1e,n2s:n2e,n3s:n3e), tmp1_ibm(n1s:n1e,n2s:n2e,n3s:n3e),xsize_ibm_ccm(1),xsize_ibm_ccm(2),xsize_ibm_ccm(3),dx1)

            !     endif

            !     if (object_in_current_proc_cpm_x) then

            !         do i=xstart_ibm_cpm(1), xend_ibm_cpm(1)
            !             do j=xstart_ibm_cpm(2), xend_ibm_cpm(2)
            !                 do k=xstart_ibm_cpm(3), xend_ibm_cpm(3)

            !                     tmp1(i,j,k) = tmp1_ibm(i,j,k)

            !                 enddo
            !             enddo
            !         enddo

            !     endif

            ! endif

            call transpose_y_to_x(tmp2, tmp2_x)

            call transpose_z_to_y(tmp3, tmp3_y)
            call transpose_y_to_x(tmp3_y, tmp3_x)

            do k=xstart(3), min(n3m, xend(3))   !do k=1,n3m
                do j=xstart(2), min(n2m, xend(2))   !do j=1,n2m

                    do i=1,n1m
                        s1tot=s1tot+(tmp1(i,j,k)+tmp2_x(i,j,k)+tmp3_x(i,j,k))*dx1*dx3*cell_size_Y(j)
                    enddo
                enddo
            enddo

            if(isnan(s1tot)) then
                write(*,*)'proc C', nrank
                write(*,*)'channel D: s1tot', s1tot
                write(*,*)'channel D: Q1', sum(q1_x), dx1, NS_Q1_BC1
                write(*,*)'channel D: x3', xstart(3), xend(3)
                write(*,*)'channel D: x2', xstart(2), xend(2)

                call exit

            endif

            call MPI_ALLREDUCE (s1tot, s1tot_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

            total_volume = 2.d0*L1*L3

            mean_grad_P1=s1tot_glob/(total_volume*ren)

            return
        end subroutine

        ! Flow rate conservation : mean of grad(P) in flow direction is performed so as to ensure conservation
        ! grad(P)=sum(laplacian(u)) in space
        subroutine correct_gradP(mean_grad_P1)

            use boundary_scheme
            use IBM
            use DRP_IBM
            use FRINGE_data, only: use_fringe
            use DNS_settings, only:dt

            implicit none

            real*8, intent(out) :: mean_grad_P1

            real*8 :: ut1, ut1_glob,total_volume
            integer k,i,j, mpi_err, n1e, n2e, n3e, n1s, n2s, n3s

            ! total_volume = 2.d0*L1*L3
            ! if (use_fringe) total_volume = total_volume/(1+fringe_length)

            ! ut1=0.d0
            ! ut1_glob=0.d0
            ! do k=xstart(3),min(xend(3),n3-1)
            !     do j=xstart(2),min(xend(2),n2-1)
            !         ut1=ut1+q1_x(1,j,k)*(Y(j+1)-Y(j))*dx3
            !     enddo
            ! enddo

            ! call MPI_ALLREDUCE(ut1, ut1_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
            ! ut1_glob=ut1_glob/(L3*L2)

            ! if (nrank.eq.0) then
            !     write(*,*) 'instantaneous flowrate', ut1_glob
            !     write(*,*) 'initial flowrate', initial_flowrate
            !     ! write(*,*) 'added dp/dx', (1/total_volume)*(abs(ut1_glob-initial_flowrate))/(dt)
            ! endif

            ! mean_grad_P1 = mean_grad_P1 - (1/total_volume)*(abs(ut1_glob-initial_flowrate))/(dt)
            mean_grad_P1 = - (ren_tau_ibm / ren )**2

            return
        end subroutine

        subroutine perform_mean_gradP1_MHD(mean_grad_P1)

            use FRINGE_data, only: n_fringe_start, n_fringe_end
            implicit none

            real*8, intent(out) :: mean_grad_P1

            real*8 s1tot,total_volume,prev_volume
            integer k,i,j, mpi_err, n1e, n2e, n3e

            s1tot=0.d0

            ! ********************** WARNING ************************* !
            ! Integral should be effectuated on CELL CENTERED variable !
            ! ************* Here, RHS3 is FACE CENTERED ************** !
            ! ******* But it is ok for periodic condition ************ !
            ! ******************************************************** !

            do k=xstart(3), min(n3m, xend(3))   !do k=1,n3m
                do j=xstart(2), min(n2m, xend(2))   !do j=1,n2m

                    do i=1,n1m
                        s1tot=s1tot + RHS1_x(i,j,k)*dx1*dx3*cell_size_Y(j)
                    enddo

                enddo
            enddo

            call MPI_ALLREDUCE (s1tot, mean_grad_P1, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

            total_volume = 2.d0*L1*L3
            mean_grad_P1=mean_grad_P1/(total_volume)

            return
        end subroutine

        ! Same as perform_mean_gradP1 but for the direction 3 as flow direction

        subroutine perform_mean_gradP3_MHD(mean_grad_P3)

            implicit none

            real*8, intent(out) :: mean_grad_P3

            real*8 s3tot
            integer k,i,j, mpi_err, n1e, n2e, n3e

            s3tot=0.d0

            ! ********************** WARNING ************************* !
            ! Integral should be effectuated on CELL CENTERED variable !
            ! ************* Here, RHS3 is FACE CENTERED ************** !
            ! ******* But it is ok for periodic condition ************ !
            ! ******************************************************** !


            do k=xstart(3), min(n3m, xend(3))   !do k=1,n3m
                do j=xstart(2), min(n2m, xend(2))   !do j=1,n2m

                    do i=1,n1m
                        s3tot=s3tot + RHS3_x(i,j,k)*dx1*dx3*cell_size_Y(j)
                    enddo

                enddo
            enddo

            call MPI_ALLREDUCE (s3tot, mean_grad_P3, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

            mean_grad_P3 = mean_grad_P3/(L1*L2*L3)

            return
        end subroutine


        ! Same as perform_mean_gradP1 but for the direction 3 as flow direction
        subroutine perform_mean_gradP3(mean_grad_P3)
            use boundary_scheme

            implicit none

            real*8, intent(out) :: mean_grad_P3

            real*8,dimension (n1)           :: don1
            real*8,dimension (n2)           :: don2
            real*8,dimension (n3)           :: don3
            real*8,dimension (xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))          :: tmp1
            real*8,dimension (ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))          :: tmp2
            real*8,dimension (zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))          :: tmp3

            real*8 s3tot, s3tot_glob
            integer k,i,j, mpi_err, n1e, n2e, n3e

            s3tot=0.d0

            ! Z ------------------------------------------------------------------------

            n1e=(min(n1m, zend(1))-zstart(1))+1
            n2e=(min(n2m, zend(2))-zstart(2))+1

            ! d²q3_z/d²x
            call D2c_3Dz(q3_z, tmp3, zsize(1),n1e,zsize(2),n2e,n3, dx3, .false., NS_Q3_BC3)

            do j = zstart(2), min(n2m, zend(2))!do j=1,n2m
                do i=zstart(1), min(n1m, zend(1))!do i=1,n1m
                    do k=1,n3m
                        s3tot=s3tot+tmp3(i,j,k)*dx1*dx3*cell_size_Y(j)
                    enddo
                enddo
            enddo

            ! Y ------------------------------------------------------------------------

            n1e=(min(n1m, yend(1))-ystart(1))+1
            n3e=(min(n3m, yend(3))-ystart(3))+1

            call D1c_MULT_3Dy(q3_y, tmp2, ysize(1),n1e,n2,ysize(3),n3e, dx2, .true., NS_Q3_BC2, Yc_to_YcTr_for_D2(:,1))
            call D2c_MULTACC_3Dy(q3_y, tmp2, ysize(1),n1e,n2,ysize(3),n3e, dx2, .true., NS_Q3_BC2, Yc_to_YcTr_for_D2(:,2))

            do k=ystart(3), min(n3m, yend(3))   !do k=1,n3m
                do i=ystart(1), min(n1m, yend(1))   !do i=1,n1m

                    tmp2(i,1,k)= ( q3_y(i,2,k)*a3_d + q3_y(i,1,k)*a2_d + q3_wall20(i,k)*a1_d )/dx2**2
                    tmp2(i,n2m,k)= ( q3_y(i,n2m-1,k)*a1_u + q3_y(i,n2m,k) *a2_u + q3_wall21(i,k)*a3_u)/dx2**2
                    do j = 1, n2m
                        s3tot=s3tot+tmp2(i,j,k)*dx1*dx3*cell_size_Y(j)
                    end do

                enddo
            enddo

            ! X ------------------------------------------------------------------------

            n2e=(min(n2m, xend(2))-xstart(2))+1
            n3e=(min(n3m, xend(3))-xstart(3))+1

            ! d²u/d²z
            call D2c_3Dx(q3_x, tmp1, n1,xsize(2),n2e,xsize(3),n3e, dx1, .true., NS_Q3_BC1)

            do k=xstart(3), min(n3m, xend(3))   !do k=1,n3m
                do j=xstart(2), min(n2m, xend(2))   !do j=1,n2m

                    if ((BC1/=UNBOUNDED).and.(BC1/=FRINGE)) then
                        tmp1(1,j,k)  =( q3_x(2,j,k)       /3.d0 - q3_x(1,j,k)      + q3_wall10(j,k)*2.d0/3.d0 ) /dx1**2
                        tmp1(n1m,j,k)=( q3_x(n1m-1,j,k)   /3.d0 - q3_x(n1m,j,k)    + q3_wall11(j,k)*2.d0/3.d0)  /dx1**2
                    endif

                    do i=1,n1m
                        s3tot=s3tot+tmp1(i,j,k)*dx1*dx3*cell_size_Y(j)
                    enddo

                enddo
            enddo

            call MPI_ALLREDUCE (s3tot, s3tot_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

            mean_grad_P3=s3tot_glob/(2.d0*L1*L3*ren)

            return
        end subroutine

        ! Compute contribution of the fringe method in the dpdx
        subroutine perform_mean_gradP1_fringe(mean_grad_P1)

            use boundary_scheme
            use FRINGE_data, only: use_fringe, n_interest_region_start, n_fringe_end
            use IBM_data, only: 

            implicit none

            real*8, intent(out) :: mean_grad_P1

            real*8,dimension (n1)                           :: don1
            real*8,dimension (n2)                           :: don2
            real*8,dimension (n3)                           :: don3
            real*8,dimension (xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))           :: tmp1, tmp2_x, tmp3_x
            real*8,dimension (ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))           :: tmp2, tmp3_y
            real*8,dimension (zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))           :: tmp3

            real*8, dimension (:,:,:), allocatable :: tmp1_ibm, tmp2_ibm, tmp3_ibm

            real*8 s1tot, s1tot_glob, total_volume, prev_volume
            integer k,i,j, mpi_err, n1e, n2e, n3e, n1s, n2s, n3s

            s1tot=0.d0
            s1tot_glob=0.d0
            tmp1=0.d0
            tmp2=0.d0
            tmp3=0.d0

            ! Z ------------------------------------------------------------------------

            n1e=(min(n1m, zend(1))-zstart(1))+1
            n2e=(min(n2m, zend(2))-zstart(2))+1

            call D2c_3Dz(q1_z, tmp3, zsize(1),n1e,zsize(2),n2e,n3, dx3, .true., NS_Q1_BC3)

            ! Y ------------------------------------------------------------------------

            n1e=(min(n1m, yend(1))-ystart(1))+1
            n3e=(min(n3m, yend(3))-ystart(3))+1

            call D1c_MULT_3Dy(q1_y, tmp2, ysize(1),n1e,n2,ysize(3),n3e, dx2, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,1))
            call D2c_MULTACC_3Dy(q1_y, tmp2, ysize(1),n1e,n2,ysize(3),n3e, dx2, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,2))

            do k=ystart(3), min(n3m, yend(3))   !do k=1,n3m
                do i=ystart(1), min(n1m, yend(1))   !do i=1,n1m

                    tmp2(i,1,k)= ( q1_y(i,2,k)*a3_d + q1_y(i,1,k)*a2_d + q1_wall20(i,k)*a1_d )/dx2**2
                    tmp2(i,n2m,k)= ( q1_y(i,n2m-1,k)*a1_u + q1_y(i,n2m,k) *a2_u + q1_wall21(i,k)*a3_u)/dx2**2

                enddo
            enddo

            ! X ------------------------------------------------------------------------

            n2e=(min(n2m, xend(2))-xstart(2))+1
            n3e=(min(n3m, xend(3))-xstart(3))+1

            call D2c_3Dx(q1_x, tmp1, n1,xsize(2),n2e,xsize(3),n3e, dx1, .false., NS_Q1_BC1)

            call transpose_y_to_x(tmp2, tmp2_x)

            call transpose_z_to_y(tmp3, tmp3_y)
            call transpose_y_to_x(tmp3_y, tmp3_x)

            do k=xstart(3), min(n3m, xend(3))   !do k=1,n3m
                do j=xstart(2), min(n2m, xend(2))   !do j=1,n2m

                    ! only in the fringe region
                    do i=n_fringe_start,n_fringe_end
                        s1tot=s1tot+(tmp1(i,j,k)+tmp2_x(i,j,k)+tmp3_x(i,j,k))*dx1*dx3*cell_size_Y(j)
                    enddo
                enddo
            enddo

            if(isnan(s1tot)) then
                write(*,*)'proc C', nrank
                write(*,*)'channel D: s1tot', s1tot
                write(*,*)'channel D: Q1', sum(q1_x), dx1, NS_Q1_BC1
                write(*,*)'channel D: x3', xstart(3), xend(3)
                write(*,*)'channel D: x2', xstart(2), xend(2)

                call exit

            endif

            call MPI_ALLREDUCE (s1tot, s1tot_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

            total_volume = 2.d0*L1*L3

            mean_grad_P1=s1tot_glob/((total_volume-ibm_volume)*ren)

            return
        end subroutine

        subroutine perform_NSterms()
            implicit none
            integer             :: i, j, k
            integer             :: n1s, n1e, n2s,n2e, n3s,n3e


            ! Transpositions of grad(P) in X configuration for the time advancement
            call transpose_z_to_y(gradP3_z, gradP3_y)
            call transpose_y_to_x(gradP3_y, gradP3_x)
            call transpose_y_to_x(gradP2_y, gradP2_x)

            n1e=min(n1-1, xend(1))
            n2e=min(n2-1, xend(2))
            n3e=min(n3-1, xend(3))

            if ((BC1==UNBOUNDED).or.(BC1==FRINGE)) n1s=xstart(1)
            if ((BC1/=UNBOUNDED).and.(BC1/=FRINGE)) n1s=max(2,xstart(1))
            n2s=xstart(2)
            n3s=xstart(3)
            do k = n3s, n3e
                do j = n2s, n2e
                    do i = n1s, n1e
!                     NSTERMS1(i,j,k)=(RHS1_x(i,j,k)*gam + previousRHS1_x(i,j,k)*rom - al*(gradP1_x(i,j,k) + mean_grad_P1) + al*(It_Factor_1*PreviousFext1(i,j,k) + It_Factor_2*Fext1(i,j,k)) )*dt
                      NSTERMS1(i,j,k)=(RHS1_x(i,j,k)*gam + previousRHS1_x(i,j,k)*rom - al*gradP1_x(i,j,k) )*dt
                    end do
                end do
            end do


            if (BC2==UNBOUNDED) n2s=xstart(2)
            if (BC2/=UNBOUNDED) n2s=max(2,xstart(2))
            ! On ne touche que q2[2..n1-2]
            if (BC1==OPEN) then
                n1s=2
                n1e=n1-2
            endif

            n3s=xstart(3)

            do k = n3s, n3e
                do j = n2s, n2e
                    do i = n1s, n1e
                    !WARNING : Integral of Fext2 MUST be 0
                        NSTERMS2(i,j,k)=(RHS2_x(i,j,k)*gam + previousRHS2_x(i,j,k)*rom - al*gradP2_x(i,j,k) )*dt
                         !NSTERMS2(i,j,k)= (RHS2_x(i,j,k)*gam + previousRHS2_x(i,j,k)*rom + gradP2_x(i,j,k) )*dt
                    end do
                end do
            end do

            if (BC3==UNBOUNDED) n3s=xstart(3)
            if (BC3/=UNBOUNDED) n3s=max(2,xstart(3))

            ! On ne touche que q3[2..n1-2]
            if (BC1==OPEN) then
                n1s=2
                n1e=n1-2
            endif

            n2s=xstart(2)

            do k = n3s, n3e
                do j = n2s, n2e
                    do i = n1s, n1e
                        NSTERMS3(i,j,k)=(RHS3_x(i,j,k)*gam + previousRHS3_x(i,j,k)*rom - al*gradP3_x(i,j,k) )*dt
!                       NSTERMS3(i,j,k)= (RHS3_x(i,j,k)*gam + previousRHS3_x(i,j,k)*rom + gradP3_x(i,j,k) + Fext3(i,j,k)*al)*dt

                    end do
                end do
            end do

        end subroutine perform_NSterms

    end subroutine

    subroutine update_velocity_ibm()

        use numerical_methods_settings
        use VELOCITY_bc_controller
        use VELOCITY_operations
        use IBM_settings
        use inflow_settings
    
        use IBM_data
        use IBM
        use mpi
        
        implicit none

        integer             :: i,j,k
        integer :: mpi_err
        real*8  :: t_bef, t_aft

        !***********************************************************
        !******************  CORE OF THE DNS ***********************
        !***********************************************************

        !   Time integration implicit viscous
        !   if nb_substep=1 Adams-Bashfort if nb_substep=3 Runge-Kutta
        call set_NS_coeffs_ibm
        t_bef = MPI_WTIME()

        ! ****************  RHS COMPUTATION ************************
        call perform_velocity_products_ibm
        ! write(*,*) 'aft perform_velocity_products_ibm'

        t_aft = MPI_WTIME()
        ! write(*,*) 'perform_velocity_products_ibm', t_aft-t_bef

        t_bef = MPI_WTIME()

        call perform_RHS_ibm
        ! write(*,*) 'aft perform_RHS_ibm'

        t_aft = MPI_WTIME()
        ! write(*,*) 'perform_RHS_ibm', t_aft-t_bef

        RHS1_x_ibm =  RHS1_x_ibm - mean_grad_P1
        RHS3_x_ibm =  RHS3_x_ibm - mean_grad_P3

        ! ******************  UPDATE VELOCITY **********************
        ! RHS=(Non linear terms + Diffusive terms)
        ! PreviousRHS = RHS(n-1)
        ! Obtention of û*, v*, w*

        ! Add external force to RHS terms

        if (object_in_current_proc_ccm_x) then

            do i=xstart_ibm_ccm(1), xend_ibm_ccm(1)
                do j=xstart_ibm_ccm(2), xend_ibm_ccm(2)
                    do k=xstart_ibm_ccm(3), xend_ibm_ccm(3)

                        ! 8 fine cells in 1 coarse cell                        
                        Fext1_ibm(2*i-1,2*j-1,2*k-1) = Fext1(i,j,k)
                        Fext1_ibm(min(2*i,n1_ibm),2*j-1,2*k-1) = Fext1(i,j,k)
                        Fext1_ibm(2*i-1,min(2*j,n2_ibm),2*k-1) = Fext1(i,j,k)
                        Fext1_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),2*k-1) = Fext1(i,j,k)
                        Fext1_ibm(2*i-1,2*j-1,min(2*k,n3_ibm)) = Fext1(i,j,k)
                        Fext1_ibm(2*i-1,min(2*j,n2_ibm),min(2*k,n3_ibm)) = Fext1(i,j,k)
                        Fext1_ibm(min(2*i,n1_ibm),2*j-1,min(2*k,n3_ibm)) = Fext1(i,j,k)
                        Fext1_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),min(2*k,n3_ibm)) = Fext1(i,j,k)

                        ! 8 fine cells in 1 coarse cell
                        Fext2_ibm(2*i-1,2*j-1,2*k-1) = Fext2(i,j,k)
                        Fext2_ibm(min(2*i,n1_ibm),2*j-1,2*k-1) = Fext2(i,j,k)
                        Fext2_ibm(2*i-1,min(2*j,n2_ibm),2*k-1) = Fext2(i,j,k)
                        Fext2_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),2*k-1) = Fext2(i,j,k)
                        Fext2_ibm(2*i-1,2*j-1,min(2*k,n3_ibm)) = Fext2(i,j,k)
                        Fext2_ibm(2*i-1,min(2*j,n2_ibm),min(2*k,n3_ibm)) = Fext2(i,j,k)
                        Fext2_ibm(min(2*i,n1_ibm),2*j-1,min(2*k,n3_ibm)) = Fext2(i,j,k)
                        Fext2_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),min(2*k,n3_ibm)) = Fext2(i,j,k)

                        ! 8 fine cells in 1 coarse cell
                        Fext3_ibm(2*i-1,2*j-1,2*k-1) = Fext3(i,j,k)
                        Fext3_ibm(min(2*i,n1_ibm),2*j-1,2*k-1) = Fext3(i,j,k)
                        Fext3_ibm(2*i-1,min(2*j,n2_ibm),2*k-1) = Fext3(i,j,k)
                        Fext3_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),2*k-1) = Fext3(i,j,k)
                        Fext3_ibm(2*i-1,2*j-1,min(2*k,n3_ibm)) = Fext3(i,j,k)
                        Fext3_ibm(2*i-1,min(2*j,n2_ibm),min(2*k,n3_ibm)) = Fext3(i,j,k)
                        Fext3_ibm(min(2*i,n1_ibm),2*j-1,min(2*k,n3_ibm)) = Fext3(i,j,k)
                        Fext3_ibm(min(2*i,n1_ibm),min(2*j,n2_ibm),min(2*k,n3_ibm)) = Fext3(i,j,k)

                    enddo
                enddo
            enddo

        endif
        
        RHS1_x_ibm = Fext1_ibm + RHS1_x_ibm
        RHS2_x_ibm = Fext2_ibm + RHS2_x_ibm
        RHS3_x_ibm = Fext3_ibm + RHS3_x_ibm

        q1_save2_ibm=q1_x_ibm
        q2_save2_ibm=q2_x_ibm
        q3_save2_ibm=q3_x_ibm

        t_bef = MPI_WTIME()
        call perform_NSterms_ibm
        ! write(*,*) 'aft perform_NSterms_ibm'

        t_aft = MPI_WTIME()
        ! write(*,*) 'perform_NSterms_ibm', t_aft-t_bef

        ! Note : here the time advancement is restricted to point where the velocity field is not constrained
        ! by the inflow/outflow condition (
        t_bef = MPI_WTIME()
        call perform_intermediate_velocity_ibm
        ! write(*,*) 'aft perform_intermediate_velocity_ibm'

        t_aft = MPI_WTIME()
        ! write(*,*) 'perform_intermediate_velocity_ibm', t_aft-t_bef

        ! Save the current RHS for the next iteration ...
        previousRHS1_x_ibm=RHS1_x_ibm
        previousRHS2_x_ibm=RHS2_x_ibm
        previousRHS3_x_ibm=RHS3_x_ibm

        return

    contains

        subroutine set_NS_coeffs_ibm()

            diff21_coef_ibm=dx2_ibm*ren

            diff1_coef_ibm=dx1_ibm*sqrt(ren)
            diff22_coef_ibm=dx2_ibm*sqrt(ren)
            diff2_coef_ibm=dx2_ibm*sqrt(ren)
            diff3_coef_ibm=dx3_ibm*sqrt(ren)

            conv1_coef_ibm=-dx1_ibm
            conv2_coef_ibm=-dx2_ibm
            conv3_coef_ibm=-dx3_ibm

        end subroutine set_NS_coeffs_ibm

        ! Adaptation of perform_NSterms, but on the finest grid while the IBM with second order box is used
        subroutine perform_NSterms_ibm()
            implicit none
            integer             :: i, j, k
            integer             :: n1s, n1e, n2s,n2e, n3s,n3e
            
            n1e=min(n1_ibm-1, xend_ibm_fpm(1))
            n2e=min(n2_ibm-1, xend_ibm_fpm(2))
            n3e=min(n3_ibm-1, xend_ibm_fpm(3))

            ! Transpositions of grad(P) in X configuration for the time advancement
            call transpose_z_to_y(gradP3_z_ibm, gradP3_y_ibm, decomp_fine)
            call transpose_y_to_x(gradP3_y_ibm, gradP3_x_ibm, decomp_fine)
            call transpose_y_to_x(gradP2_y_ibm, gradP2_x_ibm, decomp_fine)

            if ((BC1==UNBOUNDED).or.(BC1==FRINGE)) n1s=xstart_ibm_fpm(1)
            if ((BC1/=UNBOUNDED).and.(BC1/=FRINGE)) n1s=max(2,xstart_ibm_fpm(1))
            n2s=xstart_ibm_fpm(2)
            n3s=xstart_ibm_fpm(3)

            if (object_in_current_proc_fpm_x) then

                do k = n3s, n3e
                    do j = n2s, n2e
                        do i = n1s, n1e

                            NSTERMS1_ibm(i,j,k)=(RHS1_x_ibm(i,j,k)*gam + previousRHS1_x_ibm(i,j,k)*rom - al*gradP1_x_ibm(i,j,k) )*dt

                        end do
                    end do
                end do

            endif

            if (BC2==UNBOUNDED) n2s=xstart_ibm_fpm(2)
            if (BC2/=UNBOUNDED) n2s=max(2,xstart_ibm_fpm(2))
            ! On ne touche que q2[2..n1-2]
            if (BC1==OPEN) then
                n1s=max(2,xstart_ibm_fpm(1))
                n1e=min(n1_ibm-2,xstart_ibm_fpm(1))
            endif

            n3s=xstart_ibm_fpm(3)

            if (object_in_current_proc_fpm_x) then

                do k = n3s, n3e
                    do j = n2s, n2e
                        do i = n1s, n1e

                            NSTERMS2_ibm(i,j,k)=(RHS2_x_ibm(i,j,k)*gam + previousRHS2_x_ibm(i,j,k)*rom - al*gradP2_x_ibm(i,j,k) )*dt

                        end do
                    end do
                end do

            endif

            if (BC3==UNBOUNDED) n3s=xstart_ibm_fpm(3)
            if (BC3/=UNBOUNDED) n3s=max(2,xstart_ibm_fpm(3))

            ! On ne touche que q3[2..n1-2]
            if (BC1==OPEN) then
                n1s=max(2,xstart_ibm_fpm(1))
                n1e=min(n1_ibm-2,xstart_ibm_fpm(1))
            endif

            n2s=xstart_ibm_fpm(2)

            if (object_in_current_proc_fpm_x) then

                do k = n3s, n3e
                    do j = n2s, n2e
                        do i = n1s, n1e

                            NSTERMS3_ibm(i,j,k)=(RHS3_x_ibm(i,j,k)*gam + previousRHS3_x_ibm(i,j,k)*rom - al*gradP3_x_ibm(i,j,k) )*dt
                            
                        end do
                    end do
                end do

            endif

        end subroutine perform_NSterms_ibm

    end subroutine

    subroutine update_multiphysics_velocity(ntime, ns)

        implicit none

        integer, intent(in)                                     :: ntime, ns

        call perform_intermediate_velocity
        call final_velocity(ntime, ns)
    !        write (*,*) '------------------- update multi', nrank, sum(abs(previousRHS1_x))

    end subroutine update_multiphysics_velocity


    subroutine final_velocity(ntime, ns)

        use physical_fields
        use DNS_settings
        use VELOCITY_operations
        use numerical_methods_settings
        use VELOCITY_bc_controller
        use IBM_settings
        use IBM_data
        use inflow_settings

        use VELOCITY_inout_flow_old, only:    &
        OPEN_OLD_add_outflow=>add_outflow

        use VELOCITY_inout_flow, only:    &
        OPEN_add_outflow=>add_outflow

        use FRINGE_data
        use FRINGE_solver, only:    &
        FRINGE_set_inflow=>set_inflow

        use IBM_data
        use IBM
        use mpi
        use HDF5_IO

        use mesh, only: Xc, Yc, Zc

        use COMMON_workspace_view, only: COMMON_snapshot_path
        use snapshot_writer

        implicit none

        integer, intent(in)                                     :: ns, ntime

        real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: dpdy_old1
        real*8, dimension(xstart(1):xend(1),xstart(2):xend(2), xstart(3):xend(3)) :: dpdx_old, dpdy_old, dpdz_old
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2)) :: dpdy_old3
        real*8                                                  :: errorsum, errorsum_glob, divy_mean
        integer                                                 :: mpi_err, s
        integer                                                 :: i,j,k
        real*8                                                  :: max_q1, min_q1, max_q2, min_q2, max_q3, min_q3
        real*8                                                  :: max_glob_q1, min_glob_q1, max_glob_q2, min_glob_q2, max_glob_q3, min_glob_q3
        ! real*8                                                  :: flowrate_in, flowrate_in_glob, flowrate_fringe_in, flowrate_fringe_in_glob, flowrate_out, flowrate_out_glob

        real*8 :: t_bef, t_aft

        real*8 :: flowrate_in, flowrate_out, flowrate_fringe_out, flowrate_in_glob, flowrate_out_glob, flowrate_fringe_out_glob

        character(200)                :: snaps_dir, snap_dir
        character(200)    :: file_path, snap_path


        ! ******************  CORRECT VELOCITY *********************
        !  RHS3 becomes divergence(q)/(al*dt) (= laplacien (phi))
        if (BC1==OPEN) call force_flowrate_at_outflow

        ! BC and IBM are applied BEFORE correction step: => The velocity at BC and body must be adapted to be
        ! correct AFTER correction step.
        ! Use of an iterative process

        do s = 1, 1         ! Loop for iterative process

            ! Activation the O2 scheme for the pressure correction
            ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) call activate_O2_for_pressure_correction

            ! t_bef = MPI_WTIME()

            ! ! Copy of the coarsest grid on the finest grid
            ! ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) call from_coarse_to_fine_velocity(q1_x,q2_x,q3_x)
            ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) call from_coarse_to_fine_velocity(q1_save_x,q2_save_x,q3_save_x)
            ! ! write(*,*) 'aft from_coarse_to_fine_velocity'

            ! t_aft = MPI_WTIME()
            ! write(*,*) 'from_coarse_to_fine_velocity', t_aft-t_bef


            ! Correction of the velocity in CloseField regions
            ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) then 
                
            !     t_bef = MPI_WTIME()

            !     call update_velocity_ibm

            !     t_aft = MPI_WTIME()
            !     ! write(*,*) 'update_velocity_ibm', t_aft-t_bef

            !     t_bef = MPI_WTIME()

            !     call transpose_x_to_y(q2_x_ibm, q2_y_ibm, decomp_fine)
            !     call transpose_x_to_y(q3_x_ibm, q3_y_ibm, decomp_fine)
            !     call transpose_y_to_z(q3_y_ibm, q3_z_ibm, decomp_fine)

            !     t_aft = MPI_WTIME()
            !     ! write(*,*) 'apply_bc', t_aft-t_bef

            ! endif

            ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) then
            !     call correct_and_average_velocity_ibm

            !     call transpose_y_to_x(q2_y, q2_x)

            !     call transpose_z_to_y(q3_z, q3_y)
            !     call transpose_y_to_x(q3_y, q3_x)

            ! endif

            ! for IBM check
            dpdx_old = dphidx1_x
            dpdy_old = dphidx2_x
            dpdz_old = dphidx3_x

            if (streamwise==1) dpdy_old1=dphidx2_x(n1-1,:,:)

            ! applying current estimation s of grad(p) to BC and body
            ! Vc=Vc+grad(p)[s] (Vc = desired value of velocity)
            ! after correction the velocity in body and BC will be : Vc+grad(p)[N-1]-grad(p)[N]
            ! where n is the number of iteration and will be near to u
            ! call pressure_gradient_transpose_x_to_z(dphidx1_x, dphidx2_x, dphidx3_x)

            if (BC1==OPEN) call apply_gradp_inoutflow_velocity_1

            ! if (BC1==FRINGE) call apply_gradp_fringe_velocity_1

            ! if ((streamwise==1).and.(IBM_activated)) call apply_gradp_IBM_velocity_1

            call transpose_x_to_y(q2_x, q2_y)
            call transpose_x_to_y(q3_x, q3_y)
            call transpose_y_to_z(q3_y, q3_z)

            ! Applying the boundary condition
            call apply_BC3
            call apply_BC2
            call apply_BC1

            ! use the coarsest grid
            call perform_divergence(divu_z, q3_z, q2_y, q1_x)

            if ((IBM_activated).and.(interpol_type==ANTISYMMETRIC_INTERPOL)) then
                ! See Kim and Moin 2001 (doi: 10.1006/jcph.2001.6778)
                ! The source term is here to compute the divergence in flow field

                call perform_source_term(source_term,mask_kim_q3_z*q3_z,mask_kim_q2_y*q2_y,mask_kim_q1_x*q1_x)

                errorsum = 0.d0

                do i = zstart(1), zend(1)
                    do j = zstart(2), min(zend(2),n2m)
                        do k = zstart(3), zend(3)

                            errorsum = errorsum + source_term(i,j,k)*cell_size_Y(j)*dx1*dx3

                        enddo
                    enddo
                enddo

                errorsum_glob=0.d0
                call MPI_ALLREDUCE(errorsum,errorsum_glob,1,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
                if (nrank==0) write(*,*) 'SUM SOURCE TERM:', s, errorsum_glob/(L1*L2*L3)
                
                call solve_Poisson((divu_z-source_term)/(al*dt), dp_z)
            else
                call solve_Poisson(divu_z/(al*dt), dp_z)
            endif

            ! Compute pressure gradient
            call perform_gradp(dp_z)

        end do

        ! Correct velocity
        call correct_velocity

        call spread_to_all_pencil(q3_z, q2_y, q1_x, dp_z)
    

        ! ******************  PERFORM RESIDUAL DIVERGENCE **********
        call perform_divergence(divu_z, q3_z, q2_y, q1_x, divy_mean)
        if (nrank==0) write(*,*) 'divy_mean', divy_mean
        call transpose_z_to_y(divu_z, divu_y)
        call transpose_y_to_x(divu_y, divu_x)

        if (outflow_buff>0) then
            if(.not. inout_newversion) call OPEN_OLD_add_outflow(ntime, ns)
            if(inout_newversion) call OPEN_add_outflow
        endif

        if (IBM_activated) then
    
            max_q1 = maxval(q1_x)
            min_q1 = minval(q1_x)
            max_q2 = maxval(q2_x)
            min_q2 = minval(q2_x)
            max_q3 = maxval(q3_x)
            min_q3 = minval(q3_x)

            call MPI_ALLREDUCE (max_q1, max_glob_q1, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
            call MPI_ALLREDUCE (max_q2, max_glob_q2, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
            call MPI_ALLREDUCE (max_q3, max_glob_q3, 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)

            call MPI_ALLREDUCE (min_q1, min_glob_q1, 1, MPI_DOUBLE_PRECISION , MPI_MIN , MPI_COMM_WORLD , mpi_err)
            call MPI_ALLREDUCE (min_q2, min_glob_q2, 1, MPI_DOUBLE_PRECISION , MPI_MIN , MPI_COMM_WORLD , mpi_err)
            call MPI_ALLREDUCE (min_q3, min_glob_q3, 1, MPI_DOUBLE_PRECISION , MPI_MIN , MPI_COMM_WORLD , mpi_err)

            if (nrank==0) then
                write(*,*) 'q1 velocity', max_glob_q1, min_glob_q1!, maxloc(q1_x), minloc(q1_x)
                write(*,*) 'q2 velocity', max_glob_q2, min_glob_q2!, maxloc(q2_x), minloc(q2_x)
                write(*,*) 'q3 velocity', max_glob_q3, min_glob_q3!, maxloc(q3_x), minloc(q3_x)
            endif

        endif

        !!!! For debug purpose
        ! if (use_fringe) then

        !     ! FOR DEBUG
        !     flowrate_in = 0.d0
        !     flowrate_out = 0.d0
        !     flowrate_fringe_out = 0.d0

        !     flowrate_in_glob = 0.d0
        !     flowrate_out_glob = 0.d0
        !     flowrate_fringe_out_glob = 0.d0

        !     do k=xstart(3),min(xend(3),n3-1)
        !         do j=xstart(2),min(xend(2),n2-1)

        !             flowrate_in = flowrate_in + q1_x(1,j,k)*(Y(j+1)-Y(j))*dx3

        !             flowrate_out = flowrate_out + q1_x(n1m,j,k)*(Y(j+1)-Y(j))*dx3

        !             flowrate_fringe_out = flowrate_fringe_out + q1_x(n_fringe_region,j,k)*(Y(j+1)-Y(j))*dx3

        !         enddo
        !     enddo

        !     call MPI_ALLREDUCE(flowrate_in, flowrate_in_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
        !     flowrate_in_glob=flowrate_in_glob/(L3*L2)

        !     call MPI_ALLREDUCE(flowrate_out, flowrate_out_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
        !     flowrate_out_glob=flowrate_out_glob/(L3*L2)

        !     call MPI_ALLREDUCE(flowrate_fringe_out, flowrate_fringe_out_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
        !     flowrate_fringe_out_glob=flowrate_fringe_out_glob/(L3*L2)

        !     if (nrank==0) then
        !         write(*,*) 'Reb inflow', ren*flowrate_in_glob
        !         write(*,*) 'Reb out fringe', ren*flowrate_fringe_out_glob
        !         write(*,*) 'Reb outflow', ren*flowrate_out_glob
        !     endif

        ! else

        !     ! FOR DEBUG
        !     flowrate_in = 0.d0
        !     flowrate_out = 0.d0

        !     flowrate_in_glob = 0.d0
        !     flowrate_out_glob = 0.d0

        !     do k=xstart(3),min(xend(3),n3-1)
        !         do j=xstart(2),min(xend(2),n2-1)

        !             flowrate_in = flowrate_in + q1_x(1,j,k)*(Y(j+1)-Y(j))*dx3

        !             flowrate_out = flowrate_out + q1_x(n1m,j,k)*(Y(j+1)-Y(j))*dx3

        !         enddo
        !     enddo

        !     call MPI_ALLREDUCE(flowrate_in, flowrate_in_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
        !     flowrate_in_glob=flowrate_in_glob/(L3*L2)

        !     call MPI_ALLREDUCE(flowrate_out, flowrate_out_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
        !     flowrate_out_glob=flowrate_out_glob/(L3*L2)

        !     if (nrank==0) then
        !         write(*,*) 'Reb inflow', ren*flowrate_in_glob
        !         write(*,*) 'Reb outflow', ren*flowrate_out_glob
        !     endif

        ! endif

        contains

            ! The flow at the outflow is corrected so as to have the same flow rate at inflow and outflow
            ! => for global conservation
            subroutine force_flowrate_at_outflow()
                implicit none
                real*8      :: ut1, ut11
                real*8      :: ut, utt
                integer     :: j, k
                integer     :: mpi_err

                ! Perform the flowrate at the inlet...
                ut1=0.d0
                do k=xstart(3),min(xend(3),n3-1)
                    do j=xstart(2),min(xend(2),n2-1)
                        ut1=ut1+q1_x(1,j,k)*(Y(j+1)-Y(j))*dx3
                    enddo
                enddo

                call MPI_ALLREDUCE(ut1, ut11, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
                ut11=ut11/(L3*L2)

                ! ... and the flowrate at the outlet
                ut=0.d0
                do k=xstart(3),min(xend(3),n3-1)
                    do j=xstart(2),min(xend(2),n2-1)
                        ut=ut+q1_x(n1,j,k)*(Y(j+1)-Y(j))*dx3
                    enddo
                enddo

                call MPI_ALLREDUCE(ut, utt, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
                utt=utt/(L3*L2)

                if (nrank==0) print *,'A: FLOW RATE I/O',ut11,utt

                ! Adjusting the outflow rate to inflow rate
                do k=xstart(3),min(xend(3),n3-1)
                    do j=xstart(2),min(xend(2),n2-1)
                        q1_x(n1,j,k)=q1_x(n1,j,k)*ut11/utt
                    enddo
                enddo

                ! Compare the corrected outflowrate with the inflow rate to ensure they are equals
                ut=0.d0
                do k=xstart(3),min(xend(3),n3-1)
                    do j=xstart(2),min(xend(2),n2-1)
                        ut=ut+q1_x(n1,j,k)*(Y(j+1)-Y(j))*dx3
                    enddo
                enddo
                call MPI_ALLREDUCE(ut, utt, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
                utt=utt/(L3*L2)
                if (nrank==0) print *,'B: FLOW RATE I/O',ut11,utt

            end subroutine force_flowrate_at_outflow


            ! When IBM or OPEN boundary condition are activated, a correct target velocity must be applied BEFORE
            ! the correction by the pressure gradient SO AS to obtain the correct desired velocity field AFTER correction
            ! The current estimation of the target velocity is stored in the "q" array
            subroutine apply_gradp_inoutflow_velocity_1
                implicit none
                integer     :: i, j, k

                do k=xstart(3),min(xend(3),n3-1)
                    do j=xstart(2),min(xend(2),n2-1)

                        ! q1_x(n1,j,k)=q1_x(n1,j,k)+dphidx1_x(n1,j,k)*al*dt
                        ! q1_x(1,j,k)=q1_save(1,j,k)+dphidx1_x(1,j,k)*al*dt

                        q2_x(n1-1,j,k)=q2_save(n1-1,j,k)+dphidx2_x(n1-1,j,k)*al*dt
                        q2_x(1,j,k)=q2_save(1,j,k)+dphidx2_x(1,j,k)*al*dt

                        q3_x(n1-1,j,k)=q3_save(n1-1,j,k)+dphidx3_x(n1-1,j,k)*al*dt
                        q3_x(1,j,k)=q3_save(1,j,k)+dphidx3_x(1,j,k)*al*dt

                    enddo
                enddo

            end subroutine apply_gradp_inoutflow_velocity_1

            ! subroutine apply_gradp_inoutflow_velocity_3
            !     implicit none
            !     integer     :: i, j, k

            !     do k=xstart(3),min(xend(3),n3-1)
            !         do j=xstart(2),min(xend(2),n2-1)

            !             q1_z(i,j,n3-1)=q1_save2(i,j,n3-1)+dphidx1_x(i,j,n3-1)*al*dt
            !             q1_z(i,j,1)=q1_save2(i,j,1)+dphidx1_x(i,j,1)*al*dt

            !             q2_z(i,j,n3-1)=q2_save2(i,j,n3-1)+dphidx2_x(i,j,n3-1)*al*dt
            !             q2_z(i,j,1)=q2_save2(i,j,1)+dphidx2_x(i,j,1)*al*dt

            !             q3_z(i,j,n3-1)=q3_save2(i,j,n3-1)+dphidx3_x(i,j,n3-1)*al*dt
            !             q3_z(i,j,1)=q3_save2(i,j,1)+dphidx3_x(i,j,1)*al*dt

            !         enddo
            !     enddo

            ! end subroutine apply_gradp_inoutflow_velocity_3

            subroutine apply_gradp_fringe_velocity_1
                use FRINGE_data

                implicit none
                integer     :: i, j, k

                do k=xstart(3),min(xend(3),n3-1)
                    do j=xstart(2),min(xend(2),n2-1)

                        ! take the velocity after the inflow is set
                        q1_x(n_interest_region_start,j,k)=q1_save_x(n_interest_region_start,j,k)+dphidx1_x(n_interest_region_start,j,k)*al*dt
                        q2_x(n_interest_region_start,j,k)=q2_save_x(n_interest_region_start,j,k)+dphidx2_x(n_interest_region_start,j,k)*al*dt
                        q3_x(n_interest_region_start,j,k)=q3_save_x(n_interest_region_start,j,k)+dphidx3_x(n_interest_region_start,j,k)*al*dt

                        ! do i=n_fringe_start,n1-1

                        !     q1_x(i,j,k)=q1_save(i,j,k)+dphidx1_x(i,j,k)*al*dt
                        !     q2_x(i,j,k)=q2_save(i,j,k)+dphidx2_x(i,j,k)*al*dt
                        !     q3_x(i,j,k)=q3_save(i,j,k)+dphidx3_x(i,j,k)*al*dt

                        ! enddo

                    enddo
                enddo

            end subroutine apply_gradp_fringe_velocity_1

            subroutine apply_gradp_IBM_velocity_1
                use IBM_settings
                use IBM_data
                implicit none
                integer     :: i, j, k

                do i = 1, n1-1
                    do j = xstart(2),min(xend(2),n2-1)
                        do k = xstart(3),min(xend(3),n3-1)

                            q1_x(i,j,k)=q1_save(i,j,k)+mask_kim_q1_x(i,j,k)*dphidx1_x(i,j,k)*al*dt
                            q2_x(i,j,k)=q2_save(i,j,k)+mask_kim_q2_x(i,j,k)*dphidx2_x(i,j,k)*al*dt
                            q3_x(i,j,k)=q3_save(i,j,k)+mask_kim_q3_x(i,j,k)*dphidx3_x(i,j,k)*al*dt

                        end do
                    enddo
                enddo

            end subroutine apply_gradp_IBM_velocity_1

            subroutine apply_source_term_IBM_velocity_1
                use IBM_settings
                use IBM_data
                implicit none
                integer     :: i, j, k
                real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: q1_IBM
                real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: q2_IBM
                real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: q3_IBM

                ! do i = xstart(1),min(xend(1),n1-1)
                !     do j = xstart(2),min(xend(2),n2-1)
                !         do k = xstart(3),min(xend(3),n3-1)

                !             q1_IBM(i,j,k)=q1_x(i,j,k)*IBM_mask1_x(i,j,k)

                !         end do
                !     enddo
                ! enddo

                ! do i = ystart(1),min(yend(1),n1-1)
                !     do j = ystart(2),min(yend(2),n2-1)
                !         do k = ystart(3),min(yend(3),n3-1)

                !             q2_IBM(i,j,k)=q2_y(i,j,k)*IBM_mask2_y(i,j,k)

                !         end do
                !     enddo
                ! enddo

                ! do i = zstart(1),min(zend(1),n1-1)
                !     do j = zstart(2),min(zend(2),n2-1)
                !         do k = zstart(3),min(zend(3),n3-1)

                !             q3_IBM(i,j,k)=q3_z(i,j,k)*IBM_mask3_z(i,j,k)

                !         end do
                !     enddo
                ! enddo

                ! call perform_divergence(q,IBM_mask3_z*q3_z,IBM_mask2_y*q2_y,IBM_mask1_x*q1_x)
                ! call perform_source_term(q,IBM_mask3_z*q3_z,IBM_mask2_y*q2_y,IBM_mask1_x*q1_x)

            end subroutine apply_source_term_IBM_velocity_1

            ! subroutine apply_gradp_IBM_velocity_1_ibm
            !     use IBM_settings
            !     use IBM_data
            !     implicit none
            !     integer     :: i, j, k

            !     do i = 1, n1_ibm-1
            !         do j = xstart_ibm(2),min(xend_ibm(2),n2_ibm-1)
            !             do k = xstart_ibm(3),min(xend_ibm(3),n3_ibm-1)

            !                 if ((i.ge.i_start_ibm_2nd_f).and.(i.le.i_end_ibm_2nd_f).and.(j.ge.j_start_ibm_2nd_f).and.(j.le.j_end_ibm_2nd_f).and.(k.ge.k_start_ibm_2nd_f).and.(k.le.k_end_ibm_2nd_f)) then

            !                     q1_x_ibm(i,j,k)=q1_x_ibm(i,j,k)+IBM_mask1_x_ibm(i,j,k)*dphidx1_x_ibm(i,j,k)*al*dt
            !                     q2_x_ibm(i,j,k)=q2_x_ibm(i,j,k)+IBM_mask2_x_ibm(i,j,k)*dphidx2_x_ibm(i,j,k)*al*dt
            !                     q3_x_ibm(i,j,k)=q3_x_ibm(i,j,k)+IBM_mask3_x_ibm(i,j,k)*dphidx3_x_ibm(i,j,k)*al*dt

            !                 endif

            !             end do
            !         enddo
            !     enddo

            ! end subroutine apply_gradp_IBM_velocity_1_ibm


            subroutine perform_gradp(dph_z)

                use numerical_methods_settings
                use schemes_interface

                use COMMON_workspace_view, only: COMMON_snapshot_path
                use snapshot_writer
                use numerical_methods_settings, only: schemes_configuration
                use schemes_loader

                implicit none

                character(200)                :: snaps_dir, snap_dir
                character(200)    :: file_path, snap_path

                real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))      :: dph_z
                real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))      :: dph_x
                real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))      :: dph_y
                integer ::k,j,i,n, n1e, n2e, n3e

                ! Gradx(P) -----------------------------------------------------------------
                n1e=zsize(1)
                n2e=zsize(2)

                call D1ssh_3Dz(dph_z, dphidx3_z, zsize(1),n1e,zsize(2),n2e,n3, dx3, .true., POISSON_PR_BC3)

                call transpose_z_to_y(dphidx3_z, dph_y)
                call transpose_y_to_x(dph_y, dphidx3_x)         ! dph_y is used as temporary array to switch from z to x pencil

                ! Grady(P) -----------------------------------------------------------------
                call transpose_z_to_y(dph_z, dph_y)

                n1e=ysize(1)
                n3e=ysize(3)

                if ((use_generic_poisson).or.(schemes_configuration.eq.O2_SCHEMES)) then

                    call D1ssh_MULT_3Dy(dph_y, dphidx2_y, ysize(1),n1e,n2,ysize(3),n3e, dx2, .true., POISSON_PR_BC2, Y_to_YTr_for_D1)

                else    ! ONLY for channel flow. In this case, dirichlet boundary cond. are applied in Y direction

                    do k=ystart(3), min(n3m, yend(3))
                        do i=ystart(1), min(n1m, yend(1))   !do i=1,n1m
                            call D1s_Tamm_MULT(dph_y(i,:,k), dphidx2_y(i,:,k), n2, dx2, .true., Dirichlet, Y_to_YTr_for_D1)
                        enddo
                    enddo

                end if

                call transpose_y_to_x(dphidx2_y, dphidx2_x)

                ! Gradz(P) -----------------------------------------------------------------
                call transpose_y_to_x(dph_y, dph_x)

                n2e=xsize(2)
                n3e=xsize(3)

                call D1ssh_3Dx(dph_x, dphidx1_x, n1,xsize(2),n2e,xsize(3),n3e, dx1, .true., POISSON_PR_BC1)

                return

            end subroutine perform_gradp

            ! ________________________________________________________________________________________
            !  ****************************** subrout correct_velocity  **********************
            !  This subroutine calculate the solenoidal velocity field
            !       q(n+1)=qhat-grad(RHS2)*dt ,  pr=RHS2
            !  A third order runge-kutta or Adam Bashforth schemes can be used.

            !    The pressure is evaluated at the center of the box
            !    at the near boundary cells Newman B.C for RHS2 are assumed
            !    pre=pressure. RHS2=phi

            !    n+1    n
            !    p    = p  + RHS2 + alpha*dt/Reynolds * laplacien(RHS2)     {eq 8.9 p149 ORLANDI}
            ! ________________________________________________________________________________________

            subroutine correct_velocity
                use numerical_methods_settings

                implicit none
                integer ::k,j,i


                ! Z ------------------------------------------------------------------------

                ! U correction
                !do j = 1, n2m
                do j = zstart(2), min(n2m, zend(2))
                    do i=zstart(1), min(n1m, zend(1))
                        q3_z(i,j,:)=q3_z(i,j,:)-dphidx3_z(i,j,:)*al*dt
                    enddo
                enddo


                ! V correction

                !do k=1,n3m
                do k=ystart(3), min(n3m, yend(3))
                    !do i=1,n1m
                    do i=ystart(1), min(n1m, yend(1))
                        q2_y(i,:,k)=q2_y(i,:,k)-dphidx2_y(i,:,k)*al*dt
                    enddo
                enddo

                ! X ------------------------------------------------------------------------

                do k=xstart(3), min(n3m, xend(3))

                    !do j=1,n2m
                    do j=xstart(2), min(n2m, xend(2))
                        q1_x(:,j,k)=q1_x(:,j,k)-dphidx1_x(:,j,k)*al*dt

!                        pr_x(:,j,k)=0.d0!pr_x(:,j,k) + dph_x(:,j,k)/aldt

                    enddo
                enddo

                return
            end subroutine correct_velocity

            subroutine interpol_velocity
                use numerical_methods_settings
                use IBM_settings
                use IBM_data
                use IBM

                implicit none
                integer ::k,j,i
                integer :: n1s, n1e, n2s,n2e, n3s,n3e

                ! Old routine used for refined mesh with IBM

                ! call compute_antisymmetric_velocity(q1_x, q2_x, q3_x, xsize, vel_term1, vel_term2, vel_term3)

                ! if ((BC1==UNBOUNDED).or.(BC1==FRINGE)) n1s=xstart(1)
                ! if ((BC1/=UNBOUNDED).and.(BC1/=FRINGE)) n1s=max(2,xstart(1))
                ! n2s=xstart(2)
                ! n3s=xstart(3)

                ! do k = n3s, n3e
                !     do j = n2s, n2e
                !         do i = n1s, n1e
                !             q1_x(i,j,k) = q1_x(i,j,k) + IBM_mask_bounds1_x(i,j,k)*vel_term1(i,j,k)
                !         end do
                !     end do
                ! end do


                ! if (BC2==UNBOUNDED) n2s=xstart(2)
                ! if (BC2/=UNBOUNDED) n2s=max(2,xstart(2))
                ! ! On ne touche que q2[2..n1-2]
                ! if (BC1==OPEN) then
                !     n1s=2
                !     n1e=n1-2
                ! endif

                ! n3s=xstart(3)

                ! do k = n3s, n3e
                !     do j = n2s, n2e
                !         do i = n1s, n1e
                !             q2_x(i,j,k) = q2_x(i,j,k) + IBM_mask_bounds2_x(i,j,k)*vel_term2(i,j,k)
                !         end do
                !     end do
                ! end do

                ! if ((BC3==UNBOUNDED).or.(BC3==FRINGE)) n3s=xstart(3)
                ! if ((BC3/=UNBOUNDED).and.(BC3/=FRINGE)) n3s=max(2,xstart(3))

                ! ! On ne touche que q3[2..n1-2]
                ! if (BC1==OPEN) then
                !     n1s=2
                !     n1e=n1-2
                ! endif

                ! n2s=xstart(2)

                ! do k = n3s, n3e
                !     do j = n2s, n2e
                !         do i = n1s, n1e
                !             q2_x(i,j,k) = q2_x(i,j,k) + IBM_mask_bounds2_x(i,j,k)*vel_term2(i,j,k)
                !         end do
                !     end do
                ! end do

                ! call transpose_x_to_y(q1_x, q1_y)
                ! call transpose_x_to_y(q2_x, q2_y)
                ! call transpose_x_to_y(q3_x, q3_y)

                ! call transpose_y_to_z(q1_y, q1_z)
                ! call transpose_y_to_z(q2_y, q2_z)
                ! call transpose_y_to_z(q3_y, q3_z)

            end subroutine interpol_velocity

            ! Adaptation of correct_velocity, but on the finest grid while the IBM with second order box is used
            subroutine correct_and_average_velocity_ibm
                use numerical_methods_settings

                use COMMON_workspace_view, only: COMMON_snapshot_path
                use snapshot_writer
                use HDF5_IO
                use IBM_data
                use schemes_interface

                implicit none

                integer k,j,i


                ! q3_z correction ----------------------------------------------------
                if (object_in_current_proc_cpm_z) then

                    do i=zstart_ibm_cpm(1), min(n1,zend_ibm_cpm(1))
                        do j=zstart_ibm_cpm(2), min(n2,zend_ibm_cpm(2))
                            do k=zstart_ibm_cpm(3), min(n3,zend_ibm_cpm(3))

                                    q3_z(i,j,k) = q3_z_ibm(2*i-1,2*j-1,2*k-1)

                            enddo
                        enddo
                    enddo

                endif

                ! q2_y correction ----------------------------------------------------
                if (object_in_current_proc_cpm_y) then

                    do i=ystart_ibm_cpm(1), min(n1,yend_ibm_cpm(1))
                        do j=ystart_ibm_cpm(2), min(n2,yend_ibm_cpm(2))
                            do k=ystart_ibm_cpm(3), min(n3,yend_ibm_cpm(3))

                                    q2_y(i,j,k) = q2_y_ibm(2*i-1,2*j-1,2*k-1)

                            enddo
                        enddo
                    enddo

                endif

                ! write(*,*) 'q2_y', q2_y(22,1:50,32)
        
                ! q1_x correction ----------------------------------------------------
                if (object_in_current_proc_cpm_x) then

                    do i=xstart_ibm_cpm(1), min(n1,xend_ibm_cpm(1))
                        do j=xstart_ibm_cpm(2), min(n2,xend_ibm_cpm(2))
                            do k=xstart_ibm_cpm(3), min(n3,xend_ibm_cpm(3))

                                    q1_x(i,j,k) = q1_x_ibm(2*i-1,2*j-1,2*k-1)

                            enddo
                        enddo
                    enddo

                endif

                return

            end subroutine correct_and_average_velocity_ibm

    end subroutine final_velocity


    subroutine add_action_x(fx1, fx2, fx3)
        implicit none

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)  :: fx1, fx2, fx3

        Fext1 = fx1
        Fext2 = fx2
        Fext3 = fx3

    end subroutine add_action_x

    subroutine add_action_z(fz1, fz2, fz3)
        implicit none

        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in)  :: fz1, fz2, fz3
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))  :: fy1, fy2, fy3
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: fx1, fx2, fx3

        call transpose_z_to_y(fz1,fy1)
        call transpose_z_to_y(fz2,fy2)
        call transpose_z_to_y(fz3,fy3)
        call transpose_y_to_x(fy1,fx1)
        call transpose_y_to_x(fy2,fx2)
        call transpose_y_to_x(fy3,fx3)

        Fext1 = fx1
        Fext2 = fx2
        Fext3 = fx3

    end subroutine add_action_z

    subroutine substract_action(fx1, fx2, fx3)
        implicit none

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)  :: fx1, fx2, fx3

        Fext1 = Fext1 - fx1
        Fext2 = Fext2 - fx2
        Fext3 = Fext3 - fx3

    end subroutine substract_action

    subroutine perform_outflow_velocity()
        use MPI
        use DNS_settings

        use boundaries
        use irregular_derivative_coefficients
        implicit none
        integer, save   :: nb=1
        !integer, parameter  :: outflow_type=6
        real*8          :: q1max, q1min, q1max_g, q1min_g, a1,a2
        integer         :: n2s,n2e, n3s,n3e, j,k, mpi_err
        real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: cx, cx1, cx2, cx3
        real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: conv1, conv2, conv3
        real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: diff1, diff2, diff3

        diff1=0.d0
        diff2=0.d0
        diff3=0.d0

        nb=nb+1

        if (outflow_type==1) then
            cx1=0.6d0*dt/dx1
            cx2=cx1
            cx3=cx1
        endif

        if (outflow_type==2) call perform_global_cx
        if (outflow_type==3) call perform_local_cx
        if (outflow_type==4) call perform_inlet_cx

        if (outflow_type==5) then
            call perform_global_cx
            ! use d²u/dy²
            call perform_diff(0,1,0)
        endif

        if (outflow_type==6) then
            call perform_global_cx
            ! use d²u/dx²
            call perform_diff(1,0,0)
        endif

        if (outflow_type==7) then
            call perform_global_cx
            ! use d²u/dz²
            call perform_diff(0,0,1)
        endif

        if (outflow_type==8) then
            call perform_global_cx
            ! use d²u/dx² + d²u/dz²
            call perform_diff(1,0,1)
        endif

        if (outflow_type==9) then
            call perform_global_cx
            ! use d²u/dx² + d²u/dy² + d²u/dz²
            call perform_diff(1,1,1)
        endif

        ! if (outflow_type==6) then
        !     call perform_local_cx
        !     call perform_diff
        ! endif

        ! if (outflow_type==6) then
        !     call perform_inlet_cx
        !     call perform_diff
        ! endif

        ! cx1 = 10*cx1
        ! cx2 = 10*cx2
        ! cx3 = 10*cx3

        call perform_conv(1,0,0)

        n2s=xstart(2)
        n3s=xstart(3)
        n2e=xend(2)
        n3e=xend(3)

        do k=n3s,n3e
            do j=n2s,n2e
                q1_x(n1,j,k) = q1_x(n1,j,k) + conv1(j,k) + diff1(j,k)*dt
                q2_x(n1-1,j,k) = q2_x(n1-1,j,k) + conv2(j,k) + diff2(j,k)*dt
                q3_x(n1-1,j,k) = q3_x(n1-1,j,k) + conv3(j,k) + diff3(j,k)*dt
                q1_wall11(j,k) = q1_x(n1,j,k)
                q2_wall11(j,k) = q2_x(n1-1,j,k)
                q3_wall11(j,k) = q3_x(n1-1,j,k)
            enddo
        enddo

        if (nrank==debugproc) then
            open(15, file="debugoutflow", position="append")
            write(15,*)"OUTFLOW A", sum(abs(conv1(:,:)))/(xsize(2)*xsize(3)), sum(abs(conv2(:,:)))/(xsize(2)*xsize(3)), sum(abs(conv3(:,:)))/(xsize(2)*xsize(3))
            write(15,*)"OUTFLOW B", sum(abs(diff1(:,:)))/(xsize(2)*xsize(3)), sum(abs(diff2(:,:)))/(xsize(2)*xsize(3)), sum(abs(diff3(:,:)))/(xsize(2)*xsize(3))
            close(15)
        endif

    contains

        subroutine perform_diff(ax,ay,az)
            use boundary_scheme
            use second_order_schemes
            implicit none
            integer                 :: i,j,k
            integer, intent(in)     :: ax,ay,az
            real*8                  :: diff1_coef, diff2_coef, diff3_coef, diff22_coef, diff21_coef
            real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: diff11, diff21, diff31, test
            real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: test_x
            real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: diff12_y, diff22_y, diff32_y
            real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: diff13_z, diff23_z, diff33_z
            real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: diff13_y, diff23_y, diff33_y
            real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: diff12_x, diff22_x, diff32_x, diff13_x, diff23_x, diff33_x
            real*8, dimension(ystart(1):yend(1), ystart(3):yend(3))                    :: tmp

            diff1_coef=dx1*sqrt(ren)
            diff2_coef=dx2*sqrt(ren)
            diff3_coef=dx3*sqrt(ren)
            diff22_coef=dx2*sqrt(ren)
            diff21_coef=dx2*ren

            ! Initializations necessary to avoid VALGRIND errors
            diff11=0.d0
            diff21=0.d0
            diff31=0.d0
            diff13_z=0.d0
            diff23_z=0.d0
            diff33_z=0.d0
            diff12_y=0.d0
            diff22_y=0.d0
            diff32_y=0.d0
            diff13_y=0.d0
            diff23_y=0.d0
            diff33_y=0.d0
            diff12_x=0.d0
            diff22_x=0.d0
            diff32_x=0.d0
            diff13_x=0.d0
            diff23_x=0.d0
            diff33_x=0.d0

            ! for outflow_diff1c
            call D2c_3Dz(q1_z, diff13_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, diff3_coef, .true., NS_Q1_BC3)

            call D1c_MULT_3Dy(q1_y, diff12_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), diff21_coef, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,1))
            call D2c_MULTACC_3Dy(q1_y, diff12_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), diff22_coef, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,2))
            ! call D2c_3Dx(q1_x, diff11, n1, xsize(2),xsize(2),xsize(3),xsize(3), diff1_coef, .false., NS_Q1_BC1)
            ! Not necessary a priori as it does not compute what happens at n1

            do i=ystart(1),yend(1)
                do k=ystart(3),yend(3)

                    if (i.eq.n1) then
            
                        call compute_backward_diff_order_2(q1_y(i,n2m,k)/ren, q1_y(i,n2m-1,k)/ren, q1_y(i,n2m-2,k)/ren, q1_y(i,n2m-3,k)/ren, diff12_y(i,n2m,k), Yc, n2m)
                        call compute_backward_diff_order_2(q1_y(i,n2m-1,k)/ren, q1_y(i,n2m-2,k)/ren, q1_y(i,n2m-3,k)/ren, q1_y(i,n2m-4,k)/ren, diff12_y(i,n2m-1,k), Yc, n2m-1)
                        call compute_upward_diff_order_2(q1_y(i,1,k)/ren, q1_y(i,2,k)/ren, q1_y(i,3,k)/ren, q1_y(i,4,k)/ren, diff12_y(i,1,k), Yc, 1)
                        call compute_upward_diff_order_2(q1_y(i,2,k)/ren, q1_y(i,3,k)/ren, q1_y(i,4,k)/ren, q1_y(i,5,k)/ren, diff12_y(i,2,k), Yc, 2)

                    endif
                enddo
            enddo

            do j=xstart(2),xend(2)
                do k=xstart(3),xend(3)

                    call compute_backward_diff_order_2(q1_x(n1,j,k)/ren, q1_x(n1-1,j,k)/ren, q1_x(n1-2,j,k)/ren, q1_x(n1-3,j,k)/ren, diff11(j,k), Z, n1)

                enddo
            enddo

            ! for outflow_diff2c
            call D2c_3Dz(q2_z, diff23_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, diff3_coef, .true., NS_Q2_BC3)

            call D1c_MULT_3Dy(q2_y, diff22_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), diff21_coef, .false., NS_Q2_BC2, Y_to_YTr_for_D2(:,1))
            call D2c_MULTACC_3Dy(q2_y, diff22_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), diff22_coef, .false., NS_Q2_BC2, Y_to_YTr_for_D2(:,2))

            do i=ystart(1),yend(1)
                do k=ystart(3),yend(3)

                    if (i.eq.n1-1) then

                        call compute_backward_diff_order_2(q2_y(i,n2,k)/ren, q2_y(i,n2-1,k)/ren, q2_y(i,n2-2,k)/ren, q2_y(i,n2-3,k)/ren, diff22_y(i,n2m,k), Y, n2)
                        call compute_backward_diff_order_2(q2_y(i,n2-1,k)/ren, q2_y(i,n2-2,k)/ren, q2_y(i,n2-3,k)/ren, q2_y(i,n2-4,k)/ren, diff22_y(i,n2m-1,k), Y, n2-1)
                        call compute_upward_diff_order_2(q2_y(i,1,k)/ren, q2_y(i,2,k)/ren, q2_y(i,3,k)/ren, q2_y(i,4,k)/ren, diff22_y(i,1,k), Y, 1)
                        call compute_upward_diff_order_2(q2_y(i,2,k)/ren, q2_y(i,3,k)/ren, q2_y(i,4,k)/ren, q2_y(i,5,k)/ren, diff22_y(i,2,k), Y, 2)

                    endif
                enddo
            enddo

            do j=xstart(2),xend(2)
                do k=xstart(3),xend(3)

                    call compute_backward_diff_order_2(q2_x(n1-1,j,k)/ren, q2_x(n1-2,j,k)/ren, q2_x(n1-3,j,k)/ren, q2_x(n1-4,j,k)/ren, diff21(j,k), Zc, n1-1)

                enddo
            enddo

            ! for outflow_diff3c
            call D2c_3Dz(q3_z, diff33_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, diff3_coef, .false., NS_Q3_BC3)

            call D1c_MULT_3Dy(q3_y, diff32_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), diff21_coef, .true., NS_Q3_BC2, Yc_to_YcTr_for_D2(:,1))
            call D2c_MULTACC_3Dy(q3_y, diff32_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), diff22_coef, .true., NS_Q3_BC2, Yc_to_YcTr_for_D2(:,2))

            do i=ystart(1),yend(1)
                do k=ystart(3),yend(3)

                    if (i.eq.n1-1) then
                        call compute_backward_diff_order_2(q3_y(i,n2m,k)/ren, q3_y(i,n2m-1,k)/ren, q3_y(i,n2m-2,k)/ren, q3_y(i,n2m-3,k)/ren, diff32_y(i,n2m,k), Yc, n2m)
                        call compute_backward_diff_order_2(q3_y(i,n2m-1,k)/ren, q3_y(i,n2m-2,k)/ren, q3_y(i,n2m-3,k)/ren, q3_y(i,n2m-4,k)/ren, diff32_y(i,n2m,k), Yc, n2m-1)
                        call compute_upward_diff_order_2(q3_y(i,1,k)/ren, q3_y(i,2,k)/ren, q3_y(i,3,k)/ren, q3_y(i,4,k)/ren, diff32_y(i,1,k), Yc, 1)
                        call compute_upward_diff_order_2(q3_y(i,2,k)/ren, q3_y(i,3,k)/ren, q3_y(i,4,k)/ren, q3_y(i,5,k)/ren, diff32_y(i,2,k), Yc, 2)

                    endif
                enddo
            enddo

            do j=xstart(2),xend(2)
                do k=xstart(3),xend(3)

                    call compute_backward_diff_order_2(q3_x(n1-1,j,k)/ren, q3_x(n1-2,j,k)/ren, q3_x(n1-3,j,k)/ren, q3_x(n1-4,j,k)/ren, diff31(j,k), Zc, n1-1)

                enddo
            enddo

            ! transpose everything in X configuration
            call transpose_z_to_y(diff13_z, diff13_y)
            call transpose_z_to_y(diff23_z, diff23_y)
            call transpose_z_to_y(diff33_z, diff33_y)

            call transpose_y_to_x(diff13_y, diff13_x)
            call transpose_y_to_x(diff23_y, diff23_x)
            call transpose_y_to_x(diff33_y, diff33_x)

            call transpose_y_to_x(diff12_y, diff12_x)
            call transpose_y_to_x(diff22_y, diff22_x)
            call transpose_y_to_x(diff32_y, diff32_x)

            outflow_diff1c = ax*diff11(:,:) + ay*diff12_x(n1,:,:) + az*diff13_x(n1,:,:)
            outflow_diff2c = ax*diff21(:,:) + ay*diff22_x(n1-1,:,:) + az*diff23_x(n1-1,:,:)
            outflow_diff3c = ax*diff31(:,:) + ay*diff32_x(n1-1,:,:) + az*diff33_x(n1-1,:,:)

            diff1 = gam*outflow_diff1c + rom*outflow_diff1o
            diff2 = gam*outflow_diff2c + rom*outflow_diff2o
            diff3 = gam*outflow_diff3c + rom*outflow_diff3o

            outflow_diff1o = outflow_diff1c
            outflow_diff2o = outflow_diff2c
            outflow_diff3o = outflow_diff3c


        end subroutine perform_diff

        ! Update condition at outflow from the current time scheme, according an information propagation with phase velocity cx.
        subroutine perform_conv(ax, ay, az)
            use second_order_schemes
            use mesh
            use COMMON_fieldstools
            implicit none
            real*8              :: cs1, cs2, cs3
            real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: conv11, conv21, conv31
            real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: conv12_y, conv22_y, conv32_y
            real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: conv13_z, conv23_z, conv33_z
            real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: conv13_y, conv23_y, conv33_y
            real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: conv12_x, conv22_x, conv32_x, conv13_x, conv23_x, conv33_x
            integer, intent(in)     :: ax, ay, az

            cs1 = 0.d0
            cs2 = 0.d0
            cs3 = 0.d0

            conv11 = 0.d0
            conv21 = 0.d0
            conv31 = 0.d0

            conv12_y = 0.d0
            conv22_y = 0.d0
            conv32_y = 0.d0
            conv13_y = 0.d0
            conv23_y = 0.d0
            conv33_y = 0.d0

            conv13_z = 0.d0
            conv23_z = 0.d0
            conv33_z = 0.d0

            conv12_x = 0.d0
            conv22_x = 0.d0
            conv32_x = 0.d0
            conv13_x = 0.d0
            conv23_x = 0.d0
            conv33_x = 0.d0

            call perform_checksum_x(q1_x(:,:,:), cs1)
            call perform_checksum_x(q2_x(:,:,:), cs2)
            call perform_checksum_x(q3_x(:,:,:), cs3)

            ! cs1 = 10* ax * cs1 * dt
            ! cs2 = ay * cs2 * dt
            ! cs3 = az * cs3 * dt

            ! for outflow_conv1c
            ! call D1c_3Dz(q1_z, conv13_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, dx3, .true., NS_Q1_BC3)
            ! call D1c_MULT_3Dy(q1_y, conv12_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), dx2, .true., NS_Q1_BC2, Yc_to_YcTr_for_D1)
            ! call D1c_3Dx(q1_x, conv11, n1, xsize(2),xsize(2),xsize(3),xsize(3), dx1, .false., NS_Q1_BC1)
            ! Not necessary a priori as it does not compute what happens at n1

            ! do i=ystart(1),yend(1)
            !     do k=ystart(3),yend(3)

            !         if (i.eq.n1) then

            !             call compute_backward_diff_order_1(q1_y(i,n2m,k), q1_y(i,n2m-1,k), q1_y(i,n2m-2,k), conv12_y(i,n2m,k), Yc, n2m)
            !             call compute_upward_diff_order_1(q1_y(i,1,k), q1_y(i,2,k), q1_y(i,3,k), conv12_y(i,1,k), Yc, 1)

            !         endif
            !     enddo
            ! enddo

            do j=xstart(2),xend(2)
                do k=xstart(3),xend(3)

                    call compute_backward_diff_order_1(q1_x(n1,j,k), q1_x(n1-1,j,k), q1_x(n1-2,j,k), conv11(j,k), Z, n1)

                enddo
            enddo

            ! for outflow_conv2c
            ! call D1c_3Dz(q2_z, conv23_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, dx3, .true., NS_Q2_BC3)
            ! call D1c_MULT_3Dy(q2_y, conv22_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), dx2, .false., NS_Q2_BC2, Y_to_YTr_for_D1)
            ! call D1c_3Dx(q2_x, conv21, n1, xsize(2),xsize(2),xsize(3),xsize(3), dx1, .true., NS_Q2_BC1)
            ! Not necessary a priori as it does not compute what happens at n1-1

            ! do i=ystart(1),yend(1)
            !     do k=ystart(3),yend(3)

            !         if (i.eq.n1-1) then

            !             call compute_backward_diff_order_1(q2_y(i,n2,k), q2_y(i,n2-1,k), q2_y(i,n2-2,k), conv22_y(i,n2m,k), Y, n2)
            !             call compute_upward_diff_order_1(q2_y(i,1,k), q2_y(i,2,k), q2_y(i,3,k), conv22_y(i,1,k), Y, 1)

            !         endif
            !     enddo
            ! enddo

            do j=xstart(2),xend(2)
                do k=xstart(3),xend(3)

                    call compute_backward_diff_order_1(q2_x(n1-1,j,k), q2_x(n1-2,j,k), q2_x(n1-3,j,k), conv21(j,k), Zc, n1-1)

                enddo
            enddo

            ! for outflow_conv3c
            ! call D1c_3Dz(q3_z, conv33_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, dx3, .false., NS_Q3_BC3)
            ! call D1c_MULT_3Dy(q3_y, conv32_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), dx2, .true., NS_Q3_BC2, Yc_to_YcTr_for_D1)
            ! call D1c_3Dx(q3_x, conv31, n1, xsize(2),xsize(2),xsize(3),xsize(3), dx1, .true., NS_Q3_BC1)
            ! Not necessary a priori as it does not compute what happens at n1-1

            ! do i=ystart(1),yend(1)
            !     do k=ystart(3),yend(3)

            !         if (i.eq.n1-1) then

            !             call compute_backward_diff_order_1(q3_y(i,n2m,k), q3_y(i,n2m-1,k), q3_y(i,n2m-2,k), conv32_y(i,n2m,k), Yc, n2m)
            !             call compute_upward_diff_order_1(q3_y(i,1,k), q3_y(i,2,k), q3_y(i,3,k), conv32_y(i,1,k), Yc, 1)

            !         endif
            !     enddo
            ! enddo

            do j=xstart(2),xend(2)
                do k=xstart(3),xend(3)

                    call compute_backward_diff_order_1(q3_x(n1-1,j,k), q3_x(n1-2,j,k), q3_x(n1-3,j,k), conv31(j,k), Zc, n1-1)

                enddo
            enddo

            ! transpose everything in X configuration
            call transpose_z_to_y(conv13_z, conv13_y)
            call transpose_y_to_x(conv13_y, conv13_x)

            call transpose_z_to_y(conv23_z, conv23_y)
            call transpose_y_to_x(conv23_y, conv23_x)

            call transpose_z_to_y(conv33_z, conv33_y)
            call transpose_y_to_x(conv33_y, conv33_x)

            call transpose_y_to_x(conv12_y, conv12_x)

            call transpose_y_to_x(conv22_y, conv22_x)

            call transpose_y_to_x(conv32_y, conv32_x)

            outflow_conv1c(:,:) = - cx1 * conv11(:,:)
            outflow_conv2c(:,:) = - cx1 * conv21(:,:)
            outflow_conv3c(:,:) = - cx1 * conv31(:,:)

            ! Time advancement
            conv1 = gam*outflow_conv1c + rom*outflow_conv1o
            conv2 = gam*outflow_conv2c + rom*outflow_conv2o
            conv3 = gam*outflow_conv3c + rom*outflow_conv3o

            outflow_conv1o = outflow_conv1c
            outflow_conv2o = outflow_conv2c
            outflow_conv3o = outflow_conv3c

        end subroutine perform_conv

        ! Perform phase velocity by averaging the flow between 3:n1-4 (in order to avoid limit effects)
        ! cx2, cx3=cx1
        subroutine perform_global_cx()
            use COMMON_fieldstools

            implicit none
            real*8              :: cs1

            call perform_checksum_x(q1_x(:,:,:), cs1)

            cx1 = cs1*dt/dx1
            cx2 = cx1
            cx3 = cx1

        end subroutine

        ! Phase velocity is defined by the local veloctiy at the outflow
        ! cx2, cx3=cx1
        subroutine perform_local_cx()
            implicit none

            do j=xstart(2),min(xend(2), n2-1)
                do k=xstart(3),min(xend(3), n3-1)

                    cx1(j,k) = q1_x(n1-1,j,k)
                    cx1(j,k) = cx1(j,k)*dt/dx1

                enddo
            enddo

            cx2 = cx1
            cx3 = cx1

        end subroutine

        subroutine perform_inlet_cx()
            implicit none

            do j=xstart(2),min(xend(2), n2-1)
                do k=xstart(3),min(xend(3), n3-1)

                    cx1(j,k) = q1_x(1,j,k)
                    cx1(j,k) = cx1(j,k)*dt/dx1

                enddo
            enddo

            cx2 = cx1
            cx3 = cx1

        end subroutine

    end subroutine

    !GENERIC
    subroutine perform_RHS()

        use boundary_scheme
        use DRP_IBM
        use IBM
        use COMMON_workspace_view
        use snapshot_writer
        use mesh
        use object_drawer
        use object_file_reader

        implicit none

        !real*8,dimension(n3)           :: uu
        real*8,dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))           :: q1q1_x, tmp1_x, tmp2_x, tmp3_x
        real*8,dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))           :: q2q2_y, tmp1_y, tmp2_y, tmp3_y
        real*8,dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))           :: q3q3_z, tmp1_z, tmp2_z, tmp3_z

        real*8, dimension(:,:,:), allocatable :: tmp1_z_ibm, tmp2_z_ibm, tmp3_z_ibm, q3q3_z_ibm
        real*8, dimension(:,:,:), allocatable :: tmp1_y_ibm, tmp2_y_ibm, tmp3_y_ibm, q2q2_y_ibm
        real*8, dimension(:,:,:), allocatable :: tmp1_x_ibm, tmp2_x_ibm, tmp3_x_ibm, q1q1_x_ibm

        integer :: k,i,j,n, mpi_err, n1e, n2e, n3e, n1s, n2s, n3s
        real*8, dimension(:,:), allocatable                 :: vertex
        integer, dimension(:,:), allocatable                :: faces
        integer                                             :: nb_faces, nb_vertex

        real*8                                              :: dx1_bef, dx1_aft, dx3_bef, dx3_aft, dx2_bef, dx2_aft
        integer                                             :: i_bef, i_aft, k_bef, k_aft, j_bef, j_aft
        integer                                             :: i0, i1, j0, j1, k0, k1
        

        ! VALGRIND
        q1q1_x=0.d0
        q2q2_y=0.d0
        q3q3_z=0.d0

        ! _______________________________________________________________________________________________________________________
        ! ***********************************************************************************************************************
        ! ****************************** Derivative along the Z direction *******************************************************
        ! ***********************************************************************************************************************

        n1e= zsize(1)!(min(n1m, zend(1))-zstart(1))+1
        n2e=zsize(2)!(min(n2m, zend(2))-zstart(2))+1

        ! ******** CORR_JO ********* !
        ! REinitialisation des RHS_z !
        ! ************************** !

        RHS1_z=0.d0
        RHS2_z=0.d0
        RHS3_z=0.d0

        tmp1_z = 0.d0
        tmp2_z = 0.d0
        tmp3_z = 0.d0

        ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) then

        !     ! Z
        !     allocate(tmp1_z_ibm(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
        !     allocate(tmp2_z_ibm(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
        !     allocate(tmp3_z_ibm(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
        !     allocate(q3q3_z_ibm(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
            
        !     tmp1_z_ibm = 0.d0
        !     tmp2_z_ibm = 0.d0
        !     tmp3_z_ibm = 0.d0
        !     q3q3_z_ibm = 0.d0

        !     ! Y
        !     allocate(tmp1_y_ibm(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
        !     allocate(tmp2_y_ibm(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
        !     allocate(q2q2_y_ibm(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
        !     allocate(tmp3_y_ibm(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            
        !     tmp1_y_ibm = 0.d0
        !     tmp2_y_ibm = 0.d0
        !     q2q2_y_ibm = 0.d0
        !     tmp3_y_ibm = 0.d0

        !     ! X
        !     allocate(tmp1_x_ibm(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        !     allocate(q1q1_x_ibm(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        !     allocate(tmp2_x_ibm(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        !     allocate(tmp3_x_ibm(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            
        !     tmp1_x_ibm = 0.d0
        !     q1q1_x_ibm = 0.d0
        !     tmp2_x_ibm = 0.d0
        !     tmp3_x_ibm = 0.d0

        ! endif

        ! For RHS1
        call D2c_3Dz(q1_z, tmp1_z, zsize(1),n1e,zsize(2),n2e,n3, diff3_coef, .true., NS_Q1_BC3)
        call D1s_ACC_3Dz(p13_z, tmp1_z, zsize(1),n1e,zsize(2),n2e,n3, conv3_coef, .false., NS_P13_BC3)
        ! *********************************** CORR_JO ******************************* !
        ! La D2c_3Dz sur q2 est shifted --> argument .true. (ancienne valeur : false) !
        ! *************************************************************************** !

        ! For RHS2
        call D2c_3Dz(q2_z, tmp2_z, zsize(1),n1e,zsize(2),n2e,n3, diff3_coef, .true., NS_Q2_BC3)
        call D1s_ACC_3Dz(p23_z, tmp2_z, zsize(1),n1e,zsize(2),n2e,n3, conv3_coef, .false., NS_P23_BC3)

        ! For RHS3
        call D2c_3Dz(q3_z, tmp3_z, zsize(1),n1e,zsize(2),n2e,n3, diff3_coef, .false., NS_Q3_BC3)
        call D0s_3Dz(q3_z, q3q3_z, zsize(1),n1e,zsize(2),n2e,n3, NS_Q3_BC3)
        q3q3_z=q3q3_z**2
        call D1ssh_ACC_3Dz(q3q3_z, tmp3_z, zsize(1),n1e,zsize(2),n2e,n3, conv3_coef, .true., NS_P33_BC3)

        ! ******** CORR_JO ********* !
        ! Suppression de cette ligne !
        ! ************************** !
        ! WARNING: A priori inutile, mais la suppression de ce code change légèrement le résultats
        ! pour 220 (mais pas pour 020) : 2.3x10-11 apres 100 i => POURQUOI ???
!        do j = zstart(2), min(n2m, zend(2))     !do j=1,n2m
!            do i=zstart(1), min(n1m, zend(1))       !do i=1,n1m
!                uu(n3-1)=q3q3_z(i,j,n3-1)
!            enddo
!        enddo

        ! If the IBM is activated then the DRP scheme has to be corrected ...
        ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) then

        !     if (object_in_current_proc_ccm_z) then

        !         n1s=zstart_ibm_ccm(1)
        !         n2s=zstart_ibm_ccm(2)
        !         n3s=zstart_ibm_ccm(3)

        !         n1e=zend_ibm_ccm(1)
        !         n2e=zend_ibm_ccm(2)
        !         n3e=zend_ibm_ccm(3)

        !         ! For RHS1
        !         call D2c_IBM_O2_3Dz(q1_z(n1s:n1e,n2s:n2e,n3s:n3e), tmp1_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e),zsize_ibm_ccm(1),zsize_ibm_ccm(2),zsize_ibm_ccm(3),diff3_coef)
        !         call D1s_IBM_O2_ACC_3Dz(p13_z(n1s:n1e,n2s:n2e,n3s:n3e),tmp1_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e),zsize_ibm_ccm(1),zsize_ibm_ccm(2),zsize_ibm_ccm(3),conv3_coef)

        !         ! For RHS2
        !         call D2c_IBM_O2_3Dz(q2_z(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e),zsize_ibm_ccm(1),zsize_ibm_ccm(2),zsize_ibm_ccm(3),diff3_coef)
        !         call D1s_IBM_O2_ACC_3Dz(p23_z(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e),zsize_ibm_ccm(1),zsize_ibm_ccm(2),zsize_ibm_ccm(3),conv3_coef)

        !         ! For RHS1
        !         call D2c_IBM_O2_3Dz(q3_z(n1s:n1e,n2s:n2e,n3s:n3e),tmp3_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e),zsize_ibm_ccm(1),zsize_ibm_ccm(2),zsize_ibm_ccm(3),diff3_coef)

        !         call D0s_IBM_O2_3Dz(q3_z(n1s:n1e,n2s:n2e,n3s:n3e),q3q3_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e),zsize_ibm_ccm(1),zsize_ibm_ccm(2),zsize_ibm_ccm(3))
        !         q3q3_z_ibm=q3q3_z_ibm**2

        !         call D1ssh_IBM_O2_ACC_3Dz(q3q3_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp3_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e),zsize_ibm_ccm(1),zsize_ibm_ccm(2),zsize_ibm_ccm(3),conv3_coef)

        !     endif

        !     if (object_in_current_proc_cpm_z) then

        !         do i=zstart_ibm_cpm(1), zend_ibm_cpm(1)
        !             do j=zstart_ibm_cpm(2), zend_ibm_cpm(2)
        !                 do k=zstart_ibm_cpm(3), zend_ibm_cpm(3)

        !                     tmp1_z(i,j,k) = tmp1_z_ibm(i,j,k)
        !                     tmp2_z(i,j,k) = tmp2_z_ibm(i,j,k)
        !                     tmp3_z(i,j,k) = tmp3_z_ibm(i,j,k)

        !                 enddo
        !             enddo
        !         enddo

        !     endif

        ! endif


        ! if ((IBM_activated).and.(interpol_type == ANTISYMMETRIC_INTERPOL)) then

        !     do n=1,number_of_objects

        !         ! READING OBJECT FILE
        !         call read_object_size(trim(obj_file_path), nb_vertex, nb_faces)
        !         call read_object(trim(obj_file_path), vertex, faces, nb_vertex, nb_faces)

        !         ! Scale with H
        !         vertex(:,1)=vertex(:,1)*body_scale_x1(n)+body_x1(n)
        !         vertex(:,2)=vertex(:,2)*body_scale_x2(n)+body_x2(n)
        !         vertex(:,3)=vertex(:,3)*body_scale_x3(n)+body_x3(n)

        !         call get_objet_area(vertex, nb_vertex)

        !         xs = max(0.d0,xs)
        !         xe = min(L1,xe)
        !         ys = max(0.d0,ys)
        !         ye = min(L2,ye)
        !         zs = max(0.d0,zs)
        !         ze = min(L3,ze)

        !         if (object_in_current_proc_z_q3(n)) then

        !             i0 = i_start_obj_q3(n)
        !             i1 = i_end_obj_q3(n)
        !             j0 = j_start_obj_q3(n)
        !             j1 = j_end_obj_q3(n)
        !             k0 = k_start_obj_q3(n)
        !             k1 = k_end_obj_q3(n)

        !             k_bef = k0-1
        !             k_aft = k1+1

        !             if (k0.eq.1) k_bef=n3-1
        !             if (k1.gt.(n3-1)) k_aft=1

        !             dx3_bef = (zs-X(k0)-dx3)
        !             dx3_aft = (X(k1)+dx3-ze)

        !             do j=zstart_ibm_q3(n,2), zend_ibm_q3(n,2)
        !                 do i=zstart_ibm_q3(n,1), zend_ibm_q3(n,1)

        !                     if ((j.ge.j0).and.(j.le.j1).and.(i.ge.i0).and.(i.le.i1)) then

        !                         ! correct dq3q3/dz term only for k_aft (no need to do this for k_bef)
        !                         ! remove what is computed by generic O2 scheme
        !                         tmp3_z(i,j,k_aft) = tmp3_z(i,j,k_aft) - ( (1.d0/conv3_coef)*(q3q3_z(i,j,k_aft) - q3q3_z(i,j,min(k1,n3m))) )
        !                         ! add what is computed by taking into account the object position
        !                         tmp3_z(i,j,k_aft) = tmp3_z(i,j,k_aft) + (1.d0/(-dx3_aft))*(q3q3_z(i,j,k_aft) - 0.d0)

        !                         ! correct d²q3/dz² term (no need to do this for i_bef)
        !                         ! remove what is computed by generic O2 scheme
        !                         tmp3_z(i,j,k_bef) = tmp3_z(i,j,k_bef) - ( (1.d0/diff3_coef**2)*(q3_z(i,j,k_bef-1) -2.d0*q3_z(i,j,k_bef) + q3_z(i,j,k0)) )
        !                         tmp3_z(i,j,k_aft) = tmp3_z(i,j,k_aft) - ( (1.d0/diff3_coef**2)*(q3_z(i,j,k_aft+1) -2.d0*q3_z(i,j,k_aft) + q3_z(i,j,min(k1,n3m))) )
        !                         ! add what is computed by taking into account the object position
        !                         ! Caution, here we use the interpolated velocity to use the same scheme as O2 (with equal spacing)
        !                         tmp3_z(i,j,k_bef) = tmp3_z(i,j,k_bef) + (1.d0/(dx3_bef*sqrt(ren))**2)*(0.d0 + q3_z(i,j,k_bef-1) + ((dx3-dx3_bef)/dx3)*(q3_z(i,j,k_bef)-q3_z(i,j,k_bef-1))) - 2.d0*(1.d0/(dx3_bef*sqrt(ren))**2)*q3_z(i,j,k_bef)
        !                         tmp3_z(i,j,k_aft) = tmp3_z(i,j,k_aft) + (1.d0/(dx3_aft*sqrt(ren))**2)*(0.d0 + q3_z(i,j,k_aft) + (dx3_bef/dx3)*(q3_z(i,j,k_aft+1)-q3_z(i,j,k_aft))) - 2.d0*(1.d0/(dx3_aft*sqrt(ren))**2)*q3_z(i,j,k_aft)

        !                     endif

        !                 enddo
        !             enddo

        !         endif

        !         if (object_in_current_proc_z_q2(n)) then

        !             i0 = i_start_obj_q2(n)
        !             i1 = i_end_obj_q2(n)
        !             j0 = j_start_obj_q2(n)
        !             j1 = j_end_obj_q2(n)
        !             k0 = k_start_obj_q2(n)
        !             k1 = k_end_obj_q2(n)

        !             k_bef = k0-1
        !             k_aft = k1+1

        !             if (k0.eq.1) k_bef=n3-1
        !             if (k1.gt.(n3-1)) k_aft=1

        !             dx3_bef = (zs-X(k0)-dx3)
        !             dx3_aft = (X(k1)+dx3-ze)

        !             do j=zstart_ibm_q2(n,2), zend_ibm_q2(n,2)
        !                 do i=zstart_ibm_q2(n,1), zend_ibm_q2(n,1)

        !                     if ((j.ge.j0).and.(j.le.j1).and.(i.ge.i0).and.(i.le.i1)) then

        !                         ! correct dq3q3/dz term only for k_aft (no need to do this for k_bef)
        !                         ! remove what is computed by generic O2 scheme
        !                         tmp2_z(i,j,k_bef) = tmp2_z(i,j,k_bef) - ( (1.d0/conv3_coef)*(p23_z(i,j,k0) - p23_z(i,j,k_bef)) )
        !                         ! add what is computed by taking into account the object position
        !                         tmp2_z(i,j,k_bef) = tmp2_z(i,j,k_bef) + (1.d0/(-dx3_aft))*(0.d0 - p23_z(i,j,k_aft))

        !                         ! correct d²q3/dz² term (no need to do this for i_bef)
        !                         ! remove what is computed by generic O2 scheme
        !                         tmp2_z(i,j,k_bef) = tmp2_z(i,j,k_bef) - ( (1.d0/diff3_coef**2)*(q2_z(i,j,k_bef-1) -2.d0*q2_z(i,j,k_bef) + q2_z(i,j,k0)) )
        !                         tmp2_z(i,j,k_aft) = tmp2_z(i,j,k_aft) - ( (1.d0/diff3_coef**2)*(q2_z(i,j,k_aft+1) -2.d0*q2_z(i,j,k_aft) + q2_z(i,j,min(k1,n3m))) )
        !                         ! add what is computed by taking into account the object position
        !                         ! Caution, here we use the interpolated velocity to use the same scheme as O2 (with equal spacing)
        !                         tmp2_z(i,j,k_bef) = tmp2_z(i,j,k_bef) + (1.d0/(dx3_bef*sqrt(ren))**2)*(0.d0 + q2_z(i,j,k_bef-1) + ((dx3-dx3_bef)/dx3)*(q2_z(i,j,k_bef)-q2_z(i,j,k_bef-1))) - 2.d0*(1.d0/(dx3_bef*sqrt(ren))**2)*q2_z(i,j,k_bef)
        !                         tmp2_z(i,j,k_aft) = tmp2_z(i,j,k_aft) + (1.d0/(dx3_aft*sqrt(ren))**2)*(0.d0 + q2_z(i,j,k_aft) + (dx3_bef/dx3)*(q2_z(i,j,k_aft+1)-q2_z(i,j,k_aft))) - 2.d0*(1.d0/(dx3_aft*sqrt(ren))**2)*q2_z(i,j,k_aft)

        !                     endif

        !                 enddo
        !             enddo

        !         endif

        !         if (object_in_current_proc_z_q1(n)) then

        !             i0 = i_start_obj_q1(n)
        !             i1 = i_end_obj_q1(n)
        !             j0 = j_start_obj_q1(n)
        !             j1 = j_end_obj_q1(n)
        !             k0 = k_start_obj_q1(n)
        !             k1 = k_end_obj_q1(n)

        !             k_bef = k0-1
        !             k_aft = k1+1

        !             if (k0.eq.1) k_bef=n3-1
        !             if (k1.gt.(n3-1)) k_aft=1

        !             dx3_bef = (zs-X(k0)-dx3)
        !             dx3_aft = (X(k1)+dx3-ze)

        !             do j=zstart_ibm_q1(n,2), zend_ibm_q1(n,2)
        !                 do i=zstart_ibm_q1(n,1), zend_ibm_q1(n,1)

        !                     if ((j.ge.j0).and.(j.le.j1).and.(i.ge.i0).and.(i.le.i1)) then

        !                         ! correct dq3q3/dz term only for k_aft (no need to do this for k_bef)
        !                         ! remove what is computed by generic O2 scheme
        !                         tmp1_z(i,j,k_bef) = tmp1_z(i,j,k_bef) - ( (1.d0/conv3_coef)*(p13_z(i,j,k0) - p13_z(i,j,k_bef)) )
        !                         ! add what is computed by taking into account the object position
        !                         tmp1_z(i,j,k_bef) = tmp1_z(i,j,k_bef) + (1.d0/(-dx3_aft))*(0.d0 - p13_z(i,j,k_aft))

        !                         ! correct d²q3/dz² term (no need to do this for i_bef)
        !                         ! remove what is computed by generic O2 scheme
        !                         tmp1_z(i,j,k_bef) = tmp1_z(i,j,k_bef) - ( (1.d0/diff3_coef**2)*(q1_z(i,j,k_bef-1) -2.d0*q1_z(i,j,k_bef) + q1_z(i,j,k0)) )
        !                         tmp1_z(i,j,k_aft) = tmp1_z(i,j,k_aft) - ( (1.d0/diff3_coef**2)*(q1_z(i,j,k_aft+1) -2.d0*q1_z(i,j,k_aft) + q1_z(i,j,min(k1,n3m))) )
        !                         ! add what is computed by taking into account the object position
        !                         ! Caution, here we use the interpolated velocity to use the same scheme as O2 (with equal spacing)
        !                         tmp1_z(i,j,k_bef) = tmp1_z(i,j,k_bef) + (1.d0/(dx3_bef*sqrt(ren))**2)*(0.d0 + q1_z(i,j,k_bef-1) + ((dx3-dx3_bef)/dx3)*(q1_z(i,j,k_bef)-q1_z(i,j,k_bef-1))) - 2.d0*(1.d0/(dx3_bef*sqrt(ren))**2)*q1_z(i,j,k_bef)
        !                         tmp1_z(i,j,k_aft) = tmp1_z(i,j,k_aft) + (1.d0/(dx3_aft*sqrt(ren))**2)*(0.d0 + q1_z(i,j,k_aft) + (dx3_bef/dx3)*(q1_z(i,j,k_aft+1)-q1_z(i,j,k_aft))) - 2.d0*(1.d0/(dx3_aft*sqrt(ren))**2)*q1_z(i,j,k_aft)

        !                     endif

        !                 enddo
        !             enddo

        !         endif

        !     enddo

        ! endif



        if (BC3==NOSLIP) call diff_at_wall_3   ! ATTENTION

        RHS1_z = RHS1_z + tmp1_z
        RHS2_z = RHS2_z + tmp2_z
        RHS3_z = RHS3_z + tmp3_z

        ! _______________________________________________________________________________________________________________________
        ! ***********************************************************************************************************************
        ! ****************************** Derivative along the Y direction *******************************************************
        ! ***********************************************************************************************************************

        call transpose_z_to_y(p23_z, p23_y)

        call transpose_z_to_y(RHS1_z, RHS1_y)
        call transpose_z_to_y(RHS2_z, RHS2_y)
        call transpose_z_to_y(RHS3_z, RHS3_y)

        tmp1_y = 0.d0
        tmp2_y = 0.d0
        tmp3_y = 0.d0

        n1e=ysize(1)!(min(n1m, yend(1))-ystart(1))+1
        n3e=ysize(3)!(min(n3m, yend(3))-ystart(3))+1

        ! For RHS1
        call D1c_MULTACC_3Dy(q1_y, tmp1_y, ysize(1),n1e,n2,ysize(3),n3e, diff21_coef, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,1))
        call D2c_MULTACC_3Dy(q1_y, tmp1_y, ysize(1),n1e,n2,ysize(3),n3e, diff22_coef, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,2))

        ! For RHS2
        call D1c_MULTACC_3Dy(q2_y, tmp2_y, ysize(1),n1e,n2,ysize(3),n3e, diff21_coef, .false., NS_Q2_BC2, Y_to_YTr_for_D2(:,1))
        call D2c_MULTACC_3Dy(q2_y, tmp2_y, ysize(1),n1e,n2,ysize(3),n3e, diff22_coef, .false., NS_Q2_BC2, Y_to_YTr_for_D2(:,2))

        ! For RHS3
        call D1c_MULTACC_3Dy(q3_y, tmp3_y, ysize(1),n1e,n2,ysize(3),n3e, diff21_coef, .true., NS_Q3_BC2, Yc_to_YcTr_for_D2(:,1))
        call D2c_MULTACC_3Dy(q3_y, tmp3_y, ysize(1),n1e,n2,ysize(3),n3e, diff22_coef, .true., NS_Q3_BC2, Yc_to_YcTr_for_D2(:,2))

        call D1s_MULTACC_3Dy(p12_y, tmp1_y, ysize(1),n1e,n2,ysize(3),n3e, conv2_coef, .false., NS_P12_BC2, Yc_to_YcTr_for_D1)
        call D0s_3Dy(q2_y, q2q2_y, ysize(1),n1e,n2,ysize(3),n3e, NS_Q2_BC2)
        q2q2_y=q2q2_y**2
        call D1ssh_MULTACC_3Dy(q2q2_y, tmp2_y, ysize(1),n1e,n2,ysize(3),n3e, conv2_coef, .true., NS_P22_BC2, Y_to_YTr_for_D1)

        call D1s_MULTACC_3Dy(p23_y, tmp3_y, ysize(1),n1e,n2,ysize(3),n3e, conv2_coef, .false., NS_P23_BC2, Yc_to_YcTr_for_D1)

        ! If the IBM is activated then the DRP scheme has to be corrected ...
        ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) then

        !     if (object_in_current_proc_ccm_y) then

        !         n1s=ystart_ibm_ccm(1)
        !         n2s=ystart_ibm_ccm(2)
        !         n3s=ystart_ibm_ccm(3)

        !         n1e=yend_ibm_ccm(1)
        !         n2e=yend_ibm_ccm(2)
        !         n3e=yend_ibm_ccm(3)

        !         ! For RHS1
        !         call D1c_IBM_O2_MULTACC_3Dy(q1_y(n1s:n1e,n2s:n2e,n3s:n3e),tmp1_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),ysize_ibm_ccm(1),ysize_ibm_ccm(2),ysize_ibm_ccm(3),diff21_coef,Yc_to_YcTr_for_D2(:,1))
        !         call D2c_IBM_O2_MULTACC_3Dy(q1_y(n1s:n1e,n2s:n2e,n3s:n3e),tmp1_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),ysize_ibm_ccm(1),ysize_ibm_ccm(2),ysize_ibm_ccm(3),diff22_coef,Yc_to_YcTr_for_D2(:,2))

        !         call D1s_IBM_O2_MULTACC_3Dy(p12_y(n1s:n1e,n2s:n2e,n3s:n3e),tmp1_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),ysize_ibm_ccm(1),ysize_ibm_ccm(2),ysize_ibm_ccm(3),conv2_coef,Yc_to_YcTr_for_D1)

        !         ! For RHS2
        !         call D1c_IBM_O2_MULTACC_3Dy(q2_y(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),ysize_ibm_ccm(1),ysize_ibm_ccm(2),ysize_ibm_ccm(3),diff21_coef,Y_to_YTr_for_D2(:,1))
        !         call D2c_IBM_O2_MULTACC_3Dy(q2_y(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),ysize_ibm_ccm(1),ysize_ibm_ccm(2),ysize_ibm_ccm(3),diff22_coef,Y_to_YTr_for_D2(:,2))

        !         call D0s_IBM_O2_3Dy(q2_y(n1s:n1e,n2s:n2e,n3s:n3e),q2q2_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),ysize_ibm_ccm(1),ysize_ibm_ccm(2),ysize_ibm_ccm(3))
        !         q2q2_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e)=q2q2_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e)**2

        !         call D1ssh_IBM_O2_MULTACC_3Dy(q2q2_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),ysize_ibm_ccm(1),ysize_ibm_ccm(2),ysize_ibm_ccm(3),conv2_coef,Y_to_YTr_for_D1)

        !         ! For RHS3
        !         call D1c_IBM_O2_MULTACC_3Dy(q3_y(n1s:n1e,n2s:n2e,n3s:n3e),tmp3_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),ysize_ibm_ccm(1),ysize_ibm_ccm(2),ysize_ibm_ccm(3),diff21_coef,Yc_to_YcTr_for_D2(:,1))
        !         call D2c_IBM_O2_MULTACC_3Dy(q3_y(n1s:n1e,n2s:n2e,n3s:n3e),tmp3_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),ysize_ibm_ccm(1),ysize_ibm_ccm(2),ysize_ibm_ccm(3),diff22_coef,Yc_to_YcTr_for_D2(:,2))

        !         call D1s_IBM_O2_MULTACC_3Dy(p23_y(n1s:n1e,n2s:n2e,n3s:n3e),tmp3_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),ysize_ibm_ccm(1),ysize_ibm_ccm(2),ysize_ibm_ccm(3),conv2_coef,Yc_to_YcTr_for_D1)

        !     endif

        !     ! call create_snapshot(COMMON_snapshot_path, "IBM", p23_y, "p23_y", 2)
        !     ! call create_snapshot(COMMON_snapshot_path, "IBM", tmp3_y_ibm, "tmp3_y_ibm", 2)
        !     ! call create_snapshot(COMMON_snapshot_path, "IBM", tmp3_y, "tmp3_y", 2)

        !     if (object_in_current_proc_cpm_y) then

        !         do i=ystart_ibm_cpm(1), yend_ibm_cpm(1)
        !             do j=ystart_ibm_cpm(2), yend_ibm_cpm(2)
        !                 do k=ystart_ibm_cpm(3), yend_ibm_cpm(3)

        !                     tmp1_y(i,j,k) = tmp1_y_ibm(i,j,k)
        !                     tmp2_y(i,j,k) = tmp2_y_ibm(i,j,k)
        !                     tmp3_y(i,j,k) = tmp3_y_ibm(i,j,k)

        !                 enddo
        !             enddo
        !         enddo

        !     endif

        ! endif


        ! if ((IBM_activated).and.(interpol_type == ANTISYMMETRIC_INTERPOL)) then

        !     do n=1,number_of_objects

        !         ! READING OBJECT FILE
        !         call read_object_size(trim(obj_file_path), nb_vertex, nb_faces)
        !         call read_object(trim(obj_file_path), vertex, faces, nb_vertex, nb_faces)

        !         ! Scale with H
        !         vertex(:,1)=vertex(:,1)*body_scale_x1(n)+body_x1(n)
        !         vertex(:,2)=vertex(:,2)*body_scale_x2(n)+body_x2(n)
        !         vertex(:,3)=vertex(:,3)*body_scale_x3(n)+body_x3(n)

        !         call get_objet_area(vertex, nb_vertex)

        !         if (object_in_current_proc_y_q3(n)) then

        !             ! Then, fix dp13/dy and d²q3/dy²
        !             j_bef = max(j_start_obj_q3(n)-1,1)
        !             j_aft = min(j_end_obj_q3(n)+1,n2-1)

        !             dx2_bef = (ys-Yc(j_bef))
        !             dx2_aft = (Yc(j_aft)-ye)

        !             if (j_start_obj_q3(n).gt.1) then

        !                 do i=ystart_ibm_q3(n,1), yend_ibm_q3(n,1)
        !                     do k=ystart_ibm_q3(n,3), yend_ibm_q3(n,3)

        !                         ! remove dp23/dz term computed by generic O2 scheme
        !                         tmp3_y(i,j_bef,k) = tmp3_y(i,j_bef,k) - ( (1.d0/conv2_coef)*(p23_y(i,j_start_obj_q3(n),k) - p23_y(i,j_bef,k) ) ) * Yc_to_YcTr_for_D1(j_bef)

        !                         ! remove d²q2/dz² term computed by generic O2 scheme
        !                         tmp3_y(i,j_bef,k) = tmp3_y(i,j_bef,k) - ( (0.5d0/diff21_coef)*(q3_y(i,j_start_obj_q3(n),k) - q3_y(i,j_bef-1,k)) ) * Yc_to_YcTr_for_D2(j_bef,1)
        !                         tmp3_y(i,j_bef,k) = tmp3_y(i,j_bef,k) - ( (1.d0/diff22_coef**2)*(q3_y(i,j_start_obj_q3(n),k) + q3_y(i,j_bef-1,k)) - 2.d0*(1.d0/diff1_coef**2)*q3_y(i,j_bef,k) ) * Yc_to_YcTr_for_D2(j_bef,2)

        !                         ! add dp23/dz term by taking into account the object position
        !                         tmp3_y(i,j_bef,k) = tmp3_y(i,j_bef,k) + (1.d0/(-dx2_bef))*(0.d0 - p23_y(i,j_bef,k))*Yc_to_YcTr_for_D1(j_bef)

        !                         ! remove d²q2/dz² term computed by generic O2 scheme
        !                         ! Caution, here we use the interpolated velocity to use the same scheme as O2 (with equal spacing)
        !                         tmp3_y(i,j_bef,k) = tmp3_y(i,j_bef,k) + ( (0.5d0/(dx2_bef*ren))*(0.d0 - ( q3_y(i,j_bef-1,k) + ((dx2-dx2_bef)/dx2)*(q3_y(i,j_bef,k)-q3_y(i,j_bef-1,k)) ) )) * Yc_to_YcTr_for_D2(j_bef,1)
        !                         tmp3_y(i,j_bef,k) = tmp3_y(i,j_bef,k) + ( (1.d0/(dx2_bef*sqrt(ren))**2)*(0.d0 + q3_y(i,j_bef-1,k) + ((dx2-dx2_bef)/dx2)*(q3_y(i,j_bef,k)-q3_y(i,j_bef-1,k)) ) - 2.d0*(1.d0/(dx2_bef*sqrt(ren))**2)*q3_y(i,j_bef,k) ) * Yc_to_YcTr_for_D2(j_bef,2)


        !                     enddo
        !                 enddo

        !             endif

        !             if (j_end_obj_q3(n).lt.(n2-1)) then

        !                 do i=ystart_ibm_q3(n,1), yend_ibm_q3(n,1)
        !                     do k=ystart_ibm_q3(n,3), yend_ibm_q3(n,3)

        !                         ! remove dp23/dz term computed by generic O2 scheme
        !                         tmp3_y(i,j_aft,k) = tmp3_y(i,j_aft,k) - ( (1.d0/conv2_coef)*(p23_y(i,j_aft+1,k) - p23_y(i,j_aft,k)) )*Yc_to_YcTr_for_D1(j_aft)

        !                         ! remove d²q2/dz² term computed by generic O2 scheme
        !                         tmp3_y(i,j_aft,k) = tmp3_y(i,j_aft,k) - ( (0.5d0/diff21_coef)*(q3_y(i,j_aft+1,k) - q3_y(i,j_end_obj_q3(n),k)) ) * Yc_to_YcTr_for_D2(j_aft,1)
        !                         tmp3_y(i,j_aft,k) = tmp3_y(i,j_aft,k) - ( (1.d0/diff22_coef**2)*(q3_y(i,j_start_obj_q3(n),k) + q3_y(i,j_aft+1,k)) - 2.d0*(1.d0/diff1_coef**2)*q3_y(i,j_aft,k) ) * Yc_to_YcTr_for_D2(j_aft,2)

        !                         ! add dp23/dz term by taking into account the object position
        !                         tmp3_y(i,j_aft,k) = tmp3_y(i,j_aft,k) + (1.d0/(-dx2_aft))*(p23_y(i,j_aft+1,k) - p23_y(i,j_aft,k))*Yc_to_YcTr_for_D1(j_aft)

        !                         ! remove d²q2/dz² term computed by generic O2 scheme
        !                         ! Caution, here we use the interpolated velocity to use the same scheme as O2 (with equal spacing)
        !                         tmp3_y(i,j_aft,k) = tmp3_y(i,j_aft,k) + ( (0.5d0/(dx2_bef*ren))*(0.d0 - ( q3_y(i,j_aft,k) + (dx2_bef/dx2)*(q3_y(i,j_aft+1,k)-q3_y(i,j_aft,k)) ) )) * Yc_to_YcTr_for_D2(j_aft,1)
        !                         tmp3_y(i,j_aft,k) = tmp3_y(i,j_aft,k) + ( (1.d0/(dx2_bef*sqrt(ren))**2)*(0.d0 - ( q3_y(i,j_aft,k) + (dx2_bef/dx2)*(q3_y(i,j_aft+1,k)-q3_y(i,j_aft,k)) )) - 2.d0*(1.d0/(dx2_bef*sqrt(ren))**2)*q3_y(i,j_aft,k) ) * Yc_to_YcTr_for_D2(j_aft,2)


        !                     enddo
        !                 enddo

        !             endif

        !         endif

        !         ! if (object_in_current_proc_x_q2(n)) then

        !         !     ! Then, fix dp12/dx and d²q2/dx²
        !         !     i_bef = i_start_obj_q2(n)-1
        !         !     i_aft = i_end_obj_q2(n)+1

        !         !     if (i_start_obj_q2(n).eq.1) i_bef=n1-1
        !         !     if (i_end_obj_q2(n).gt.(n1-1)) i_aft=1

        !         !     dx1_bef = (xs-Zc(i_bef))
        !         !     dx1_aft = (Zc(i_aft)-xe)

        !         !     if (i_start_obj_q2(n).eq.1) dx1_bef = (xs-Zc(1)-dx1)
        !         !     if (i_end_obj_q2(n).gt.(n1-1)) dx1_aft = (Zc(n1-1)+dx1-xe)

        !         !     do j=xstart_ibm_q2(n,2), xend_ibm_q2(n,2)
        !         !         do k=xstart_ibm_q2(n,3), xend_ibm_q2(n,3)

        !         !             ! remove dp23/dz term computed by generic O2 scheme
        !         !             tmp2_x(i_bef,j,k) = tmp2_x(i_bef,j,k) - ( (1.d0/conv1_coef)*(p12_x(i_start_obj_q2(n),j,k) - p12_x(i_bef,j,k)) )
        !         !             tmp2_x(i_aft,j,k) = tmp2_x(i_aft,j,k) - ( (1.d0/conv1_coef)*(p12_x(i_aft+1,j,k) - p12_x(i_aft,j,k)) )

        !         !             ! remove d²q2/dz² term computed by generic O2 scheme
        !         !             tmp2_x(i_bef,j,k) = tmp2_x(i_bef,j,k) - ( (1.d0/diff1_coef**2)*(q2_x(i_start_obj_q2(n),j,k) + q2_x(i_bef-1,j,k)) - 2.d0*(1.d0/diff1_coef**2)*q2_x(i_bef,j,k) )
        !         !             tmp2_x(i_aft,j,k) = tmp2_x(i_aft,j,k) - ( (1.d0/diff1_coef**2)*(q2_x(i_end_obj_q2(n),j,k) + q2_x(i_aft+1,j,k)) - 2.d0*(1.d0/diff1_coef**2)*q2_x(i_aft,j,k) )

        !         !             ! add dp23/dz term by taking into account the object position
        !         !             tmp2_x(i_bef,j,k) = tmp2_x(i_bef,j,k) + (1.d0/(-dx1_bef))*(0.d0 - p12_x(i_bef,j,k))
        !         !             tmp2_x(i_aft,j,k) = tmp2_x(i_aft,j,k) + (1.d0/(-dx1_aft))*(p12_x(i_aft+1,j,k) - p12_x(i_aft,j,k))

        !         !             ! remove d²q2/dz² term computed by generic O2 scheme
        !         !             ! Caution, here we use the interpolated velocity to use the same scheme as O2 (with equal spacing)
        !         !             tmp2_x(i_bef,j,k) = tmp2_x(i_bef,j,k) + (1.d0/(dx1_bef*sqrt(ren))**2)*(0.d0 + q2_x(i_bef-1,j,k) + ((dx1-dx1_bef)/dx1)*(q2_x(i_bef,j,k)-q2_x(i_bef-1,j,k))) - 2.d0*(1.d0/(dx1_bef*sqrt(ren))**2)*q2_x(i_bef,j,k)
        !         !             tmp2_x(i_aft,j,k) = tmp2_x(i_aft,j,k) + (1.d0/(dx1_aft*sqrt(ren))**2)*(0.d0 + q2_x(i_aft,j,k) + (dx1_bef/dx1)*(q2_x(i_aft+1,j,k)-q2_x(i_aft,j,k))) - 2.d0*(1.d0/(dx1_aft*sqrt(ren))**2)*q2_x(i_aft,j,k)


        !         !         enddo
        !         !     enddo

        !         ! endif

        !         if (object_in_current_proc_y_q1(n)) then

        !             ! Then, fix dp13/dy and d²q3/dy²
        !             j_bef = max(j_start_obj_q3(n)-1,1)
        !             j_aft = min(j_end_obj_q3(n)+1,n2-1)

        !             dx2_bef = (ys-Yc(j_bef))
        !             dx2_aft = (Yc(j_aft)-ye)

        !             if (j_start_obj_q1(n).gt.1) then

        !                 do i=ystart_ibm_q1(n,1), yend_ibm_q1(n,1)
        !                     do k=ystart_ibm_q1(n,3), yend_ibm_q1(n,3)

        !                         ! remove dp23/dz term computed by generic O2 scheme
        !                         tmp1_y(i,j_bef,k) = tmp1_y(i,j_bef,k) - ( (1.d0/conv2_coef)*(p23_y(i,j_start_obj_q1(n),k) - p23_y(i,j_bef,k)) )*Yc_to_YcTr_for_D1(j_bef)

        !                         ! remove d²q2/dz² term computed by generic O2 scheme
        !                         tmp1_y(i,j_bef,k) = tmp1_y(i,j_bef,k) - ( (0.5d0/diff21_coef)*(q1_y(i,j_start_obj_q1(n),k) - q1_y(i,j_bef-1,k)) ) * Yc_to_YcTr_for_D2(j_bef,1)
        !                         tmp1_y(i,j_bef,k) = tmp1_y(i,j_bef,k) - ( (1.d0/diff22_coef**2)*(q1_y(i,j_start_obj_q1(n),k) + q1_y(i,j_bef-1,k)) - 2.d0*(1.d0/diff1_coef**2)*q1_y(i,j_bef,k) ) * Yc_to_YcTr_for_D2(j_bef,2)

        !                         ! add dp23/dz term by taking into account the object position
        !                         tmp1_y(i,j_bef,k) = tmp1_y(i,j_bef,k) + (1.d0/(-dx2_bef))*(0.d0 - p23_y(i,j_bef,k))*Yc_to_YcTr_for_D1(j_bef)

        !                         ! remove d²q2/dz² term computed by generic O2 scheme
        !                         ! Caution, here we use the interpolated velocity to use the same scheme as O2 (with equal spacing)
        !                         tmp1_y(i,j_bef,k) = tmp1_y(i,j_bef,k) + ( (0.5d0/(dx2_bef*ren))*(0.d0 - ( q1_y(i,j_bef-1,k) + ((dx2-dx2_bef)/dx2)*(q1_y(i,j_bef,k)-q1_y(i,j_bef-1,k)) ) )) * Yc_to_YcTr_for_D2(j_bef,1)
        !                         tmp1_y(i,j_bef,k) = tmp1_y(i,j_bef,k) + ( (1.d0/(dx2_bef*sqrt(ren))**2)*(0.d0 + q1_y(i,j_bef-1,k) + ((dx2-dx2_bef)/dx2)*(q1_y(i,j_bef,k)-q1_y(i,j_bef-1,k)) ) - 2.d0*(1.d0/(dx2_bef*sqrt(ren))**2)*q1_y(i,j_bef,k) ) * Yc_to_YcTr_for_D2(j_bef,2)


        !                     enddo
        !                 enddo

        !             endif

        !             if (j_end_obj_q1(n).lt.(n2-1)) then

        !                 do i=ystart_ibm_q1(n,1), yend_ibm_q1(n,1)
        !                     do k=ystart_ibm_q1(n,3), yend_ibm_q1(n,3)

        !                         ! remove dp23/dz term computed by generic O2 scheme
        !                         tmp1_y(i,j_aft,k) = tmp1_y(i,j_aft,k) - ( (1.d0/conv2_coef)*(p23_y(i,j_aft+1,k) - p23_y(i,j_aft,k)) )*Yc_to_YcTr_for_D1(j_aft)

        !                         ! remove d²q2/dz² term computed by generic O2 scheme
        !                         tmp1_y(i,j_aft,k) = tmp1_y(i,j_aft,k) - ( (0.5d0/diff21_coef)*(q1_y(i,j_aft+1,k) - q1_y(i,j_end_obj_q1(n),k)) ) * Yc_to_YcTr_for_D2(j_aft,1)
        !                         tmp1_y(i,j_aft,k) = tmp1_y(i,j_aft,k) - ( (1.d0/diff22_coef**2)*(q1_y(i,j_start_obj_q1(n),k) + q1_y(i,j_aft+1,k)) - 2.d0*(1.d0/diff1_coef**2)*q1_y(i,j_aft,k) ) * Yc_to_YcTr_for_D2(j_aft,2)

        !                         ! add dp23/dz term by taking into account the object position
        !                         tmp1_y(i,j_aft,k) = tmp1_y(i,j_aft,k) + (1.d0/(-dx2_aft))*(p23_y(i,j_aft+1,k) - p23_y(i,j_aft,k))*Yc_to_YcTr_for_D1(j_aft)

        !                         ! remove d²q2/dz² term computed by generic O2 scheme
        !                         ! Caution, here we use the interpolated velocity to use the same scheme as O2 (with equal spacing)
        !                         tmp1_y(i,j_aft,k) = tmp1_y(i,j_aft,k) + ( (0.5d0/(dx2_bef*ren))*(0.d0 - ( q1_y(i,j_aft,k) + (dx2_bef/dx2)*(q1_y(i,j_aft+1,k)-q1_y(i,j_aft,k)) ) )) * Yc_to_YcTr_for_D2(j_aft,1)
        !                         tmp1_y(i,j_aft,k) = tmp1_y(i,j_aft,k) + ( (1.d0/(dx2_bef*sqrt(ren))**2)*(0.d0 - ( q1_y(i,j_aft,k) + (dx2_bef/dx2)*(q1_y(i,j_aft+1,k)-q1_y(i,j_aft,k)) )) - 2.d0*(1.d0/(dx2_bef*sqrt(ren))**2)*q1_y(i,j_aft,k) ) * Yc_to_YcTr_for_D2(j_aft,2)


        !                     enddo
        !                 enddo

        !             endif

        !         endif

        !     enddo

        ! endif

        if(BC2==NOSLIP) call diff_at_wall_2 ! implicit: NS_Q3_BC2=Dirichlet

        RHS1_y = RHS1_y + tmp1_y
        RHS2_y = RHS2_y + tmp2_y
        RHS3_y = RHS3_y + tmp3_y

        ! write(*,*) 'RHS1_y = ', RHS1_y(:,20,32)

        ! _______________________________________________________________________________________________________________________
        ! ***********************************************************************************************************************
        ! ****************************** Derivative along the X direction *******************************************************
        ! ***********************************************************************************************************************

        call transpose_z_to_y(p13_z, p13_y)
        call transpose_y_to_x(p13_y, p13_x)

        call transpose_y_to_x(p12_y, p12_x)

        call transpose_y_to_x(RHS1_y, RHS1_x)
        call transpose_y_to_x(RHS2_y, RHS2_x)
        call transpose_y_to_x(RHS3_y, RHS3_x)

        tmp1_x = 0.d0
        tmp2_x = 0.d0
        tmp3_x = 0.d0

        n2e=xsize(2)!(min(n2m, xend(2))-xstart(2))+1
        n3e=xsize(3)!(min(n3m, xend(3))-xstart(3))+1

        ! For RHS1
        call D0s_3Dx(q1_x, q1q1_x, n1,xsize(2),n2e,xsize(3),n3e, NS_Q1_BC1)
        q1q1_x=q1q1_x**2
        call D1ssh_ACC_3Dx(q1q1_x, tmp1_x, n1,xsize(2),n2e,xsize(3),n3e, conv1_coef, .true., NS_P11_BC1)
        call D2c_ACC_3Dx(q1_x, tmp1_x, n1,xsize(2),n2e,xsize(3),n3e, diff1_coef, .false., NS_Q1_BC1)

        ! For RHS2
        call D1s_ACC_3Dx(p12_x, tmp2_x, n1,xsize(2),n2e,xsize(3),n3e, conv1_coef, .false., NS_P12_BC1)
        call D2c_ACC_3Dx(q2_x, tmp2_x, n1,xsize(2),n2e,xsize(3),n3e, diff1_coef, .true., NS_Q2_BC1)

        ! For RHS3
        call D1s_ACC_3Dx(p13_x, tmp3_x, n1,xsize(2),n2e,xsize(3),n3e, conv1_coef, .false., NS_P13_BC1)
        call D2c_ACC_3Dx(q3_x, tmp3_x, n1,xsize(2),n2e,xsize(3),n3e, diff1_coef, .true., NS_Q3_BC1)

        ! If the IBM is activated then the DRP scheme has to be corrected ...
        ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) then

        !     if (object_in_current_proc_ccm_x) then

        !         n1s=xstart_ibm_ccm(1)
        !         n2s=xstart_ibm_ccm(2)
        !         n3s=xstart_ibm_ccm(3)

        !         n1e=xend_ibm_ccm(1)
        !         n2e=xend_ibm_ccm(2)
        !         n3e=xend_ibm_ccm(3)

        !         ! For RHS1
        !         call D2c_IBM_O2_3Dx(q1_x(n1s:n1e,n2s:n2e,n3s:n3e),tmp1_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e),xsize_ibm_ccm(1),xsize_ibm_ccm(2),xsize_ibm_ccm(3),diff1_coef)

        !         call D0s_IBM_O2_3Dx(q1_x(n1s:n1e,n2s:n2e,n3s:n3e),q1q1_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e),xsize_ibm_ccm(1),xsize_ibm_ccm(2),xsize_ibm_ccm(3))
        !         q1q1_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e)=q1q1_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e)**2

        !         call D1ssh_IBM_O2_ACC_3Dx(q1q1_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp1_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e),xsize_ibm_ccm(1),xsize_ibm_ccm(2),xsize_ibm_ccm(3),conv1_coef)

        !         ! For RHS2
        !         call D2c_IBM_O2_3Dx(q2_x(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e),xsize_ibm_ccm(1),xsize_ibm_ccm(2),xsize_ibm_ccm(3),diff1_coef)
        !         call D1s_IBM_O2_ACC_3Dx(p12_x(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e),xsize_ibm_ccm(1),xsize_ibm_ccm(2),xsize_ibm_ccm(3),conv1_coef)

        !         ! For RHS3
        !         call D2c_IBM_O2_3Dx(q3_x(n1s:n1e,n2s:n2e,n3s:n3e),tmp3_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e),xsize_ibm_ccm(1),xsize_ibm_ccm(2),xsize_ibm_ccm(3),diff1_coef)
        !         call D1s_IBM_O2_ACC_3Dx(p13_x(n1s:n1e,n2s:n2e,n3s:n3e),tmp3_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e),xsize_ibm_ccm(1),xsize_ibm_ccm(2),xsize_ibm_ccm(3),conv1_coef)

        !     endif

        !     if (object_in_current_proc_cpm_x) then

        !         do i=xstart_ibm_cpm(1), xend_ibm_cpm(1)
        !             do j=xstart_ibm_cpm(2), xend_ibm_cpm(2)
        !                 do k=xstart_ibm_cpm(3), xend_ibm_cpm(3)

        !                     tmp1_x(i,j,k) = tmp1_x_ibm(i,j,k)
        !                     tmp2_x(i,j,k) = tmp2_x_ibm(i,j,k)
        !                     tmp3_x(i,j,k) = tmp3_x_ibm(i,j,k)

        !                 enddo
        !             enddo
        !         enddo

        !     endif

        ! endif


        ! if ((IBM_activated).and.(interpol_type == ANTISYMMETRIC_INTERPOL)) then

        !     do n=1,number_of_objects

        !         ! READING OBJECT FILE
        !         call read_object_size(trim(obj_file_path), nb_vertex, nb_faces)
        !         call read_object(trim(obj_file_path), vertex, faces, nb_vertex, nb_faces)

        !         ! Scale with H
        !         vertex(:,1)=vertex(:,1)*body_scale_x1(n)+body_x1(n)
        !         vertex(:,2)=vertex(:,2)*body_scale_x2(n)+body_x2(n)
        !         vertex(:,3)=vertex(:,3)*body_scale_x3(n)+body_x3(n)

        !         call get_objet_area(vertex, nb_vertex)

        !         xs = max(0.d0,xs)
        !         xe = min(L1,xe)
        !         ys = max(0.d0,ys)
        !         ye = min(L2,ye)
        !         zs = max(0.d0,zs)
        !         ze = min(L3,ze)

        !         if (object_in_current_proc_x_q3(n)) then

        !             i0 = i_start_obj_q3(n)
        !             i1 = i_end_obj_q3(n)
        !             j0 = j_start_obj_q3(n)
        !             j1 = j_end_obj_q3(n)
        !             k0 = k_start_obj_q3(n)
        !             k1 = k_end_obj_q3(n)

        !             i_bef = i0-1
        !             i_aft = i1+1

        !             if (i0.eq.1) i_bef=n1-1
        !             if (i1.gt.(n1-1)) i_aft=1

        !             dx1_bef = abs(xs-Zc(i0)-dx1)
        !             dx1_aft = abs(Zc(i1)+dx1-xe)

        !             do j=xstart_ibm_q3(n,2), xend_ibm_q3(n,2)
        !                 do k=xstart_ibm_q3(n,3), xend_ibm_q3(n,3)

        !                     if ((j.ge.j0).and.(j.le.j1).and.(k.ge.k0).and.(k.le.k1)) then

        !                         ! correct dp12/dx term only for i_bef (no need to do this for i_aft)
        !                         ! remove what is computed by generic O2 scheme
        !                         tmp3_x(i_bef,j,k) = tmp3_x(i_bef,j,k) - ( (1.d0/conv1_coef)*(p13_x(i0,j,k) - p13_x(i_bef,j,k)) )
        !                         ! add what is computed by taking into account the object position
        !                         tmp3_x(i_bef,j,k) = tmp3_x(i_bef,j,k) + (1.d0/(-dx1_bef))*(0.d0 - p13_x(i_bef,j,k))

        !                         ! correct d²u/dx² term (no need to do this for i_bef)
        !                         ! remove what is computed by generic O2 scheme
        !                         tmp3_x(i_bef,j,k) = tmp3_x(i_bef,j,k) - ( (1.d0/diff1_coef**2)*(q3_x(i_bef-1,j,k) -2.d0*q3_x(i_bef,j,k) + q3_x(i0,j,k)) )
        !                         tmp3_x(i_aft,j,k) = tmp3_x(i_aft,j,k) - ( (1.d0/diff1_coef**2)*(q3_x(i_aft+1,j,k) -2.d0*q3_x(i_aft,j,k) + q3_x(i1,j,k)) )
        !                         ! add what is computed by taking into account the object position
        !                         ! Caution, here we use the interpolated velocity to use the same scheme as O2 (with equal spacing)
        !                         tmp3_x(i_bef,j,k) = tmp3_x(i_bef,j,k) + (1.d0/(dx1_bef*sqrt(ren))**2)*(0.d0 + q3_x(i_bef-1,j,k) + ((dx1-dx1_bef)/dx1)*(q3_x(i_bef,j,k)-q3_x(i_bef-1,j,k))) - 2.d0*(1.d0/(dx1_bef*sqrt(ren))**2)*q3_x(i_bef,j,k)
        !                         tmp3_x(i_aft,j,k) = tmp3_x(i_aft,j,k) + (1.d0/(dx1_aft*sqrt(ren))**2)*(0.d0 + q3_x(i_aft,j,k) + (dx1_bef/dx1)*(q3_x(i_aft+1,j,k)-q3_x(i_aft,j,k))) - 2.d0*(1.d0/(dx1_aft*sqrt(ren))**2)*q3_x(i_aft,j,k)

        !                     endif

        !                 enddo
        !             enddo

        !         endif

        !         if (object_in_current_proc_x_q2(n)) then

        !             i0 = i_start_obj_q2(n)
        !             i1 = i_end_obj_q2(n)
        !             j0 = j_start_obj_q2(n)
        !             j1 = j_end_obj_q2(n)
        !             k0 = k_start_obj_q2(n)
        !             k1 = k_end_obj_q2(n)

        !             i_bef = i0-1
        !             i_aft = i1+1

        !             if (i0.eq.1) i_bef=n1-1
        !             if (i1.gt.(n1-1)) i_aft=1

        !             dx1_bef = abs(xs-Zc(i0)-dx1)
        !             dx1_aft = abs(Zc(i1)+dx1-xe)

        !             do j=xstart_ibm_q2(n,2), xend_ibm_q2(n,2)
        !                 do k=xstart_ibm_q2(n,3), xend_ibm_q2(n,3)

        !                     if ((j.ge.j0).and.(j.le.j1).and.(k.ge.k0).and.(k.le.k1)) then

        !                         ! correct dp12/dx term only for i_bef (no need to do this for i_aft)
        !                         ! remove what is computed by generic O2 scheme
        !                         tmp2_x(i_bef,j,k) = tmp2_x(i_bef,j,k) - ( (1.d0/conv1_coef)*(p12_x(i0,j,k) - p12_x(i_bef,j,k)) )
        !                         ! add what is computed by taking into account the object position
        !                         tmp2_x(i_bef,j,k) = tmp2_x(i_bef,j,k) + (1.d0/(-dx1_bef))*(0.d0 - p12_x(i_bef,j,k))

        !                         ! correct d²u/dx² term (no need to do this for i_bef)
        !                         ! remove what is computed by generic O2 scheme
        !                         tmp2_x(i_bef,j,k) = tmp2_x(i_bef,j,k) - ( (1.d0/diff1_coef**2)*(q2_x(i_bef-1,j,k) -2.d0*q2_x(i_bef,j,k) + q2_x(i0,j,k)) )
        !                         tmp2_x(i_aft,j,k) = tmp2_x(i_aft,j,k) - ( (1.d0/diff1_coef**2)*(q2_x(i_aft+1,j,k) -2.d0*q2_x(i_aft,j,k) + q2_x(i1,j,k)) )
        !                         ! add what is computed by taking into account the object position
        !                         ! Caution, here we use the interpolated velocity to use the same scheme as O2 (with equal spacing)
        !                         tmp2_x(i_bef,j,k) = tmp2_x(i_bef,j,k) + (1.d0/(dx1_bef*sqrt(ren))**2)*(0.d0 + q2_x(i_bef-1,j,k) + ((dx1-dx1_bef)/dx1)*(q2_x(i_bef,j,k)-q2_x(i_bef-1,j,k))) - 2.d0*(1.d0/(dx1_bef*sqrt(ren))**2)*q2_x(i_bef,j,k)
        !                         tmp2_x(i_aft,j,k) = tmp2_x(i_aft,j,k) + (1.d0/(dx1_aft*sqrt(ren))**2)*(0.d0 + q2_x(i_aft,j,k) + (dx1_bef/dx1)*(q2_x(i_aft+1,j,k)-q2_x(i_aft,j,k))) - 2.d0*(1.d0/(dx1_aft*sqrt(ren))**2)*q2_x(i_aft,j,k)

        !                     endif

        !                 enddo
        !             enddo

        !         endif

        !         if (object_in_current_proc_x_q1(n)) then

        !             i0 = i_start_obj_q1(n)
        !             i1 = i_end_obj_q1(n)
        !             j0 = j_start_obj_q1(n)
        !             j1 = j_end_obj_q1(n)
        !             k0 = k_start_obj_q1(n)
        !             k1 = k_end_obj_q1(n)

        !             i_bef = i0-1
        !             i_aft = i1+1

        !             if (i0.eq.1) i_bef=n1-1
        !             if (i1.gt.(n1-1)) i_aft=1

        !             dx1_bef = abs(xs-Z(i0)-dx1)
        !             dx1_aft = abs(Z(i1)+dx1-xe)

        !             do j=xstart_ibm_q1(n,2), xend_ibm_q1(n,2)
        !                 do k=xstart_ibm_q1(n,3), xend_ibm_q1(n,3)

        !                     if ((j.ge.j0).and.(j.le.j1).and.(k.ge.k0).and.(k.le.k1)) then

        !                         ! correct dq1q1/dx term only for i_aft (no need to do this for i_bef)
        !                         ! remove what is computed by generic O2 scheme
        !                         tmp1_x(i_aft,j,k) = tmp1_x(i_aft,j,k) - ( (1.d0/conv1_coef)*(q1q1_x(i_aft,j,k) - q1q1_x(i1,j,k)) )
        !                         ! add what is computed by taking into account the object position
        !                         tmp1_x(i_aft,j,k) = tmp1_x(i_aft,j,k) + (1.d0/(-dx1_aft))*(q1q1_x(i_aft,j,k) - 0.d0)

        !                         ! correct d²u/dx² term
        !                         ! remove what is computed by generic O2 scheme
        !                         tmp1_x(i_bef,j,k) = tmp1_x(i_bef,j,k) - ( (1.d0/diff1_coef**2)*(q1_x(i_bef-1,j,k) -2.d0*q1_x(i_bef,j,k) + q1_x(i0,j,k)) )
        !                         tmp1_x(i_aft,j,k) = tmp1_x(i_aft,j,k) - ( (1.d0/diff1_coef**2)*(q1_x(i_aft+1,j,k) -2.d0*q1_x(i_aft,j,k) + q1_x(i1,j,k)) )
        !                         ! add what is computed by taking into account the object position
        !                         ! Caution, here we use the interpolated velocity to use the same scheme as O2 (with equal spacing)
        !                         tmp1_x(i_bef,j,k) = tmp1_x(i_bef,j,k) + (1.d0/(dx1_bef*sqrt(ren))**2)*(0.d0 + q1_x(i_bef-1,j,k) + ((dx1-dx1_bef)/dx1)*(q1_x(i_bef,j,k)-q1_x(i_bef-1,j,k))) - 2.d0*(1.d0/(dx1_bef*sqrt(ren))**2)*q1_x(i_bef,j,k)
        !                         tmp1_x(i_aft,j,k) = tmp1_x(i_aft,j,k) + (1.d0/(dx1_aft*sqrt(ren))**2)*(0.d0 + q1_x(i_aft,j,k) + (dx1_bef/dx1)*(q1_x(i_aft+1,j,k)-q1_x(i_aft,j,k))) - 2.d0*(1.d0/(dx1_aft*sqrt(ren))**2)*q1_x(i_aft,j,k)

        !                     endif

        !                 enddo
        !             enddo

        !         endif

        !     enddo

        ! endif

        if((BC1==NOSLIP).or.(BC1==PSEUDO_PERIODIC).or.(BC1==OPEN)) call diff_at_wall_1   ! ATTENTION

        RHS1_x = RHS1_x + tmp1_x
        RHS2_x = RHS2_x + tmp2_x
        RHS3_x = RHS3_x + tmp3_x

        return

    contains

        subroutine diff_at_wall_1()

            ! d2f=a1*f0 + a2*f1 + a3*f2
            ! Regular case: down    a1=2/3, a2=-1 a3=1/3
            ! Regular case: up      a1=1/3, a2=-1 a3=2/3
            ! ATTENTION
            do k=xstart(3), min(n3m, xend(3))       !do k=1,n3m
                do j=xstart(2), min(n2m, xend(2))       !do j=1,n2m
                    tmp2_x(1,j,k)   =tmp2_x(1,j,k)  +( q2_x(2,j,k)       /3.d0 - q2_x(1,j,k)      + q2_wall10(j,k)*2.d0/3.d0)  /diff1_coef**2
                    tmp2_x(n1m,j,k) =tmp2_x(n1m,j,k)+( q2_x(n1m-1,j,k)   /3.d0 - q2_x(n1m,j,k)    + q2_wall11(j,k)*2.d0/3.d0 ) /diff1_coef**2

                    tmp3_x(1,j,k)   =tmp3_x(1,j,k)  +( q3_x(2,j,k)       /3.d0 - q3_x(1,j,k)      + q3_wall10(j,k)*2.d0/3.d0 ) /diff1_coef**2
                    tmp3_x(n1m,j,k) =tmp3_x(n1m,j,k)+( q3_x(n1m-1,j,k)   /3.d0 - q3_x(n1m,j,k)    + q3_wall11(j,k)*2.d0/3.d0)  /diff1_coef**2
                enddo
            enddo

        end subroutine diff_at_wall_1

        subroutine diff_at_wall_1_v2()

            ! d2f=a1*f0 + a2*f1 + a3*f2
            ! Regular case: down    a1=2/3, a2=-1 a3=1/3
            ! Regular case: up      a1=1/3, a2=-1 a3=2/3
            ! ATTENTION
            RHS2_x(1,j,k)   =RHS2_x(1,j,k)  +( q2_x(2,j,k) - 2.d0*q2_x(1,j,k)   + q2_x(n1m,j,k)  )     /diff1_coef**2
            RHS2_x(n1m,j,k) =RHS2_x(n1m,j,k)+( q2_x(1,j,k) - 2.d0*q2_x(n1m,j,k) + q2_x(n1m-1,j,k))     /diff1_coef**2

            RHS3_x(1,j,k)   =RHS3_x(1,j,k)  +( q3_x(2,j,k) - 2.d0*q3_x(1,j,k)   + q3_x(n1m,j,k)  )     /diff1_coef**2
            RHS3_x(n1m,j,k) =RHS3_x(n1m,j,k)+( q3_x(1,j,k) - 2.d0*q3_x(n1m,j,k) + q3_x(n1m-1,j,k))     /diff1_coef**2

        end subroutine diff_at_wall_1_v2

        subroutine diff_at_wall_2()

            do k=ystart(3), min(n3m, yend(3))       !do k=1,n3m
                do i=ystart(1), min(n1m, yend(1))       !do i=1,n1m

                    tmp1_y(i,1,k)=tmp1_y(i,1,k)+( q1_y(i,2,k)*a3_d + q1_y(i,1,k)*a2_d + q1_wall20(i,k)*a1_d)            /diff2_coef**2
                    tmp1_y(i,n2m,k)=tmp1_y(i,n2m,k)+( q1_y(i,n2m-1,k)*a1_u + q1_y(i,n2m,k)*a2_u + q1_wall21(i,k)*a3_u ) /diff2_coef**2

                    tmp3_y(i,1,k)=tmp3_y(i,1,k)+( q3_y(i,2,k)*a3_d + q3_y(i,1,k)*a2_d + q3_wall20(i,k)*a1_d )           /diff2_coef**2
                    tmp3_y(i,n2m,k)=tmp3_y(i,n2m,k)+( q3_y(i,n2m-1,k)*a1_u + q3_y(i,n2m,k) *a2_u + q3_wall21(i,k)*a3_u) /diff2_coef**2
                enddo
            enddo

        end subroutine diff_at_wall_2

        subroutine diff_at_wall_3()

            ! d2f=a1*f0 + a2*f1 + a3*f2
            ! Regular case: down    a1=2/3, a2=-1 a3=1/3
            ! Regular case: up      a1=1/3, a2=-1 a3=2/3
            ! ATTENTION

            ! ************* CORR_JO ******************* !
            ! n3m au lieu de n2m dans RHS_z(i,j,n3m) ** !
            ! ***************************************** !

            do j = zstart(2), min(n2m, zend(2))     !do j=1,n2m
                do i=zstart(1), min(n1m, zend(1))       !do i=1,n1m
                    tmp1_z(i,j,1)   =tmp1_z(i,j,1)  +( q1_z(i,j,2)       /3.d0 - q1_z(i,j,1)      + q1_wall30(i,j)*2.d0/3.d0 ) /diff3_coef**2
                    tmp1_z(i,j,n3m) =tmp1_z(i,j,n3m)+( q1_z(i,j,n3m-1)   /3.d0 - q1_z(i,j,n3m)    + q1_wall31(i,j)*2.d0/3.d0)  /diff3_coef**2

                    tmp2_z(i,j,1)   =tmp2_z(i,j,1)  +( q2_z(i,j,2)       /3.d0 - q2_z(i,j,1)      + q2_wall30(i,j)*2.d0/3.d0)  /diff3_coef**2
                    tmp2_z(i,j,n3m) =tmp2_z(i,j,n3m)+( q2_z(i,j,n3m-1)   /3.d0 - q2_z(i,j,n3m)    + q2_wall31(i,j)*2.d0/3.d0 ) /diff3_coef**2
                enddo
            enddo

        end subroutine diff_at_wall_3

    end subroutine

    !GENERIC
    subroutine perform_RHS_ibm()

        use boundary_scheme
        use DRP_IBM
        ! use DRP_IBM_fine
        use IBM

        implicit none

        !real*8,dimension(n3)           :: uu
        real*8,dimension(:,:,:), allocatable           :: q1q1_x_ibm, tmp1_x_ibm, tmp2_x_ibm, tmp3_x_ibm
        real*8,dimension(:,:,:), allocatable           :: q2q2_y_ibm, tmp1_y_ibm, tmp2_y_ibm, tmp3_y_ibm
        real*8,dimension(:,:,:), allocatable           :: q3q3_z_ibm, tmp1_z_ibm, tmp2_z_ibm, tmp3_z_ibm

        integer :: k,i,j, mpi_err, n1e, n2e, n3e, n1s, n2s, n3s

        real*8 :: t_bef, t_aft
        
        t_bef = MPI_WTIME()

        ! _______________________________________________________________________________________________________________________
        ! ***********************************************************************************************************************
        ! ****************************** Derivative along the Z direction *******************************************************
        ! ***********************************************************************************************************************

        ! ******** CORR_JO ********* !
        ! REinitialisation des RHS_z !
        ! ************************** !
        RHS1_z_ibm=0.d0
        RHS2_z_ibm=0.d0
        RHS3_z_ibm=0.d0

        if (object_in_current_proc_fcm_z) then

            ! For RHS1
            allocate(tmp1_z_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2), zstart_ibm_fcm(3):zend_ibm_fcm(3)))
            tmp1_z_ibm = 0.d0
            ! For RHS2
            allocate(tmp2_z_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2), zstart_ibm_fcm(3):zend_ibm_fcm(3)))
            tmp2_z_ibm = 0.d0
            ! For RHS3
            allocate(tmp3_z_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2), zstart_ibm_fcm(3):zend_ibm_fcm(3)))
            tmp3_z_ibm = 0.d0

            allocate(q3q3_z_ibm(zstart_ibm_fcm(1):zend_ibm_fcm(1), zstart_ibm_fcm(2):zend_ibm_fcm(2), zstart_ibm_fcm(3):zend_ibm_fcm(3)))
            q3q3_z_ibm = 0.d0

            n1s=zstart_ibm_fcm(1)
            n2s=zstart_ibm_fcm(2)
            n3s=zstart_ibm_fcm(3)

            n1e=zend_ibm_fcm(1)
            n2e=zend_ibm_fcm(2)
            n3e=zend_ibm_fcm(3)

            ! For RHS1
            call D2c_IBM_O2_3Dz(q1_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp1_z_ibm,zsize_ibm_fcm(1),zsize_ibm_fcm(2),zsize_ibm_fcm(3),diff3_coef_ibm)
            ! call D1s_IBM_O2_ACC_3Dz(p13_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp1_z_ibm,n1s,n1e,n2s,n2e,n3s,n3e,conv3_coef_ibm)
            call D1c_IBM_O2_ACC_3Dz(p13_z_q1_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp1_z_ibm,zsize_ibm_fcm(1),zsize_ibm_fcm(2),zsize_ibm_fcm(3),conv3_coef_ibm)

            ! For RHS2
            call D2c_IBM_O2_3Dz(q2_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_z_ibm,zsize_ibm_fcm(1),zsize_ibm_fcm(2),zsize_ibm_fcm(3),diff3_coef_ibm)
            ! call D1s_IBM_O2_ACC_3Dz(p23_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_z_ibm,n1s,n1e,n2s,n2e,n3s,n3e,conv3_coef_ibm)
            call D1c_IBM_O2_ACC_3Dz(p23_z_q2_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_z_ibm,zsize_ibm_fcm(1),zsize_ibm_fcm(2),zsize_ibm_fcm(3),conv3_coef_ibm)

            ! For RHS1
            call D2c_IBM_O2_3Dz(q3_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp3_z_ibm,zsize_ibm_fcm(1),zsize_ibm_fcm(2),zsize_ibm_fcm(3),diff3_coef_ibm)

            call D0s_IBM_O2_3Dz(q3_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e),q3q3_z_ibm,zsize_ibm_fcm(1),zsize_ibm_fcm(2),zsize_ibm_fcm(3))
            q3q3_z_ibm=q3q3_z_ibm**2

            call D1ssh_IBM_O2_ACC_3Dz(q3q3_z_ibm,tmp3_z_ibm,zsize_ibm_fcm(1),zsize_ibm_fcm(2),zsize_ibm_fcm(3),conv3_coef_ibm)

            if (BC3==NOSLIP) call diff_at_wall_3_ibm   ! ATTENTION

            do i=zstart_ibm_fcm(1),zend_ibm_fcm(1)
                do j=zstart_ibm_fcm(2),zend_ibm_fcm(2)
                    do k=zstart_ibm_fcm(3),zend_ibm_fcm(3)

                        RHS1_z_ibm(i,j,k) = RHS1_z_ibm(i,j,k) + tmp1_z_ibm(i,j,k)
                        RHS2_z_ibm(i,j,k) = RHS2_z_ibm(i,j,k) + tmp2_z_ibm(i,j,k)
                        RHS3_z_ibm(i,j,k) = RHS3_z_ibm(i,j,k) + tmp3_z_ibm(i,j,k)

                    enddo
                enddo
            enddo

        endif



        t_aft = MPI_WTIME()
        ! write(*,*) 'finish dz', t_aft-t_bef

        t_bef = MPI_WTIME()

        ! _______________________________________________________________________________________________________________________
        ! ***********************************************************************************************************************
        ! ****************************** Derivative along the Y direction *******************************************************
        ! ***********************************************************************************************************************

        call transpose_z_to_y(RHS1_z_ibm, RHS1_y_ibm, decomp_fine)
        call transpose_z_to_y(RHS2_z_ibm, RHS2_y_ibm, decomp_fine)
        call transpose_z_to_y(RHS3_z_ibm, RHS3_y_ibm, decomp_fine)

        if (object_in_current_proc_fcm_y) then

            ! For RHS1
            allocate(tmp1_y_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(2):yend_ibm_fcm(2), ystart_ibm_fcm(3):yend_ibm_fcm(3)))
            tmp1_y_ibm = 0.d0

            ! For RHS2
            allocate(tmp2_y_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(2):yend_ibm_fcm(2), ystart_ibm_fcm(3):yend_ibm_fcm(3)))
            tmp2_y_ibm = 0.d0

            allocate(q2q2_y_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(2):yend_ibm_fcm(2), ystart_ibm_fcm(3):yend_ibm_fcm(3)))
            q2q2_y_ibm = 0.d0

            ! For RHS3
            allocate(tmp3_y_ibm(ystart_ibm_fcm(1):yend_ibm_fcm(1), ystart_ibm_fcm(2):yend_ibm_fcm(2), ystart_ibm_fcm(3):yend_ibm_fcm(3)))
            tmp3_y_ibm = 0.d0

            n1s=ystart_ibm_fcm(1)
            n2s=ystart_ibm_fcm(2)
            n3s=ystart_ibm_fcm(3)

            n1e=yend_ibm_fcm(1)
            n2e=yend_ibm_fcm(2)
            n3e=yend_ibm_fcm(3)

            ! For RHS1
            call D1c_IBM_O2_MULTACC_3Dy(q1_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp1_y_ibm,ysize_ibm_fcm(1),ysize_ibm_fcm(2),ysize_ibm_fcm(3),diff21_coef_ibm,Y_to_YTr_for_D2_ibm(:,1))
            call D2c_IBM_O2_MULTACC_3Dy(q1_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp1_y_ibm,ysize_ibm_fcm(1),ysize_ibm_fcm(2),ysize_ibm_fcm(3),diff22_coef_ibm,Y_to_YTr_for_D2_ibm(:,2))

            ! call D1s_IBM_O2_MULTACC_3Dy(p12_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp1_y_ibm,n1s,n1e,n2s,n2e,n3s,n3e,conv2_coef_ibm,Y_to_YTr_for_D1_ibm)
            call D1c_IBM_O2_MULTACC_3Dy(p12_y_q1_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp1_y_ibm,ysize_ibm_fcm(1),ysize_ibm_fcm(2),ysize_ibm_fcm(3),conv2_coef_ibm,Y_to_YTr_for_D1_ibm)

            ! For RHS2
            call D1c_IBM_O2_MULTACC_3Dy(q2_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_y_ibm,ysize_ibm_fcm(1),ysize_ibm_fcm(2),ysize_ibm_fcm(3),diff21_coef_ibm,Y_to_YTr_for_D2_ibm(:,1))
            call D2c_IBM_O2_MULTACC_3Dy(q2_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_y_ibm,ysize_ibm_fcm(1),ysize_ibm_fcm(2),ysize_ibm_fcm(3),diff22_coef_ibm,Y_to_YTr_for_D2_ibm(:,2))

            call D0s_IBM_O2_3Dy(q2_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),q2q2_y_ibm,ysize_ibm_fcm(1),ysize_ibm_fcm(2),ysize_ibm_fcm(3))
            q2q2_y_ibm=q2q2_y_ibm**2

            call D1ssh_IBM_O2_MULTACC_3Dy(q2q2_y_ibm,tmp2_y_ibm,ysize_ibm_fcm(1),ysize_ibm_fcm(2),ysize_ibm_fcm(3),conv2_coef_ibm,Y_to_YTr_for_D1_ibm)

            ! For RHS3
            call D1c_IBM_O2_MULTACC_3Dy(q3_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp3_y_ibm,ysize_ibm_fcm(1),ysize_ibm_fcm(2),ysize_ibm_fcm(3),diff21_coef_ibm,Y_to_YTr_for_D2_ibm(:,1))
            call D2c_IBM_O2_MULTACC_3Dy(q3_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp3_y_ibm,ysize_ibm_fcm(1),ysize_ibm_fcm(2),ysize_ibm_fcm(3),diff22_coef_ibm,Y_to_YTr_for_D2_ibm(:,2))

            ! call D1s_IBM_O2_MULTACC_3Dy(p23_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp3_y_ibm,n1s,n1e,n2s,n2e,n3s,n3e,conv2_coef_ibm,Y_to_YTr_for_D1_ibm)
            call D1c_IBM_O2_MULTACC_3Dy(p23_y_q3_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp3_y_ibm,ysize_ibm_fcm(1),ysize_ibm_fcm(2),ysize_ibm_fcm(3),conv2_coef_ibm,Y_to_YTr_for_D1_ibm)

            if(BC2==NOSLIP) call diff_at_wall_2_ibm ! implicit: NS_Q3_BC2=Dirichlet

            do i=ystart_ibm_fcm(1),yend_ibm_fcm(1)
                do j=ystart_ibm_fcm(2),yend_ibm_fcm(2)
                    do k=ystart_ibm_fcm(3),yend_ibm_fcm(3)

                        RHS1_y_ibm(i,j,k) = RHS1_y_ibm(i,j,k) + tmp1_y_ibm(i,j,k)
                        RHS2_y_ibm(i,j,k) = RHS2_y_ibm(i,j,k) + tmp2_y_ibm(i,j,k)
                        RHS3_y_ibm(i,j,k) = RHS3_y_ibm(i,j,k) + tmp3_y_ibm(i,j,k)

                    enddo
                enddo
            enddo

        endif

        t_aft = MPI_WTIME()
        ! write(*,*) 'finish dy', t_aft-t_bef

        t_bef = MPI_WTIME()

        ! _______________________________________________________________________________________________________________________
        ! ***********************************************************************************************************************
        ! ****************************** Derivative along the X direction *******************************************************
        ! ***********************************************************************************************************************

        call transpose_y_to_x(RHS1_y_ibm, RHS1_x_ibm, decomp_fine)
        call transpose_y_to_x(RHS2_y_ibm, RHS2_x_ibm, decomp_fine)
        call transpose_y_to_x(RHS3_y_ibm, RHS3_x_ibm, decomp_fine)

        if (object_in_current_proc_fcm_x) then

            ! For RHS1
            allocate(tmp1_x_ibm(xstart_ibm_fcm(1):xend_ibm_fcm(1), xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)))
            tmp1_x_ibm = 0.d0

            allocate(q1q1_x_ibm(xstart_ibm_fcm(1):xend_ibm_fcm(1), xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)))
            q1q1_x_ibm = 0.d0

            ! For RHS2
            allocate(tmp2_x_ibm(xstart_ibm_fcm(1):xend_ibm_fcm(1), xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)))
            tmp2_x_ibm = 0.d0

            ! For RHS3
            allocate(tmp3_x_ibm(xstart_ibm_fcm(1):xend_ibm_fcm(1), xstart_ibm_fcm(2):xend_ibm_fcm(2), xstart_ibm_fcm(3):xend_ibm_fcm(3)))
            tmp3_x_ibm = 0.d0

            n1s=xstart_ibm_fcm(1)
            n2s=xstart_ibm_fcm(2)
            n3s=xstart_ibm_fcm(3)

            n1e=xend_ibm_fcm(1)
            n2e=xend_ibm_fcm(2)
            n3e=xend_ibm_fcm(3)

            ! For RHS1
            call D2c_IBM_O2_3Dx(q1_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp1_x_ibm,xsize_ibm_fcm(1),xsize_ibm_fcm(2),xsize_ibm_fcm(3),diff1_coef_ibm)

            call D0s_IBM_O2_3Dx(q1_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e),q1q1_x_ibm,xsize_ibm_fcm(1),xsize_ibm_fcm(2),xsize_ibm_fcm(3))
            q1q1_x_ibm=q1q1_x_ibm**2

            call D1ssh_IBM_O2_ACC_3Dx(q1q1_x_ibm,tmp1_x_ibm,xsize_ibm_fcm(1),xsize_ibm_fcm(2),xsize_ibm_fcm(3),conv1_coef_ibm)

            ! For RHS2
            call D2c_IBM_O2_3Dx(q2_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_x_ibm,xsize_ibm_fcm(1),xsize_ibm_fcm(2),xsize_ibm_fcm(3),diff1_coef_ibm)
            ! call D1s_IBM_O2_ACC_3Dx(p12_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_x_ibm,n1s,n1e,n2s,n2e,n3s,n3e,conv1_coef_ibm)
            call D1c_IBM_O2_ACC_3Dx(p12_x_q2_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp2_x_ibm,xsize_ibm_fcm(1),xsize_ibm_fcm(2),xsize_ibm_fcm(3),conv1_coef_ibm)

            ! For RHS3
            call D2c_IBM_O2_3Dx(q3_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp3_x_ibm,xsize_ibm_fcm(1),xsize_ibm_fcm(2),xsize_ibm_fcm(3),diff1_coef_ibm)
            ! call D1s_IBM_O2_ACC_3Dx(p13_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp3_x_ibm,n1s,n1e,n2s,n2e,n3s,n3e,conv1_coef_ibm)
            call D1c_IBM_O2_ACC_3Dx(p13_x_q3_ibm(n1s:n1e,n2s:n2e,n3s:n3e),tmp3_x_ibm,xsize_ibm_fcm(1),xsize_ibm_fcm(2),xsize_ibm_fcm(3),conv1_coef_ibm)
        

            if((BC1==NOSLIP).or.(BC1==PSEUDO_PERIODIC).or.(BC1==OPEN)) call diff_at_wall_1_ibm   ! ATTENTION

            do i=xstart_ibm_fcm(1),xend_ibm_fcm(1)
                do j=xstart_ibm_fcm(2),xend_ibm_fcm(2)
                    do k=xstart_ibm_fcm(3),xend_ibm_fcm(3)

                        RHS1_x_ibm(i,j,k) = RHS1_x_ibm(i,j,k) + tmp1_x_ibm(i,j,k)
                        RHS2_x_ibm(i,j,k) = RHS2_x_ibm(i,j,k) + tmp2_x_ibm(i,j,k)
                        RHS3_x_ibm(i,j,k) = RHS3_x_ibm(i,j,k) + tmp3_x_ibm(i,j,k)

                    enddo
                enddo
            enddo

        endif

        t_aft = MPI_WTIME()
        ! write(*,*) 'finish dx', t_aft-t_bef

        return

    contains

        subroutine diff_at_wall_1_ibm()

            ! d2f=a1*f0 + a2*f1 + a3*f2
            ! Regular case: down    a1=2/3, a2=-1 a3=1/3
            ! Regular case: up      a1=1/3, a2=-1 a3=2/3
            ! ATTENTION
            if ((xstart_ibm_fcm(1).eq.1).and.(object_in_current_proc_fcm_x)) then

                do k=xstart_ibm_fcm(3), min(n3_ibm-1, xend_ibm_fcm(3))       !do k=1,n3m
                    do j=xstart_ibm_fcm(2), min(n2_ibm-1, xend_ibm_fcm(2))       !do j=1,n2m
                        
                            tmp2_x_ibm(1,j,k)   =tmp2_x_ibm(1,j,k)  +( q2_x_ibm(2,j,k)       /3.d0 - q2_x_ibm(1,j,k)      + q2_wall10_ibm(j,k)*2.d0/3.d0)  /diff1_coef_ibm**2

                            tmp3_x_ibm(1,j,k)   =tmp3_x_ibm(1,j,k)  +( q3_x_ibm(2,j,k)       /3.d0 - q3_x_ibm(1,j,k)      + q3_wall10_ibm(j,k)*2.d0/3.d0 ) /diff1_coef_ibm**2

                    enddo
                enddo

            endif

            if ((xend_ibm_fcm(1).eq.(n1_ibm-1)).and.(object_in_current_proc_fcm_x)) then

                do k=xstart_ibm_fcm(3), min(n3_ibm-1, xend_ibm_fcm(3))       !do k=1,n3m
                    do j=xstart_ibm_fcm(2), min(n2_ibm-1, xend_ibm_fcm(2))       !do j=1,n2m
        
                            tmp2_x_ibm(n1_ibm-1,j,k) =tmp2_x_ibm(n1_ibm-1,j,k)+( q2_x_ibm(n1_ibm-2,j,k)   /3.d0 - q2_x_ibm(n1_ibm-1,j,k)    + q2_wall11_ibm(j,k)*2.d0/3.d0 ) /diff1_coef_ibm**2

                            tmp3_x_ibm(n1_ibm-1,j,k) =tmp3_x_ibm(n1_ibm-1,j,k)+( q3_x_ibm(n1_ibm-2,j,k)   /3.d0 - q3_x_ibm(n1_ibm-1,j,k)    + q3_wall11_ibm(j,k)*2.d0/3.d0)  /diff1_coef_ibm**2

                    enddo
                enddo

            endif

        end subroutine diff_at_wall_1_ibm

        subroutine diff_at_wall_2_ibm()

            if ((ystart_ibm_fcm(2).eq.1).and.(object_in_current_proc_fcm_y)) then

                do k=ystart_ibm_fcm(3), min(n3_ibm-1, yend_ibm_fcm(3))       !do k=1,n3m
                    do i=ystart_ibm_fcm(1), min(n1_ibm-1, yend_ibm_fcm(1))       !do i=1,n1m

                            tmp1_y_ibm(i,1,k)=tmp1_y_ibm(i,1,k)+( q1_y_ibm(i,2,k)*a3_d + q1_y_ibm(i,1,k)*a2_d + q1_wall20_ibm(i,k)*a1_d)            /diff2_coef_ibm**2

                            tmp3_y_ibm(i,1,k)=tmp3_y_ibm(i,1,k)+( q3_y_ibm(i,2,k)*a3_d + q3_y_ibm(i,1,k)*a2_d + q3_wall20_ibm(i,k)*a1_d )           /diff2_coef_ibm**2
                
                    enddo
                enddo

            endif

            if ((yend_ibm_fcm(2).eq.(n2_ibm-1)).and.(object_in_current_proc_fcm_y)) then

                do k=ystart_ibm_fcm(3), min(n3_ibm-1, yend_ibm_fcm(3))       !do k=1,n3m
                    do i=ystart_ibm_fcm(1), min(n1_ibm-1, yend_ibm_fcm(1))       !do i=1,n1m

                            tmp1_y_ibm(i,n2_ibm-1,k)=tmp1_y_ibm(i,n2_ibm-1,k)+( q1_y_ibm(i,n2m-1,k)*a1_u + q1_y_ibm(i,n2m,k)*a2_u + q1_wall21_ibm(i,k)*a3_u ) /diff2_coef_ibm**2

                            tmp3_y_ibm(i,n2_ibm-1,k)=tmp3_y_ibm(i,n2_ibm-1,k)+( q3_y_ibm(i,n2_ibm-2,k)*a1_u + q3_y_ibm(i,n2_ibm-1,k) *a2_u + q3_wall21_ibm(i,k)*a3_u) /diff2_coef_ibm**2
                
                    enddo
                enddo

            endif

        end subroutine diff_at_wall_2_ibm

        subroutine diff_at_wall_3_ibm()

            ! d2f=a1*f0 + a2*f1 + a3*f2
            ! Regular case: down    a1=2/3, a2=-1 a3=1/3
            ! Regular case: up      a1=1/3, a2=-1 a3=2/3
            ! ATTENTION

            ! ************* CORR_JO ******************* !
            ! n3m au lieu de n2m dans RHS_z(i,j,n3m) ** !
            ! ***************************************** !

            if ((zstart_ibm_fcm(3).eq.1).and.(object_in_current_proc_fcm_z)) then

                do j = zstart_ibm_fcm(2), min(n2_ibm-1, zend_ibm_fcm(2))     !do j=1,n2m
                    do i=zstart_ibm_fcm(1), min(n1_ibm-1, zend_ibm_fcm(1))       !do i=1,n1m
                        
                            tmp1_z_ibm(i,j,1)   =tmp1_z_ibm(i,j,1)  +( q1_z_ibm(i,j,2)       /3.d0 - q1_z_ibm(i,j,1)      + q1_wall30_ibm(i,j)*2.d0/3.d0 ) /diff3_coef_ibm**2

                            tmp2_z_ibm(i,j,1)   =tmp2_z_ibm(i,j,1)  +( q2_z_ibm(i,j,2)       /3.d0 - q2_z_ibm(i,j,1)      + q2_wall30_ibm(i,j)*2.d0/3.d0)  /diff3_coef_ibm**2
                    
                    enddo
                enddo

            endif

            if ((zend_ibm_fcm(3).eq.(n3_ibm-1)).and.(object_in_current_proc_fcm_z)) then

                do j = zstart_ibm_fcm(2), min(n2_ibm-1, zend_ibm_fcm(2))     !do j=1,n2m
                    do i=zstart_ibm_fcm(1), min(n1_ibm-1, zend_ibm_fcm(1))       !do i=1,n1m
                        
                            tmp1_z_ibm(i,j,n3_ibm-1) =tmp1_z_ibm(i,j,n3_ibm-1)+( q1_z_ibm(i,j,n3_ibm-2)   /3.d0 - q1_z_ibm(i,j,n3_ibm-1)    + q1_wall31_ibm(i,j)*2.d0/3.d0)  /diff3_coef_ibm**2

                            tmp2_z_ibm(i,j,n3_ibm-1) =tmp2_z_ibm(i,j,n3m)+( q2_z_ibm(i,j,n3_ibm-2)   /3.d0 - q2_z_ibm(i,j,n3_ibm-1)    + q2_wall31_ibm(i,j)*2.d0/3.d0 ) /diff3_coef_ibm**2
                    
                    enddo
                enddo

            endif

        end subroutine diff_at_wall_3_ibm

    end subroutine


    subroutine perform_intermediate_velocity()

        use IBM_settings
        use IBM_data
        use IBM

        implicit none

        real*8                                                                                  :: RHS1i, RHS2i, RHS3i
        integer                                                                                 :: i,j,k
        integer                                                                                 :: n1s, n1e, n2s,n2e, n3s,n3e

        n1e=min(n1-1, xend(1))
        n2e=min(n2-1, xend(2))
        n3e=min(n3-1, xend(3))

        if (IBM_activated) then

            if ((interpol_type == IBM_INTERPOL_NONE).or.(interpol_type == SECOND_ORDER_INTERPOL)) then
                call force_velocity(q1_save_x, q2_save_x, q3_save_x, xsize, vel_term1, vel_term2, vel_term3)
            else if (interpol_type == ANTISYMMETRIC_INTERPOL) then
                call compute_antisymmetric_velocity(xsize, vel_term1, vel_term2, vel_term3)
            endif

            if (interpol_type == ANTISYMMETRIC_INTERPOL) then

                if ((BC1==UNBOUNDED).or.(BC1==FRINGE)) n1s=xstart(1)
                if ((BC1/=UNBOUNDED).and.(BC1/=FRINGE)) n1s=max(2,xstart(1))
                n2s=xstart(2)
                n3s=xstart(3)

                do k = n3s, n3e
                    do j = n2s, n2e
                        do i = n1s, n1e
                            RHS1i = NSTERMS1(i,j,k)
                            q1_x(i,j,k) = q1_save_x(i,j,k) + RHS1i + mask_kim_q1_x(i,j,k)*(-RHS1i + vel_term1(i,j,k)) 
                        end do
                    end do
                end do


                if (BC2==UNBOUNDED) n2s=xstart(2)
                if (BC2/=UNBOUNDED) n2s=max(2,xstart(2))
                ! On ne touche que q2[2..n1-2]
                if (BC1==OPEN) then
                    n1s=2
                    n1e=n1-2
                endif

                n3s=xstart(3)

                do k = n3s, n3e
                    do j = n2s, n2e
                        do i = n1s, n1e
                            RHS2i = NSTERMS2(i,j,k)
                            q2_x(i,j,k) = q2_save_x(i,j,k) + RHS2i + mask_kim_q2_x(i,j,k)*(-RHS2i + vel_term2(i,j,k)) 
                        end do
                    end do
                end do

                if ((BC3==UNBOUNDED).or.(BC3==FRINGE)) n3s=xstart(3)
                if ((BC3/=UNBOUNDED).and.(BC3/=FRINGE)) n3s=max(2,xstart(3))

                ! On ne touche que q3[2..n1-2]
                if (BC1==OPEN) then
                    n1s=2
                    n1e=n1-2
                endif

                n2s=xstart(2)

                do k = n3s, n3e
                    do j = n2s, n2e
                        do i = n1s, n1e
                            RHS3i = NSTERMS3(i,j,k)
                            q3_x(i,j,k) = q3_save_x(i,j,k) + RHS3i + mask_kim_q3_x(i,j,k)*(-RHS3i + vel_term3(i,j,k)) 
                        end do
                    end do
                end do

            else
                ! By default, set velocity to 0 in the objects ...

                if ((BC1==UNBOUNDED).or.(BC1==FRINGE)) n1s=xstart(1)
                if ((BC1/=UNBOUNDED).and.(BC1/=FRINGE)) n1s=max(2,xstart(1))
                n2s=xstart(2)
                n3s=xstart(3)

                do k = n3s, n3e
                    do j = n2s, n2e
                        do i = n1s, n1e
                            RHS1i = NSTERMS1(i,j,k)
                            q1_x(i,j,k) = q1_save_x(i,j,k) + RHS1i + mask_objects_q1_x(i,j,k)*(-RHS1i + vel_term1(i,j,k)) 
                        end do
                    end do
                end do


                if (BC2==UNBOUNDED) n2s=xstart(2)
                if (BC2/=UNBOUNDED) n2s=max(2,xstart(2))
                ! On ne touche que q2[2..n1-2]
                if (BC1==OPEN) then
                    n1s=2
                    n1e=n1-2
                endif

                n3s=xstart(3)

                do k = n3s, n3e
                    do j = n2s, n2e
                        do i = n1s, n1e
                            RHS2i = NSTERMS2(i,j,k)
                            q2_x(i,j,k) = q2_save_x(i,j,k) + RHS2i + mask_objects_q2_x(i,j,k)*(-RHS2i + vel_term2(i,j,k)) 
                        end do
                    end do
                end do

                if ((BC3==UNBOUNDED).or.(BC3==FRINGE)) n3s=xstart(3)
                if ((BC3/=UNBOUNDED).and.(BC3/=FRINGE)) n3s=max(2,xstart(3))

                ! On ne touche que q3[2..n1-2]
                if (BC1==OPEN) then
                    n1s=2
                    n1e=n1-2
                endif

                n2s=xstart(2)

                do k = n3s, n3e
                    do j = n2s, n2e
                        do i = n1s, n1e
                            RHS3i = NSTERMS3(i,j,k)
                            q3_x(i,j,k) = q3_save_x(i,j,k) + RHS3i + mask_objects_q3_x(i,j,k)*(-RHS3i + vel_term3(i,j,k)) 
                        end do
                    end do
                end do

            endif


        else if (MHD_state.eq.0) then ! MHD is not used or is activated (either active or passive)

            if ((BC1==UNBOUNDED).or.(BC1==FRINGE)) n1s=xstart(1)
            if ((BC1/=UNBOUNDED).and.(BC1/=FRINGE)) n1s=max(2,xstart(1))
            n2s=xstart(2)
            n3s=xstart(3)
            do k = n3s, n3e
                do j = n2s, n2e
                    do i = n1s, n1e
                        q1_x(i,j,k) = q1_save_x(i,j,k) + NSTERMS1(i,j,k)
                    end do
                end do
            end do


            if (BC2==UNBOUNDED) n2s=xstart(2)
            if (BC2/=UNBOUNDED) n2s=max(2,xstart(2))
            ! On ne touche que q2[2..n1-2]
            if (BC1==OPEN) then
                n1s=2
                n1e=n1-2
            endif

            n3s=xstart(3)

            do k = n3s, n3e
                do j = n2s, n2e
                    do i = n1s, n1e
                        q2_x(i,j,k) = q2_save_x(i,j,k) + NSTERMS2(i,j,k)
                    end do
                end do
            end do

            if (BC3==UNBOUNDED) n3s=xstart(3)
            if (BC3/=UNBOUNDED) n3s=max(2,xstart(3))

            ! On ne touche que q3[2..n1-2]
            if (BC1==OPEN) then
                n1s=2
                n1e=n1-2
            endif

            n2s=xstart(2)

            do k = n3s, n3e
                do j = n2s, n2e
                    do i = n1s, n1e
                        q3_x(i,j,k) = q3_save_x(i,j,k) + NSTERMS3(i,j,k)
                    end do
                end do
            end do

        end if
        ! q3_x, q2_x and w_c contains u*,v* and w*

    end subroutine perform_intermediate_velocity

    subroutine perform_intermediate_velocity_ibm()

        use IBM_settings
        use IBM_data
        use IBM

        implicit none

        real*8                                                                                  :: RHS1i_ibm, RHS2i_ibm, RHS3i_ibm
        integer                                                                                 :: i,j,k
        integer                                                                                 :: n1s, n1e, n2s,n2e, n3s,n3e

        ! Old routine used for refined mesh for IBM

        ! n1e=min(n1_ibm-1, xend_ibm_fpm(1))
        ! n2e=min(n2_ibm-1, xend_ibm_fpm(2))
        ! n3e=min(n3_ibm-1, xend_ibm_fpm(3))

        ! call force_velocity(q1_save2_ibm, q2_save2_ibm, q3_save2_ibm, xsize_fine, vel_term1_ibm, vel_term2_ibm, vel_term3_ibm)

        ! if ((BC1==UNBOUNDED).or.(BC1==FRINGE)) n1s=xstart_ibm_fpm(1)
        ! if ((BC1/=UNBOUNDED).and.(BC1/=FRINGE)) n1s=max(2,xstart_ibm_fpm(1))
        ! n2s=xstart_ibm_fpm(2)
        ! n3s=xstart_ibm_fpm(3)

        ! if (object_in_current_proc_fpm_x) then

        !     do k = n3s, n3e
        !         do j = n2s, n2e
        !             do i = n1s, n1e

        !                 RHS1i_ibm = NSTERMS1_ibm(i,j,k) + Fext1_ibm(i, j, k)
        !                 q1_x_ibm(i,j,k) = q1_save2_ibm(i,j,k) + RHS1i_ibm + IBM_mask1_x_ibm(i,j,k)*(-RHS1i_ibm + vel_term1_ibm(i,j,k))

        !                 ! q1_x_ibm(i,j,k) = q1_save2_ibm(i,j,k) + NSTERMS1_ibm(i,j,k) + Fext1_ibm(i,j,k)

        !             end do
        !         end do
        !     end do

        ! endif

        ! if (BC2==UNBOUNDED) n2s=xstart_ibm_fpm(2)
        ! if (BC2/=UNBOUNDED) n2s=max(2,xstart_ibm_fpm(2))
        ! ! On ne touche que q2[2..n1-2]
        ! if (BC1==OPEN) then
        !     n1s=2
        !     n1e=n1_ibm-2
        ! endif

        ! n3s=xstart_ibm_fpm(3)

        ! if (object_in_current_proc_fpm_x) then

        !     do k = n3s, n3e
        !         do j = n2s, n2e
        !             do i = n1s, n1e

        !                 RHS2i_ibm = NSTERMS2_ibm(i,j,k) + Fext2_ibm(i, j, k)
        !                 q2_x_ibm(i,j,k) = q2_save2_ibm(i,j,k) + RHS2i_ibm + IBM_mask2_x_ibm(i,j,k)*(-RHS2i_ibm + vel_term2_ibm(i,j,k))

        !                 ! q2_x_ibm(i,j,k) = q2_save2_ibm(i,j,k) + NSTERMS2_ibm(i,j,k) + Fext2_ibm(i,j,k)

        !             end do
        !         end do
        !     end do

        ! endif

        ! if ((BC3==UNBOUNDED).or.(BC3==FRINGE)) n3s=xstart_ibm_fpm(3)
        ! if ((BC3/=UNBOUNDED).and.(BC3/=FRINGE)) n3s=max(2,xstart_ibm_fpm(3))

        ! ! On ne touche que q3[2..n1-2]
        ! if (BC1==OPEN) then
        !     n1s=2
        !     n1e=n1_ibm-2
        ! endif

        ! n2s=xstart_ibm_fpm(2)

        ! if (object_in_current_proc_fpm_x) then

        !     do k = n3s, n3e
        !         do j = n2s, n2e
        !             do i = n1s, n1e

        !                 RHS3i_ibm = NSTERMS3_ibm(i,j,k) + Fext3_ibm(i, j, k)
        !                 q3_x_ibm(i,j,k) = q3_save2_ibm(i,j,k) + RHS3i_ibm + IBM_mask3_x_ibm(i,j,k)*(-RHS3i_ibm + vel_term3_ibm(i,j,k))

        !                 ! q3_x_ibm(i,j,k) = q3_save2_ibm(i,j,k) + NSTERMS3_ibm(i,j,k) + Fext3_ibm(i,j,k)

        !             end do
        !         end do
        !     end do

        ! endif

    end subroutine perform_intermediate_velocity_ibm

    ! Perform UV, UW, VW at the center of edges
    subroutine perform_velocity_products()
        use mesh
        use IBM
        use DRP_IBM

        implicit none

        integer k,j,i
        real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3))   :: tmp12_y
        real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3))   :: tmp23_z, tmp13_z

        real*8, dimension(:,:,:), allocatable  :: tmp12_x_ibm, tmp13_x_ibm
        real*8, dimension(:,:,:), allocatable  :: tmp12_y_ibm, tmp23_y_ibm
        real*8, dimension(:,:,:), allocatable  :: tmp23_z_ibm, tmp13_z_ibm

        integer :: n1e, n2e, n3e, n1s, n2s, n3s

        ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) then

        !     allocate(tmp12_x_ibm(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
        !     allocate(tmp13_x_ibm(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))

        !     tmp12_x_ibm = 0.d0
        !     tmp13_x_ibm = 0.d0

        !     allocate(tmp12_y_ibm(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
        !     allocate(tmp23_y_ibm(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
            
        !     tmp12_y_ibm = 0.d0
        !     tmp23_y_ibm = 0.d0

        !     allocate(tmp23_z_ibm(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
        !     allocate(tmp13_z_ibm(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
            
        !     tmp23_z_ibm = 0.d0
        !     tmp13_z_ibm = 0.d0

        ! endif

        ! X ------------------------------------------------------------------------

        call D0ssh_3Dx(q3_x, p13_x, n1,xsize(2),xsize(2),xsize(3),xsize(3), NS_Q3_BC1)
        call D0ssh_3Dx(q2_x, p12_x, n1,xsize(2),xsize(2),xsize(3),xsize(3), NS_Q2_BC1)

        ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) then

        !     if (object_in_current_proc_ccm_x) then

        !         n1s=xstart_ibm_ccm(1)
        !         n2s=xstart_ibm_ccm(2)
        !         n3s=xstart_ibm_ccm(3)

        !         n1e=xend_ibm_ccm(1)
        !         n2e=xend_ibm_ccm(2)
        !         n3e=xend_ibm_ccm(3)

        !         ! for p13_ibm
        !         call D0ssh_IBM_O2_3Dx(q3_x(n1s:n1e,n2s:n2e,n3s:n3e), tmp13_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e),xsize_ibm_ccm(1),xsize_ibm_ccm(2),xsize_ibm_ccm(3))
                
        !         ! for p12_ibm
        !         call D0ssh_IBM_O2_3Dx(q2_x(n1s:n1e,n2s:n2e,n3s:n3e), tmp12_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e),xsize_ibm_ccm(1),xsize_ibm_ccm(2),xsize_ibm_ccm(3))

        !     endif

        !     if (object_in_current_proc_cpm_x) then

        !         do i=xstart_ibm_cpm(1), xend_ibm_cpm(1)
        !             do j=xstart_ibm_cpm(2), xend_ibm_cpm(2)
        !                 do k=xstart_ibm_cpm(3), xend_ibm_cpm(3)

        !                     p13_x(i,j,k) = tmp13_x_ibm(i,j,k)

        !                     p12_x(i,j,k) = tmp12_x_ibm(i,j,k)

        !                 enddo
        !             enddo
        !         enddo

        !     endif

        ! endif

        if ((BC1==NOSLIP).or.(BC1==FREESLIP).or.(BC1==PSEUDO_PERIODIC).or.(BC1==OPEN)) then

            do k=xstart(3), xend(3)       !do k=1,n3
                do j=xstart(2), xend(2)       !do j=1,n2
                    p13_x(1, j,  k)   =q3_wall10(j, k)
                    p13_x(n1, j, k)   =q3_wall11(j, k)

                    p12_x(1, j,  k)   =q2_wall10(j, k)
                    p12_x(n1, j, k)   =q2_wall11(j, k)
                enddo
            enddo

        end if

        ! Y ------------------------------------------------------------------------
        call transpose_x_to_y(p12_x, p12_y)

        call D0ssh_3Dy(q3_y, p23_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_Q3_BC2)
        ! call D0ssh_MULT_3Dy(q1_y, p12_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_Q1_BC2)
        call D0ssh_3Dy(q1_y, tmp12_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_Q1_BC2)

        ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) then

        !     if (object_in_current_proc_ccm_y) then

        !         n1s=ystart_ibm_ccm(1)
        !         n2s=ystart_ibm_ccm(2)
        !         n3s=ystart_ibm_ccm(3)

        !         n1e=yend_ibm_ccm(1)
        !         n2e=yend_ibm_ccm(2)
        !         n3e=yend_ibm_ccm(3)

        !         ! for p13_ibm
        !         call D0ssh_IBM_O2_3Dy(q3_y(n1s:n1e,n2s:n2e,n3s:n3e), tmp23_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),ysize_ibm_ccm(1),ysize_ibm_ccm(2),ysize_ibm_ccm(3))
                
        !         ! for p12_ibm
        !         call D0ssh_IBM_O2_3Dy(q1_y(n1s:n1e,n2s:n2e,n3s:n3e), tmp12_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e),ysize_ibm_ccm(1),ysize_ibm_ccm(2),ysize_ibm_ccm(3))

        !     endif

        !     if (object_in_current_proc_cpm_y) then

        !         do i=ystart_ibm_cpm(1), yend_ibm_cpm(1)
        !             do j=ystart_ibm_cpm(2), yend_ibm_cpm(2)
        !                 do k=ystart_ibm_cpm(3), yend_ibm_cpm(3)

        !                     p23_y(i,j,k) = tmp23_y_ibm(i,j,k)

        !                     tmp12_y(i,j,k) = tmp12_y_ibm(i,j,k)

        !                 enddo
        !             enddo
        !         enddo

        !     endif

        ! endif

        if ((BC2==NOSLIP).or.(BC2==FREESLIP)) then

            do k=ystart(3), yend(3)      !do k=1,n3
                do i=ystart(1),yend(1)       !do i=1,n1
                    p23_y(i, 1,  k)   =q3_wall20(i, k)
                    p23_y(i, n2, k)   =q3_wall21(i, k)

                    ! p12_y(i, 1,  k)   =p12_y(i, 1,  k) * q1_wall20(i, k)
                    ! p12_y(i, n2, k)   =p12_y(i, n2, k) * q1_wall21(i, k)

                    tmp12_y(i, 1,  k)   = q1_wall20(i, k)
                    tmp12_y(i, n2, k)   = q1_wall21(i, k)
                enddo
            enddo

        end if

        p12_y = p12_y * tmp12_y

        ! Z ------------------------------------------------------------------------
        call transpose_y_to_z(p23_y, p23_z)

        call transpose_x_to_y(p13_x, p13_y)
        call transpose_y_to_z(p13_y, p13_z)

        ! call D0ssh_MULT_3Dz(q2_z, p23_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_Q2_BC3)
        ! call D0ssh_MULT_3Dz(q1_z, p13_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_Q1_BC3)

        call D0ssh_3Dz(q2_z, tmp23_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_Q2_BC3)
        call D0ssh_3Dz(q1_z, tmp13_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_Q1_BC3)

        ! if ((IBM_activated).and.(interpol_type == SECOND_ORDER_INTERPOL)) then

        !     if (object_in_current_proc_ccm_z) then

        !         n1s=zstart_ibm_ccm(1)
        !         n2s=zstart_ibm_ccm(2)
        !         n3s=zstart_ibm_ccm(3)

        !         n1e=zend_ibm_ccm(1)
        !         n2e=zend_ibm_ccm(2)
        !         n3e=zend_ibm_ccm(3)

        !         ! for p13_ibm
        !         call D0ssh_IBM_O2_3Dz(q2_z(n1s:n1e,n2s:n2e,n3s:n3e), tmp23_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e),zsize_ibm_ccm(1),zsize_ibm_ccm(2),zsize_ibm_ccm(3))
                
        !         ! for p12_ibm
        !         call D0ssh_IBM_O2_3Dz(q1_z(n1s:n1e,n2s:n2e,n3s:n3e), tmp13_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e),zsize_ibm_ccm(1),zsize_ibm_ccm(2),zsize_ibm_ccm(3))

        !     endif

        !     if (object_in_current_proc_cpm_z) then

        !         do i=zstart_ibm_cpm(1), zend_ibm_cpm(1)
        !             do j=zstart_ibm_cpm(2), zend_ibm_cpm(2)
        !                 do k=zstart_ibm_cpm(3), zend_ibm_cpm(3)

        !                     tmp23_z(i,j,k) = tmp23_z_ibm(i,j,k)

        !                     tmp13_z(i,j,k) = tmp13_z_ibm(i,j,k)

        !                 enddo
        !             enddo
        !         enddo

        !     endif

        ! endif

        ! ATTENTION
        if ((BC3==NOSLIP).or.(BC3==FREESLIP)) then

            do i=zstart(1), zend(1)
                do j=zstart(2), zend(2)
                    ! p23_z(i, j,  1)   =p23_z(i, j, 1)   * q2_wall30(i, j)
                    ! p23_z(i, j, n3)   =p23_z(i, j, n3)  * q2_wall31(i, j)

                    ! p13_z(i, j,  1)   =p13_z(i, j, 1)   * q1_wall30(i, j)
                    ! p13_z(i, j, n3)   =p13_z(i, j, n3)  * q1_wall31(i, j)

                    tmp23_z(i, j,  1)   = q2_wall30(i, j)
                    tmp23_z(i, j, n3)   = q2_wall31(i, j)

                    tmp13_z(i, j,  1)   = q1_wall30(i, j)
                    tmp13_z(i, j, n3)   = q1_wall31(i, j)
                end do
            enddo

        end if
        ! ENDATTENTION

        p23_z = p23_z * tmp23_z
        p13_z = p13_z * tmp13_z

        return

    end subroutine

    ! Adaptation of perform_velocity_products, but on the finest grid while the IBM with second order box is used
    ! Perform UV, UW, VW at the center of edges
    subroutine perform_velocity_products_ibm()
        use mesh
        use IBM
        use IBM_data
        use DRP_IBM

        implicit none

        integer k,j,i

        real*8, dimension(xstart_fine(1):xend_fine(1), xstart_fine(2):xend_fine(2), xstart_fine(3):xend_fine(3)) :: upshift_q3_x, upshift_q2_x, upshift_q1_x
        real*8, dimension(ystart_fine(1):yend_fine(1), ystart_fine(2):yend_fine(2), ystart_fine(3):yend_fine(3)) :: upshift_q3_y, upshift_q2_y, upshift_q1_y
        real*8, dimension(zstart_fine(1):zend_fine(1), zstart_fine(2):zend_fine(2), zstart_fine(3):zend_fine(3)) :: upshift_q3_z, upshift_q2_z, upshift_q1_z

        real*8, dimension(ystart_fine(1):yend_fine(1), ystart_fine(2):yend_fine(2), ystart_fine(3):yend_fine(3)) :: tmp_y

        integer :: n1e, n2e, n3e, n1s, n2s, n3s

        ! X ------------------------------------------------------------------------

        if (object_in_current_proc_fcm_z) then

            n1s=zstart_ibm_fcm(1)
            n2s=zstart_ibm_fcm(2)
            n3s=zstart_ibm_fcm(3)

            n1e=zend_ibm_fcm(1)
            n2e=zend_ibm_fcm(2)
            n3e=zend_ibm_fcm(3)

            call Forward_shift_IBM_3Dz(q3_z_ibm(n1s:n1e,n2s:n2e,n3s:n3e), upshift_q3_z,zsize_ibm_fcm(1),zsize_ibm_fcm(2),zsize_ibm_fcm(3))
        endif

        call transpose_z_to_y(upshift_q3_z,upshift_q3_y, decomp_fine)
        call transpose_y_to_x(upshift_q3_y,upshift_q3_x, decomp_fine)
            
        if (object_in_current_proc_fcm_y) then

            n1s=ystart_ibm_fcm(1)
            n2s=ystart_ibm_fcm(2)
            n3s=ystart_ibm_fcm(3)

            n1e=yend_ibm_fcm(1)
            n2e=yend_ibm_fcm(2)
            n3e=yend_ibm_fcm(3)

            call Forward_shift_IBM_3Dy(q2_y_ibm(n1s:n1e,n2s:n2e,n3s:n3e), upshift_q2_y,ysize_ibm_fcm(1),ysize_ibm_fcm(2),ysize_ibm_fcm(3))
        endif

        call transpose_y_to_x(upshift_q2_y,upshift_q2_x, decomp_fine)

        if (object_in_current_proc_fcm_x) then
            do i=max(i_start_ibm_2nd_fc+1,xstart_ibm_fcm(1)), xend_ibm_fcm(1)
                do j=xstart_ibm_fcm(2), xend_ibm_fcm(2)
                    do k=xstart_ibm_fcm(3), xend_ibm_fcm(3)

                        p13_x_q1_ibm(i,j,k) = upshift_q3_x(i-1,j,k)
                        p13_x_q3_ibm(i,j,k) = q3_x_ibm(i,j,k)

                        p12_x_q1_ibm(i,j,k) = upshift_q2_x(i-1,j,k)
                        p12_x_q2_ibm(i,j,k) = q2_x_ibm(i,j,k)

                    enddo
                enddo
            enddo

        endif

        call transpose_x_to_y(p13_x_q1_ibm, tmp_y, decomp_fine)
        call transpose_y_to_z(tmp_y, p13_z_q1_ibm, decomp_fine)

        call transpose_x_to_y(p13_x_q3_ibm, tmp_y, decomp_fine)
        call transpose_y_to_z(tmp_y, p13_z_q3_ibm, decomp_fine)

        call transpose_x_to_y(p12_x_q1_ibm, p12_y_q1_ibm, decomp_fine)
        call transpose_x_to_y(p12_x_q2_ibm, p12_y_q2_ibm, decomp_fine)

        ! Y ------------------------------------------------------------------------

        if (object_in_current_proc_fcm_x) then

            n1s=xstart_ibm_fcm(1)
            n2s=xstart_ibm_fcm(2)
            n3s=xstart_ibm_fcm(3)

            n1e=xend_ibm_fcm(1)
            n2e=xend_ibm_fcm(2)
            n3e=xend_ibm_fcm(3)

            call Forward_shift_IBM_3Dx(q1_x_ibm(n1s:n1e,n2s:n2e,n3s:n3e), upshift_q1_x,xsize_ibm_fcm(1),xsize_ibm_fcm(2),xsize_ibm_fcm(3))
        endif

        call transpose_x_to_y(upshift_q1_x,upshift_q1_y, decomp_fine)

            ! if ((BC2==NOSLIP).or.(BC2==FREESLIP)) then

            !     if (ystart_ibm_fcm(2).eq.1) then

            !         do k=ystart_ibm_fcm(3), yend_ibm_fcm(3)      !do k=1,n3
            !             do i=ystart_ibm_fcm(1),yend_ibm_fcm(1)       !do i=1,n1

            !                 tmp23_y_ibm(i, 1,  k)   =q3_wall20_ibm(i, k)

            !                 tmp12_y_ibm(i, 1,  k)   = q1_wall20_ibm(i, k)

            !             enddo
            !         enddo

            !     endif

            !     if (yend_ibm_fcm(2).eq.n2_ibm) then

            !         do k=ystart_ibm_fcm(3), yend_ibm_fcm(3)      !do k=1,n3
            !             do i=ystart_ibm_fcm(1),yend_ibm_fcm(1)       !do i=1,n1

            !                 tmp23_y_ibm(i, n2_ibm, k)   =q3_wall21_ibm(i, k)

            !                 tmp12_y_ibm(i, n2_ibm, k)   = q1_wall21_ibm(i, k)

            !             enddo
            !         enddo

            !     endif

            ! end if

            ! do i=ystart_ibm_fcm(1), yend_ibm_fcm(1)
            !     do j=ystart_ibm_fcm(2), yend_ibm_fcm(2)
            !         do k=ystart_ibm_fcm(3), yend_ibm_fcm(3)

            !             p23_y_ibm(i,j,k) = tmp23_y_ibm(i,j,k)

            !             p12_y_ibm(i,j,k) = p12_y_ibm(i,j,k) * tmp12_y_ibm(i,j,k)

            !         enddo
            !     enddo
            ! enddo

        if (object_in_current_proc_fcm_y) then
            do i=ystart_ibm_fcm(1), yend_ibm_fcm(1)
                do j=max(j_start_ibm_2nd_fc+1,ystart_ibm_fcm(2)), yend_ibm_fcm(2)
                    do k=ystart_ibm_fcm(3), yend_ibm_fcm(3)

                        p23_y_q2_ibm(i,j,k) = upshift_q3_y(i,j-1,k)
                        p23_y_q3_ibm(i,j,k) = q3_y_ibm(i,j,k)

                        p12_y_q1_ibm(i,j,k) = p12_y_q1_ibm(i,j,k) * q1_y_ibm(i,j,k)
                        p12_y_q2_ibm(i,j,k) = p12_y_q2_ibm(i,j,k) * upshift_q1_y(i,j-1,k)

                    enddo
                enddo
            enddo

        endif

        call transpose_y_to_z(p23_y_q2_ibm, p23_z_q2_ibm, decomp_fine)
        call transpose_y_to_z(p23_y_q3_ibm, p23_z_q3_ibm, decomp_fine)

        call transpose_y_to_z(upshift_q1_y,upshift_q1_z, decomp_fine)
        call transpose_y_to_z(upshift_q2_y,upshift_q2_z, decomp_fine)

        ! Z ------------------------------------------------------------------------

        if (object_in_current_proc_fcm_z) then

            do i=zstart_ibm_fcm(1), zend_ibm_fcm(1)
                do j=zstart_ibm_fcm(2), zend_ibm_fcm(2)
                    do k=max(k_start_ibm_2nd_fc+1,zstart_ibm_fcm(3)), zend_ibm_fcm(3)
                        
                        p23_z_q2_ibm(i,j,k) = p23_z_q2_ibm(i,j,k) * q2_z_ibm(i,j,k)
                        p23_z_q3_ibm(i,j,k) = p23_z_q3_ibm(i,j,k) * upshift_q2_z(i,j,k-1)

                        p13_z_q1_ibm(i,j,k) = p13_z_q1_ibm(i,j,k) * q1_z_ibm(i,j,k)
                        p13_z_q3_ibm(i,j,k) = p13_z_q3_ibm(i,j,k) * upshift_q1_z(i,j,k-1)

                    enddo
                enddo
            enddo

        endif

        call transpose_z_to_y(p23_z_q3_ibm, p23_y_q3_ibm, decomp_fine)

        call transpose_y_to_x(p12_y_q2_ibm, p12_x_q2_ibm, decomp_fine)

        call transpose_z_to_y(p13_z_q3_ibm, tmp_y, decomp_fine)
        call transpose_y_to_x(tmp_y, p13_x_q3_ibm, decomp_fine)

        return

    end subroutine

    subroutine save_state()

        use VELOCITY_workspace_view, only: recovery_RHS_dir
        use HDF5_IO

        implicit none

        character(200)    :: file_path

        file_path=trim(recovery_RHS_dir)//"/RHS1"
        if(nrank==0)  call hdf_create_emptyfile(file_path)
        call hdf_add_3Dfield(file_path, previousRHS1_x(:,:,:), "RHS1", nx_global, ny_global, nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))

        file_path=trim(recovery_RHS_dir)//"/RHS2"
        if(nrank==0)  call hdf_create_emptyfile(file_path)
        call hdf_add_3Dfield(file_path, previousRHS2_x(:,:,:), "RHS2", nx_global, ny_global, nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))


        file_path=trim(recovery_RHS_dir)//"/RHS3"
        if(nrank==0)  call hdf_create_emptyfile(file_path)
        call hdf_add_3Dfield(file_path, previousRHS3_x(:,:,:), "RHS3", nx_global, ny_global, nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))


    end subroutine save_state

    subroutine finalize()
        implicit none

        call deallocate_data

        contains

            subroutine deallocate_data()

            implicit none

            !!! VELOCITY_SOLVER DATA !!!
            deallocate(RHS1_x)
            deallocate(RHS2_x)
            deallocate(RHS3_x)
            deallocate(gradP1_x)
            deallocate(gradP2_x)
            deallocate(gradP3_x)

            deallocate(NSTERMS1)
            deallocate(NSTERMS2)
            deallocate(NSTERMS3)
            deallocate(Fext1)
            deallocate(Fext2)
            deallocate(Fext3)

            deallocate(RHS1_y)
            deallocate(RHS2_y)
            deallocate(RHS3_y)
            deallocate(gradP2_y)
            deallocate(gradP3_y)

            deallocate(RHS1_z)
            deallocate(RHS2_z)
            deallocate(RHS3_z)
            deallocate(gradP3_z)

            deallocate(previousRHS1_x)
            deallocate(previousRHS2_x)
            deallocate(previousRHS3_x)

            deallocate(q1_save)
            deallocate(q2_save)
            deallocate(q3_save)

            deallocate(q1_save_x)
            deallocate(q2_save_x)
            deallocate(q3_save_x)

            deallocate(q1_save_y)
            deallocate(q2_save_y)
            deallocate(q3_save_y)

            deallocate(q1_save_z)
            deallocate(q2_save_z)
            deallocate(q3_save_z)

            deallocate(p13_z)
            deallocate(p23_z)

            deallocate(p23_y)
            deallocate(p13_y)
            deallocate(p12_y)

            deallocate(p13_x)
            deallocate(p12_x)

            deallocate(outflow_conv1o)
            deallocate(outflow_conv2o)
            deallocate(outflow_conv3o)

            deallocate(outflow_conv1c)
            deallocate(outflow_conv2c)
            deallocate(outflow_conv3c)

            deallocate(outflow_diff1o)
            deallocate(outflow_diff2o)
            deallocate(outflow_diff3o)

            deallocate(outflow_diff1c)
            deallocate(outflow_diff2c)
            deallocate(outflow_diff3c)

        end subroutine deallocate_data

    end subroutine finalize

end module VELOCITY_solver



