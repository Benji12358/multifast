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

    implicit none

    ! Only update_velocity is callable from the outside
    public :: update_velocity, init, save_state, update_multiphysics_velocity, add_action_x, add_action_z, substract_action
    public :: Fext1, Fext2, Fext3
    public :: finalize

    private

    procedure(LAMB_solve_Poisson), pointer  :: solve_Poisson

    logical     :: previousRHS_are_available=.false.

    ! Navier-stokes equation coefficients
    real*8      :: diff21_coef, diff22_coef
    real*8      :: diff1_coef, diff2_coef , diff3_coef
    real*8      :: conv1_coef, conv2_coef, conv3_coef

    real*8, dimension(:,:,:), allocatable       :: p13_z, p13_y, p13_x
    real*8, dimension(:,:,:), allocatable       :: p12_y, p12_x
    real*8, dimension(:,:,:), allocatable       :: p23_z, p23_y

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
    real*8, dimension(:,:,:), allocatable     :: q1_save2, q2_save2, q3_save2

    real*8, dimension(:,:), allocatable     :: outflow_conv1o, outflow_conv2o, outflow_conv3o
    real*8, dimension(:,:), allocatable     :: outflow_conv1c, outflow_conv2c, outflow_conv3c
    real*8, dimension(:,:), allocatable     :: outflow_diff1o, outflow_diff2o, outflow_diff3o
    real*8, dimension(:,:), allocatable     :: outflow_diff1c, outflow_diff2c, outflow_diff3c


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

        implicit none

        ! Poisson solver initialization
        if (use_generic_poisson) then
            call poisson_init(n1,n2,n3, BC1, BC2, BC3, L1, L2, L3, istret, alpha, beta, 1)
            call generic_poisson_infos
            solve_Poisson=>LAMB_solve_Poisson
        else
            call poisson_020_init
            solve_Poisson=>ORL_solve_Poisson
        end if
        !call generic_poisson_infos

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

            allocate(q1_save2(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            q1_save2=0.d0
            allocate(q2_save2(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            q2_save2=0.d0
            allocate(q3_save2(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
            q3_save2=0.d0

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
        integer             :: j
        real*8              :: mean_grad_P1, mean_grad_P3
        integer :: mpi_err

        debugproc=-1

        !***********************************************************
        !******************  CORE OF THE DNS ***********************
        !***********************************************************

        !   Time integration implicit viscous
        !   if nb_substep=1 Adams-Bashfort if nb_substep=3 Runge-Kutta
        call set_time_coeffs
        call set_NS_coeffs


        ! ****************  RHS COMPUTATION ************************
        call perform_velocity_products
        call perform_RHS

        ! The mean pressure gradient is performed so as to fullfill the flow rate conservation
        if (streamwise==3)  then
            call perform_mean_gradP3(mean_grad_P3)
            mean_grad_P1=0.d0
        elseif (streamwise==1)  then
            call perform_mean_gradP1(mean_grad_P1)
            mean_grad_P3=0.d0
        endif

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
        ! Reinject the flow at the outlet at the inlet for fringe
        if (use_fringe) then
            !if (streamwise==3) then
            !    call velocity_transpose_x_to_z(q1_x, q2_x, q3_x)
            !endif
            call FRINGE_set_inflow
            !if (streamwise==3) then
            !    call velocity_transpose_z_to_x(q1_z, q2_z, q3_z)
            !endif
        endif

        ! ******************  UPDATE VELOCITY **********************
        ! RHS=(Non linear terms + Diffusive terms)
        ! PreviousRHS = RHS(n-1)
        ! Obtention of û*, v*, w*

        ! Store velocity at the current time step for support of an external force.
        ! The arrays "q" will be erased by the current estimation of velocity at n+1 (in the coupling process)
        ! So its value is saved in the q*_save arrays
        q1_save2=q1_x
        q2_save2=q2_x
        q3_save2=q3_x

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

            implicit none

            real*8, intent(out) :: mean_grad_P1

            real*8,dimension (n1)           :: don1
            real*8,dimension (n2)           :: don2
            real*8,dimension (n3)           :: don3
            real*8,dimension (xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))          :: tmp1
            real*8,dimension (ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))          :: tmp2
            real*8,dimension (zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))          :: tmp3

            real*8 s1tot, s1tot_glob
            integer k,i,j, mpi_err, n1e, n2e, n3e

            s1tot=0.d0
            tmp1=0.d0
            tmp2=0.d0
            tmp3=0.d0

            ! _______________________________________________________________________________________________________________________
            ! ***********************************************************************************************************************
            ! ****************************** Derivative along the Z direction *******************************************************
            ! ***********************************************************************************************************************

            !do j=1,n2m
            if (NS_Q1_BC3/=periodic) then

                n1e=(min(n1m, zend(1))-zstart(1))+1
                n2e=(min(n2m, zend(2))-zstart(2))+1

                call D2c_3Dz(q1_z, tmp3, zsize(1),n1e,zsize(2),n2e,n3, dx3, .true., NS_Q1_BC3)

                do j = zstart(2), min(n2m, zend(2))
                    do i=zstart(1), min(n1m, zend(1))   !do i=1,n1m

                        do k=1,n3m
                            s1tot=s1tot+don3(k)*dx1*dx3*cell_size_Y(j)
                        enddo

                    enddo
                enddo

            endif



            n1e=(min(n1m, yend(1))-ystart(1))+1
            n3e=(min(n3m, yend(3))-ystart(3))+1

            call D1c_MULT_3Dy(q1_y, tmp2, ysize(1),n1e,n2,ysize(3),n3e, dx2, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,1))
            call D2c_MULTACC_3Dy(q1_y, tmp2, ysize(1),n1e,n2,ysize(3),n3e, dx2, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,2))

            do k=ystart(3), min(n3m, yend(3))   !do k=1,n3m
                do i=ystart(1), min(n1m, yend(1))   !do i=1,n1m

                    tmp2(i,1,k)= ( q1_y(i,2,k)*a3_d + q1_y(i,1,k)*a2_d + q1_wall20(i,k)*a1_d )/dx2**2
                    tmp2(i,n2m,k)= ( q1_y(i,n2m-1,k)*a1_u + q1_y(i,n2m,k) *a2_u + q1_wall21(i,k)*a3_u)/dx2**2
                    do j = 1, n2m
                        s1tot=s1tot+tmp2(i,j,k)*dx1*dx3*cell_size_Y(j)
                    end do

                enddo
            enddo

            ! X ------------------------------------------------------------------------

            n2e=(min(n2m, xend(2))-xstart(2))+1
            n3e=(min(n3m, xend(3))-xstart(3))+1

            call D2c_3Dx(q1_x, tmp1, n1,xsize(2),n2e,xsize(3),n3e, dx1, .false., NS_Q1_BC1)

            do k=xstart(3), min(n3m, xend(3))   !do k=1,n3m
                do j=xstart(2), min(n2m, xend(2))   !do j=1,n2m

                    do i=1,n1m
                        s1tot=s1tot+tmp1(i,j,k)*dx1*dx3*cell_size_Y(j)
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

            mean_grad_P1=s1tot_glob/(2.d0*L1*L3*ren)


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

            ! _______________________________________________________________________________________________________________________
            ! ***********************************************************************************************************************
            ! ****************************** Derivative along the Z direction *******************************************************
            ! ***********************************************************************************************************************

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
            if (BC1==OPEN) then
                n1s=2
                n1e=n1-2
            endif
            n2s=xstart(2)
            n3s=xstart(3)
            do k = n3s, n3e
                do j = n2s, n2e
                    do i = n1s, n1e
!                     NSTERMS1(i,j,k)=(RHS1_x(i,j,k)*gam + previousRHS1_x(i,j,k)*rom - al*(gradP1_x(i,j,k) + mean_grad_P1) + al*(It_Factor_1*PreviousFext1(i,j,k) + It_Factor_2*Fext1(i,j,k)) )*dt
                      NSTERMS1(i,j,k)=(RHS1_x(i,j,k)*gam + previousRHS1_x(i,j,k)*rom - al*(gradP1_x(i,j,k) + mean_grad_P1) )*dt
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
                        NSTERMS3(i,j,k)=(RHS3_x(i,j,k)*gam + previousRHS3_x(i,j,k)*rom - al*(gradP3_x(i,j,k) + mean_grad_P3) )*dt
!                       NSTERMS3(i,j,k)= (RHS3_x(i,j,k)*gam + previousRHS3_x(i,j,k)*rom + gradP3_x(i,j,k) + Fext3(i,j,k)*al)*dt

                    end do
                end do
            end do

        end subroutine perform_NSterms

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

        implicit none

        integer, intent(in)                                     :: ns, ntime

        real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: dpdy_old1
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2)) :: dpdy_old3
        real*8                                                  :: errorsum, errorsum_glob, divy_mean
        integer                                                 :: mpi_err, s
        real*8, dimension(3)                                    :: fringe_error_velocity, fringe_error_velocity_glob
        integer                                                 :: i,j,k


        ! ******************  CORRECT VELOCITY *********************
        !  RHS3 becomes divergence(q)/(al*dt) (= laplacien (phi))
        if (BC1==OPEN) call force_flowrate_at_outflow

        ! BC and IBM are applied BEFORE correction step: => The velocity at BC and body must be adapted to be
        ! correct AFTER correction step.
        ! Use of an iterative process

        do s = 1, 1         ! Loop for iterative process

            if (streamwise==1) dpdy_old1=dphidx2_x(n1-1,:,:)
            if (streamwise==3) then
                call transpose_x_to_y(dphidx2_x, dphidx2_y)
                call transpose_y_to_z(dphidx2_y, dphidx2_z)
                dpdy_old3=dphidx2_z(:,:,n3-1)
            endif

            ! applying current estimation s of grad(p) to BC and body
            ! Vc=Vc+grad(p)[s] (Vc = desired value of velocity)
            ! after correction the velocity in body and BC will be : Vc+grad(p)[N-1]-grad(p)[N]
            ! where n is the number of iteration and will be near to u
            call pressure_gradient_transpose_x_to_z(dphidx1_x, dphidx2_x, dphidx3_x)

            if (BC1==OPEN) call apply_gradp_inoutflow_velocity_1
            if (BC3==OPEN) call apply_gradp_inoutflow_velocity_3

            if (BC1==FRINGE) call apply_gradp_fringe_velocity_1
            !if (BC3==OPEN) call apply_gradp_inoutflow_velocity_3

            if ((streamwise==1).and.(IBM_activated)) call apply_gradp_IBM_velocity_1
            if ((streamwise==3).and.(IBM_activated)) then
                call velocity_transpose_x_to_z(q1_x, q2_x, q3_x)
                ! pressure correction
                call apply_gradp_IBM_velocity_3
                call velocity_transpose_z_to_x(q1_z, q2_z, q3_z)
            endif                

            call transpose_x_to_y(q2_x, q2_y)
            call transpose_x_to_y(q3_x, q3_y)
            call transpose_y_to_z(q3_y, q3_z)

            ! Applying the boundary condition
            call apply_BC3
            call apply_BC2
            call apply_BC1

            call perform_divergence(divu_z, q3_z, q2_y, q1_x)

            call solve_Poisson(divu_z/(al*dt), dp_z)

            ! Correct velocity and update pressure
            call perform_gradp(dp_z)

            ! Convergence control of iterative process :
            if (streamwise==1) errorsum=sum(abs(dpdy_old1-dphidx2_x(n1-1,:,:)))
            if (streamwise==3) then
                call transpose_x_to_y(dphidx2_x, dphidx2_y)
                call transpose_y_to_z(dphidx2_y, dphidx2_z)
                errorsum=sum(abs(dpdy_old3-dphidx2_z(:,:,n3-1)))
            endif

            call MPI_ALLREDUCE(errorsum,errorsum_glob,1,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
            if (nrank==0) write(*,*) 'ERROR GRADP:', s, errorsum_glob

        end do

        ! Correction of velocity so as ensuring the free-divergence condition
        call correct_velocity

        ! ******************  PERFORM RESIDUAL DIVERGENCE **********
        ! RHS3_z=P*aldt
        call perform_divergence(divu_z, q3_z, q2_y, q1_x, divy_mean)
        if (nrank==0) write(*,*) 'divy_mean', divy_mean
        call transpose_z_to_y(divu_z, divu_y)
        call transpose_y_to_x(divu_y, divu_x)


        call spread_to_all_pencil(q3_z, q2_y, q1_x, dp_z)

        if (outflow_buff>0) then
            if(.not. inout_newversion) call OPEN_OLD_add_outflow(ntime, ns)
            if(inout_newversion) call OPEN_add_outflow
        endif

        if (use_fringe) then
            do j=xstart(2),min(xend(2),n2m)
                do k=xstart(3),min(xend(3),n3m)
                    fringe_error_velocity(1) = fringe_error_velocity(1) + (abs(q1_inflow(j,k))-abs(q1_x(n1-1,j,k))) / ((n2-1)*(n3-1))
                    fringe_error_velocity(2) = fringe_error_velocity(1) + (abs(q2_inflow(j,k))-abs(q2_x(n1-1,j,k))) / ((n2-1)*(n3-1))
                    fringe_error_velocity(3) = fringe_error_velocity(1) + (abs(q3_inflow(j,k))-abs(q3_x(n1-1,j,k))) / ((n2-1)*(n3-1))
                enddo
            enddo

            call MPI_ALLREDUCE(fringe_error_velocity,fringe_error_velocity_glob,1,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
            if (nrank==0) write(*,*) 'ERROR VELOCITY IN/OUT:', fringe_error_velocity_glob
        endif

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

                        q1_x(n1,j,k)=q1_save2(n1,j,k)+dphidx1_x(n1,j,k)*al*dt
                        q1_x(1,j,k)=q1_save2(1,j,k)+dphidx1_x(1,j,k)*al*dt

                        q2_x(n1,j,k)=q2_save2(n1,j,k)+dphidx2_x(n1,j,k)*al*dt
                        q2_x(1,j,k)=q2_save2(1,j,k)+dphidx2_x(1,j,k)*al*dt

                        q3_x(n1,j,k)=q3_save2(n1,j,k)+dphidx3_x(n1,j,k)*al*dt
                        q3_x(1,j,k)=q3_save2(1,j,k)+dphidx3_x(1,j,k)*al*dt

                    enddo
                enddo

            end subroutine apply_gradp_inoutflow_velocity_1

            subroutine apply_gradp_inoutflow_velocity_3
                implicit none
                integer     :: i, j, k

                do k=xstart(3),min(xend(3),n3-1)
                    do j=xstart(2),min(xend(2),n2-1)

                        q1_z(i,j,n3-1)=q1_save2(i,j,n3-1)+dphidx1_x(i,j,n3-1)*al*dt
                        q1_z(i,j,1)=q1_save2(i,j,1)+dphidx1_x(i,j,1)*al*dt

                        q2_z(i,j,n3-1)=q2_save2(i,j,n3-1)+dphidx2_x(i,j,n3-1)*al*dt
                        q2_z(i,j,1)=q2_save2(i,j,1)+dphidx2_x(i,j,1)*al*dt

                        q3_z(i,j,n3-1)=q3_save2(i,j,n3-1)+dphidx3_x(i,j,n3-1)*al*dt
                        q3_z(i,j,1)=q3_save2(i,j,1)+dphidx3_x(i,j,1)*al*dt

                    enddo
                enddo

            end subroutine apply_gradp_inoutflow_velocity_3

            subroutine apply_gradp_fringe_velocity_1
                use FRINGE_data

                implicit none
                integer     :: i, j, k

                do k=xstart(3),min(xend(3),n3-1)
                    do j=xstart(2),min(xend(2),n2-1)

                        ! take the velocity after the inflow is set
                        q1_x(1,j,k)=q1_save2(1,j,k)+dphidx1_x(1,j,k)*al*dt
                        q2_x(1,j,k)=q2_save2(1,j,k)+dphidx2_x(1,j,k)*al*dt
                        q3_x(1,j,k)=q3_save2(1,j,k)+dphidx3_x(1,j,k)*al*dt

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

                            q1_x(i,j,k)=q1_save(i,j,k)+IBM_mask1_x(i,j,k)*dphidx1_x(i,j,k)*al*dt
                            q2_x(i,j,k)=q2_save(i,j,k)+IBM_mask2_x(i,j,k)*dphidx2_x(i,j,k)*al*dt
                            q3_x(i,j,k)=q3_save(i,j,k)+IBM_mask3_x(i,j,k)*dphidx3_x(i,j,k)*al*dt

                        end do
                    enddo
                enddo

            end subroutine apply_gradp_IBM_velocity_1

            subroutine apply_gradp_IBM_velocity_3
                use IBM_settings
                use IBM_data
                implicit none
                integer     :: i, j, k

                do i = zstart(1),min(zend(1),n1-1)
                    do j = zstart(2),min(zend(2),n2-1)
                        do k = 1, n3-1

                            q1_z(i,j,k)=q1_save(i,j,k)+IBM_mask1_z(i,j,k)*dphidx1_z(i,j,k)*al*dt
                            q2_z(i,j,k)=q2_save(i,j,k)+IBM_mask2_z(i,j,k)*dphidx2_z(i,j,k)*al*dt
                            q3_z(i,j,k)=q3_save(i,j,k)+IBM_mask3_z(i,j,k)*dphidx3_z(i,j,k)*al*dt

                        end do
                    enddo
                enddo

            end subroutine apply_gradp_IBM_velocity_3


            subroutine perform_gradp(dph_z)

                use numerical_methods_settings
                use schemes_interface
                implicit none

                real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))      :: dph_z
                real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))      :: dph_x
                real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))      :: dph_y
                integer ::k,j,i, n1e, n2e, n3e

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

                if (use_generic_poisson) then

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
            end subroutine


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
            end subroutine

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
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: diff13, diff23, diff33
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: diff12, diff22, diff32
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: diff13y, diff23y, diff33y
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: diff12x, diff22x, diff32x, diff13x, diff23x, diff33x

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
            call perform_diff
        endif

        if (outflow_type==6) then
            call perform_local_cx
            call perform_diff
        endif

        if (outflow_type==6) then
            call perform_inlet_cx
            call perform_diff
        endif

        call perform_conv

        n2s=xstart(2)
        n3s=xstart(3)
        n2e=xend(2)
        n3e=xend(3)

        do k=n3s,n3e
            do j=n2s,n2e
                q1_x(n1,j,k) = q1_x(n1,j,k) + conv1(j,k) + diff1(j,k)*dt
                q1_wall11(j,k) = q1_x(n1,j,k)
                q2_x(n1-1,j,k) = q2_x(n1-1,j,k) + conv2(j,k) + diff2(j,k)*dt
                q3_x(n1-1,j,k) = q3_x(n1-1,j,k) + conv3(j,k) + diff3(j,k)*dt
            enddo
        enddo

        if (nrank==debugproc) then
            open(15, file="debugoutflow", position="append")
            write(15,*)"OUTFLOW A", sum(abs(conv1(:,:)))/(xsize(2)*xsize(3)), sum(abs(conv2(:,:)))/(xsize(2)*xsize(3)), sum(abs(conv3(:,:)))/(xsize(2)*xsize(3))
            write(15,*)"OUTFLOW B", sum(abs(diff1(:,:)))/(xsize(2)*xsize(3)), sum(abs(diff2(:,:)))/(xsize(2)*xsize(3)), sum(abs(diff3(:,:)))/(xsize(2)*xsize(3))
            close(15)
        endif

    contains

        subroutine perform_diff()
            use boundary_scheme
            use snapshot_writer
            implicit none
            integer :: i,j,k
            real*8      :: diff1_coef, diff2_coef, diff3_coef, diff22_coef, diff21_coef



            diff1_coef=dx1*sqrt(ren)
            diff2_coef=dx2*sqrt(ren)
            diff3_coef=dx3*sqrt(ren)
            diff22_coef=dx2*sqrt(ren)
            diff21_coef=dx2*ren

            ! Initializations necessary to avoid VALGRIND errors
            diff13=0.d0
            diff23=0.d0
            diff33=0.d0
            diff12=0.d0
            diff22=0.d0
            diff32=0.d0

            if (.not. allocated(outflow_diff1o)) then
            endif

            call D2c_3Dz(q1_z, diff13, zsize(1),zsize(1),zsize(2),zsize(2),n3, diff3_coef, .true., NS_Q1_BC3)
            call D2c_3Dz(q2_z, diff23, zsize(1),zsize(1),zsize(2),zsize(2),n3, diff3_coef, .true., NS_Q2_BC3)
            call D2c_3Dz(q3_z, diff33, zsize(1),zsize(1),zsize(2),zsize(2),n3, diff3_coef, .true., NS_Q3_BC3)

            do k = ystart(3), yend(3)
                do i = ystart(1), yend(1)
                    diff12(i,1,k)=( q1_y(i,2,k)*a3_d + q1_y(i,1,k)*a2_d + q1_wall20(i,k)*a1_d)            /diff2_coef**2
                    diff12(i,n2m,k)=( q1_y(i,n2m-1,k)*a1_u + q1_y(i,n2m,k)*a2_u + q1_wall21(i,k)*a3_u ) /diff2_coef**2

                    diff32(i,1,k)=( q3_y(i,2,k)*a3_d + q3_y(i,1,k)*a2_d + q3_wall20(i,k)*a1_d )           /diff2_coef**2
                    diff32(i,n2m,k)=( q3_y(i,n2m-1,k)*a1_u + q3_y(i,n2m,k) *a2_u + q3_wall21(i,k)*a3_u) /diff2_coef**2
                end do
            end do

            call D1c_MULT_3Dy(q1_y, diff12, ysize(1),ysize(1),n2,ysize(3),ysize(3), diff21_coef, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,1))
            call D2c_MULTACC_3Dy(q1_y, diff12, ysize(1),ysize(1),n2,ysize(3),ysize(3), diff22_coef, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,2))

            call D1c_MULT_3Dy(q2_y, diff22, ysize(1),ysize(1),n2,ysize(3),ysize(3), diff21_coef, .false., NS_Q2_BC2, Y_to_YTr_for_D2(:,1))
            call D2c_MULTACC_3Dy(q2_y, diff22, ysize(1),ysize(1),n2,ysize(3),ysize(3), diff22_coef, .false., NS_Q2_BC2, Y_to_YTr_for_D2(:,2))

            call D1c_MULT_3Dy(q3_y, diff32, ysize(1),ysize(1),n2,ysize(3),ysize(3), diff21_coef, .true., NS_Q3_BC2, Yc_to_YcTr_for_D2(:,1))
            call D2c_MULTACC_3Dy(q3_y, diff32, ysize(1),ysize(1),n2,ysize(3),ysize(3), diff22_coef, .true., NS_Q3_BC2, Yc_to_YcTr_for_D2(:,2))

            call transpose_z_to_y(diff13, diff13y)
            call transpose_z_to_y(diff23, diff23y)
            call transpose_z_to_y(diff33, diff33y)


            call transpose_y_to_x(diff13y, diff13x)
            call transpose_y_to_x(diff23y, diff23x)
            call transpose_y_to_x(diff33y, diff33x)

            call transpose_y_to_x(diff12, diff12x)
            call transpose_y_to_x(diff22, diff22x)
            call transpose_y_to_x(diff32, diff32x)

            outflow_diff1c = diff13x(n1,:,:) + diff12x(n1,:,:)
            outflow_diff2c = diff23x(n1-1,:,:) + diff22x(n1-1,:,:)
            outflow_diff3c = diff33x(n1-1,:,:) + diff32x(n1-1,:,:)

            diff1 = gam*outflow_diff1c + rom*outflow_diff1o
            diff2 = gam*outflow_diff2c + rom*outflow_diff2o
            diff3 = gam*outflow_diff3c + rom*outflow_diff3o

            outflow_diff1o = outflow_diff1c
            outflow_diff2o = outflow_diff2c
            outflow_diff3o = outflow_diff3c


        end subroutine perform_diff

        ! Update condition at outflow from the current time scheme, according an information propagation with phase velocity cx.
        subroutine perform_conv()
            use snapshot_writer
            implicit none

            do k=xstart(3),min(xend(3), n3-1)
                do j=xstart(2),min(xend(2), n2-1)
                    outflow_conv1c(j,k) = -cx1(j,k)*( q1_x(n1,j,k) - q1_x(n1-1,j,k) )
                    outflow_conv2c(j,k) = -cx2(j,k)*( q2_x(n1-1,j,k) - q2_x(n1-2,j,k) )
                    outflow_conv3c(j,k) = -cx3(j,k)*( q3_x(n1-1,j,k) - q3_x(n1-2,j,k) )
                enddo
            enddo

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

        implicit none

        !real*8,dimension(n3)           :: uu
        real*8,dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))           :: q2q2_y
        real*8,dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))           :: q3q3_z
        real*8,dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))           :: q1q1_x

        integer k,i,j, mpi_err, n1e, n2e, n3e

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

        if (BC3==NOSLIP) call diff_at_wall_3   ! ATTENTION

        call D2c_3Dz(q1_z, RHS1_z, zsize(1),n1e,zsize(2),n2e,n3, diff3_coef, .true., NS_Q1_BC3)
        call D1s_ACC_3Dz(p13_z, RHS1_z, zsize(1),n1e,zsize(2),n2e,n3, conv3_coef, .false., NS_P13_BC3)
        ! *********************************** CORR_JO ******************************* !
        ! La D2c_3Dz sur q2 est shifted --> argument .true. (ancienne valeur : false) !
        ! *************************************************************************** !
        call D2c_3Dz(q2_z, RHS2_z, zsize(1),n1e,zsize(2),n2e,n3, diff3_coef, .true., NS_Q2_BC3)
        call D2c_3Dz(q3_z, RHS3_z, zsize(1),n1e,zsize(2),n2e,n3, diff3_coef, .false., NS_Q3_BC3)
        call D1s_ACC_3Dz(p23_z, RHS2_z, zsize(1),n1e,zsize(2),n2e,n3, conv3_coef, .false., NS_P23_BC3)

        call D0s_3Dz(q3_z, q3q3_z, zsize(1),n1e,zsize(2),n2e,n3, NS_Q3_BC3)
        q3q3_z=q3q3_z**2
        call D1ssh_ACC_3Dz(q3q3_z, RHS3_z, zsize(1),n1e,zsize(2),n2e,n3, conv3_coef, .true., NS_P33_BC3)

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

        ! _______________________________________________________________________________________________________________________
        ! ***********************************************************************************************************************
        ! ****************************** Derivative along the Y direction *******************************************************
        ! ***********************************************************************************************************************

        call transpose_z_to_y(p23_z, p23_y)

        call transpose_z_to_y(RHS1_z, RHS1_y)
        call transpose_z_to_y(RHS2_z, RHS2_y)
        call transpose_z_to_y(RHS3_z, RHS3_y)

        n1e=ysize(1)!(min(n1m, yend(1))-ystart(1))+1
        n3e=ysize(3)!(min(n3m, yend(3))-ystart(3))+1

        if(BC2==NOSLIP) call diff_at_wall_2 ! implicit: NS_Q3_BC2=Dirichlet

        call D1c_MULTACC_3Dy(q1_y, RHS1_y, ysize(1),n1e,n2,ysize(3),n3e, diff21_coef, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,1))
        call D2c_MULTACC_3Dy(q1_y, RHS1_y, ysize(1),n1e,n2,ysize(3),n3e, diff22_coef, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,2))

        call D1c_MULTACC_3Dy(q2_y, RHS2_y, ysize(1),n1e,n2,ysize(3),n3e, diff21_coef, .false., NS_Q2_BC2, Y_to_YTr_for_D2(:,1))
        call D2c_MULTACC_3Dy(q2_y, RHS2_y, ysize(1),n1e,n2,ysize(3),n3e, diff22_coef, .false., NS_Q2_BC2, Y_to_YTr_for_D2(:,2))

        call D1c_MULTACC_3Dy(q3_y, RHS3_y, ysize(1),n1e,n2,ysize(3),n3e, diff21_coef, .true., NS_Q3_BC2, Yc_to_YcTr_for_D2(:,1))
        call D2c_MULTACC_3Dy(q3_y, RHS3_y, ysize(1),n1e,n2,ysize(3),n3e, diff22_coef, .true., NS_Q3_BC2, Yc_to_YcTr_for_D2(:,2))

        call D1s_MULTACC_3Dy(p12_y, RHS1_y, ysize(1),n1e,n2,ysize(3),n3e, conv2_coef, .false., NS_P12_BC2, Yc_to_YcTr_for_D1)
        call D0s_3Dy(q2_y, q2q2_y, ysize(1),n1e,n2,ysize(3),n3e, NS_Q2_BC2)
        q2q2_y=q2q2_y**2
        call D1ssh_MULTACC_3Dy(q2q2_y, RHS2_y, ysize(1),n1e,n2,ysize(3),n3e, conv2_coef, .true., NS_P22_BC2, Y_to_YTr_for_D1)

        call D1s_MULTACC_3Dy(p23_y, RHS3_y, ysize(1),n1e,n2,ysize(3),n3e, conv2_coef, .false., NS_P23_BC2, Yc_to_YcTr_for_D1)
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

        n2e=xsize(2)!(min(n2m, xend(2))-xstart(2))+1
        n3e=xsize(3)!(min(n3m, xend(3))-xstart(3))+1

        if((BC1==NOSLIP).or.(BC1==PSEUDO_PERIODIC).or.(BC1==OPEN)) call diff_at_wall_1   ! ATTENTION

        call D0s_3Dx(q1_x, q1q1_x, n1,xsize(2),n2e,xsize(3),n3e, NS_Q1_BC1)
        q1q1_x=q1q1_x**2
        call D1ssh_ACC_3Dx(q1q1_x, RHS1_x, n1,xsize(2),n2e,xsize(3),n3e, conv1_coef, .true., NS_P11_BC1)
        call D1s_ACC_3Dx(p12_x, RHS2_x, n1,xsize(2),n2e,xsize(3),n3e, conv1_coef, .false., NS_P12_BC1)
        call D1s_ACC_3Dx(p13_x, RHS3_x, n1,xsize(2),n2e,xsize(3),n3e, conv1_coef, .false., NS_P13_BC1)

        call D2c_ACC_3Dx(q1_x, RHS1_x, n1,xsize(2),n2e,xsize(3),n3e, diff1_coef, .false., NS_Q1_BC1)
        call D2c_ACC_3Dx(q2_x, RHS2_x, n1,xsize(2),n2e,xsize(3),n3e, diff1_coef, .true., NS_Q2_BC1)
        call D2c_ACC_3Dx(q3_x, RHS3_x, n1,xsize(2),n2e,xsize(3),n3e, diff1_coef, .true., NS_Q3_BC1)

        return

    contains

        subroutine diff_at_wall_1()

            ! d2f=a1*f0 + a2*f1 + a3*f2
            ! Regular case: down    a1=2/3, a2=-1 a3=1/3
            ! Regular case: up      a1=1/3, a2=-1 a3=2/3
            ! ATTENTION
            do k=xstart(3), min(n3m, xend(3))       !do k=1,n3m
                do j=xstart(2), min(n2m, xend(2))       !do j=1,n2m
                    RHS2_x(1,j,k)   =RHS2_x(1,j,k)  +( q2_x(2,j,k)       /3.d0 - q2_x(1,j,k)      + q2_wall10(j,k)*2.d0/3.d0)  /diff1_coef**2
                    RHS2_x(n1m,j,k) =RHS2_x(n1m,j,k)+( q2_x(n1m-1,j,k)   /3.d0 - q2_x(n1m,j,k)    + q2_wall11(j,k)*2.d0/3.d0 ) /diff1_coef**2

                    RHS3_x(1,j,k)   =RHS3_x(1,j,k)  +( q3_x(2,j,k)       /3.d0 - q3_x(1,j,k)      + q3_wall10(j,k)*2.d0/3.d0 ) /diff1_coef**2
                    RHS3_x(n1m,j,k) =RHS3_x(n1m,j,k)+( q3_x(n1m-1,j,k)   /3.d0 - q3_x(n1m,j,k)    + q3_wall11(j,k)*2.d0/3.d0)  /diff1_coef**2
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

                    RHS1_y(i,1,k)=RHS1_y(i,1,k)+( q1_y(i,2,k)*a3_d + q1_y(i,1,k)*a2_d + q1_wall20(i,k)*a1_d)            /diff2_coef**2
                    RHS1_y(i,n2m,k)=RHS1_y(i,n2m,k)+( q1_y(i,n2m-1,k)*a1_u + q1_y(i,n2m,k)*a2_u + q1_wall21(i,k)*a3_u ) /diff2_coef**2

                    RHS3_y(i,1,k)=RHS3_y(i,1,k)+( q3_y(i,2,k)*a3_d + q3_y(i,1,k)*a2_d + q3_wall20(i,k)*a1_d )           /diff2_coef**2
                    RHS3_y(i,n2m,k)=RHS3_y(i,n2m,k)+( q3_y(i,n2m-1,k)*a1_u + q3_y(i,n2m,k) *a2_u + q3_wall21(i,k)*a3_u) /diff2_coef**2
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
                    RHS1_z(i,j,1)   =RHS1_z(i,j,1)  +( q1_z(i,j,2)       /3.d0 - q1_z(i,j,1)      + q1_wall30(i,j)*2.d0/3.d0 ) /diff3_coef**2
                    RHS1_z(i,j,n3m) =RHS1_z(i,j,n3m)+( q1_z(i,j,n3m-1)   /3.d0 - q1_z(i,j,n3m)    + q1_wall31(i,j)*2.d0/3.d0)  /diff3_coef**2

                    RHS2_z(i,j,1)   =RHS2_z(i,j,1)  +( q2_z(i,j,2)       /3.d0 - q2_z(i,j,1)      + q2_wall30(i,j)*2.d0/3.d0)  /diff3_coef**2
                    RHS2_z(i,j,n3m) =RHS2_z(i,j,n3m)+( q2_z(i,j,n3m-1)   /3.d0 - q2_z(i,j,n3m)    + q2_wall31(i,j)*2.d0/3.d0 ) /diff3_coef**2
                enddo
            enddo

        end subroutine diff_at_wall_3

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

            call force_velocity(q1_save2, q2_save2, q3_save2, xsize, vel_term1, vel_term2, vel_term3)

            if ((BC1==UNBOUNDED).or.(BC1==FRINGE)) n1s=xstart(1)
            if ((BC1/=UNBOUNDED).and.(BC1/=FRINGE)) n1s=max(2,xstart(1))
            n2s=xstart(2)
            n3s=xstart(3)

            do k = n3s, n3e
                do j = n2s, n2e
                    do i = n1s, n1e
                        RHS1i = NSTERMS1(i,j,k) + Fext1(i, j, k)
                        q1_x(i,j,k) = q1_save2(i,j,k) + RHS1i + IBM_mask1_x(i,j,k)*(-RHS1i + vel_term1(i,j,k))
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
                        RHS2i = NSTERMS2(i,j,k) + Fext2(i, j, k)
                        q2_x(i,j,k) = q2_save2(i,j,k) + RHS2i + IBM_mask2_x(i,j,k)*(-RHS2i + vel_term2(i,j,k))
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
                        RHS3i = NSTERMS3(i,j,k) + Fext3(i, j, k)
                        q3_x(i,j,k) = q3_save2(i,j,k) + RHS3i + IBM_mask3_x(i,j,k)*(-RHS3i + vel_term3(i,j,k))
                    end do
                end do
            end do

        else

            if ((BC1==UNBOUNDED).or.(BC1==FRINGE)) n1s=xstart(1)
            if ((BC1/=UNBOUNDED).and.(BC1/=FRINGE)) n1s=max(2,xstart(1))
            if (BC1==OPEN) then
                n1s=2
                n1e=n1-2
            endif
            n2s=xstart(2)
            n3s=xstart(3)
            do k = n3s, n3e
                do j = n2s, n2e
                    do i = n1s, n1e
                        q1_x(i,j,k) = q1_save2(i,j,k) + NSTERMS1(i,j,k) + Fext1(i,j,k)
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
                        q2_x(i,j,k) = q2_save2(i,j,k) + NSTERMS2(i,j,k) + Fext2(i,j,k)
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
                        q3_x(i,j,k) = q3_save2(i,j,k) + NSTERMS3(i,j,k) + Fext3(i,j,k)
                    end do
                end do
            end do

        end if
        ! q3_x, q2_x and w_c contains u*,v* and w*

    end subroutine perform_intermediate_velocity

    ! Perform UV, UW, VW at the center of edges
    subroutine perform_velocity_products()

        implicit none

        integer k,j,i

        ! X ------------------------------------------------------------------------

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

        call D0ssh_3Dx(q3_x, p13_x, n1,xsize(2),xsize(2),xsize(3),xsize(3), NS_Q3_BC1)
        call D0ssh_3Dx(q2_x, p12_x, n1,xsize(2),xsize(2),xsize(3),xsize(3), NS_Q2_BC1)

        ! Y ------------------------------------------------------------------------
        call transpose_x_to_y(p12_x, p12_y)

        if ((BC2==NOSLIP).or.(BC2==FREESLIP)) then

            do k=ystart(3), yend(3)      !do k=1,n3
                do i=ystart(1),yend(1)       !do i=1,n1
                    p23_y(i, 1,  k)   =q3_wall20(i, k)
                    p23_y(i, n2, k)   =q3_wall21(i, k)

                    p12_y(i, 1,  k)   =p12_y(i, 1,  k) * q1_wall20(i, k)
                    p12_y(i, n2, k)   =p12_y(i, n2, k) * q1_wall21(i, k)
                enddo
            enddo

        end if

        call D0ssh_3Dy(q3_y, p23_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_Q3_BC2)
        call D0ssh_MULT_3Dy(q1_y, p12_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_Q1_BC2)


        ! Z ------------------------------------------------------------------------
        call transpose_y_to_z(p23_y, p23_z)

        call transpose_x_to_y(p13_x, p13_y)
        call transpose_y_to_z(p13_y, p13_z)


        ! ATTENTION
        if ((BC3==NOSLIP).or.(BC3==FREESLIP)) then

            do i=zstart(1), zend(1)
                do j=zstart(2), zend(2)
                    p23_z(i, j,  1)   =p23_z(i, j, 1)   * q2_wall30(i, j)
                    p23_z(i, j, n3)   =p23_z(i, j, n3)  * q2_wall31(i, j)

                    p13_z(i, j,  1)   =p13_z(i, j, 1)   * q1_wall30(i, j)
                    p13_z(i, j, n3)   =p13_z(i, j, n3)  * q1_wall31(i, j)
                end do
            enddo

        end if
        ! ENDATTENTION

        call D0ssh_MULT_3Dz(q2_z, p23_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_Q2_BC3)
        call D0ssh_MULT_3Dz(q1_z, p13_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_Q2_BC3)

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

            deallocate(q1_save2)
            deallocate(q2_save2)
            deallocate(q3_save2)

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



