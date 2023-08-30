module following_velocity_solver

    use schemes3D_interface
    use time_schemes
    use following_physical_fields
    use boundaries
    use following_mesh

    use following_irregular_derivative_coefficients
    use DNS_settings
    use following_poisson_020_solver, poisson_020_init=>init, solve_Poisson=>solve_Poisson

    use decomp_2d
    use mpi

    use following_data

    implicit none

    ! Only update_velocity is callable from the outside
    public :: update_velocity, init, add_action_x, substract_action
    public :: Fext1, Fext2, Fext3, mean_grad_P1, mean_grad_P3, mean_grad_P1_viscous
    public :: finalize

    private

    logical     :: previousRHS_are_available=.false.

    ! Mean pressure gradient, used for both cases, with IBM and without IBM
    real*8              :: mean_grad_P1, mean_grad_P3, mean_grad_P1_viscous

    ! Navier-stokes equation coefficients
    real*8      :: diff21_coef, diff22_coef
    real*8      :: diff1_coef, diff2_coef , diff3_coef
    real*8      :: conv1_coef, conv2_coef, conv3_coef

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

    real*8  :: al,gam,rom
    logical :: first_time_advancement=.true.

    integer :: NNN, NN
    integer :: debugproc

contains

    subroutine init()

        use run_ctxt_data
        use numerical_methods_settings
        use following_fringe_data
        ! use FRINGE_dao, only:fringe_infos

        implicit none

        call poisson_020_init

        ! if (use_fringe) call fringe_infos

        ! Allocate AND initialize the array used by the VELOCITY_solver (RHS, previousRHS etc.)
        call allocate_data

    contains

        subroutine allocate_data()

            implicit none

            ! RHS arrays allocations
            allocate(RHS1_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            RHS1_x=0.d0
            allocate(RHS2_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            RHS2_x=0.d0
            allocate(RHS3_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            RHS3_x=0.d0
            allocate(gradP1_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            gradP1_x=0.d0
            allocate(gradP2_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            gradP2_x=0.d0
            allocate(gradP3_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            gradP3_x=0.d0

            allocate(NSTERMS1(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            NSTERMS1=0.d0
            allocate(NSTERMS2(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            NSTERMS2=0.d0
            allocate(NSTERMS3(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            NSTERMS3=0.d0
            allocate(Fext1(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            Fext1=0.d0
            allocate(Fext2(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            Fext2=0.d0
            allocate(Fext3(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            Fext3=0.d0

            allocate(RHS1_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
            RHS1_y=0.d0
            allocate(RHS2_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
            RHS2_y=0.d0
            allocate(RHS3_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
            RHS3_y=0.d0
            allocate(gradP2_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
            gradP2_y=0.d0
            allocate(gradP3_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
            gradP3_y=0.d0

            allocate(RHS1_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
            RHS1_z=0.d0
            allocate(RHS2_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
            RHS2_z=0.d0
            allocate(RHS3_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
            RHS3_z=0.d0
            allocate(gradP3_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
            gradP3_z=0.d0


            allocate(previousRHS1_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            previousRHS1_x=0.d0
            allocate(previousRHS2_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            previousRHS2_x=0.d0
            allocate(previousRHS3_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            previousRHS3_x=0.d0

            allocate(q1_save_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            q1_save_x=0.d0
            allocate(q2_save_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            q2_save_x=0.d0
            allocate(q3_save_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            q3_save_x=0.d0

            allocate(q1_save(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            q1_save=0.d0
            allocate(q2_save(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            q2_save=0.d0
            allocate(q3_save(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            q3_save=0.d0

            ! Cross quantities array allocation
            allocate(p13_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
            p13_z=0.d0
            allocate(p23_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
            p23_z=0.d0
            allocate(p12_z(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3)))
            p12_z=0.d0

            allocate(p23_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
            p23_y=0.d0
            allocate(p13_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
            p13_y=0.d0
            allocate(p12_y(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3)))
            p12_y=0.d0

            allocate(p13_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            p13_x=0.d0
            allocate(p12_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            p12_x=0.d0
            allocate(p23_x(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)))
            p23_x=0.d0


        end subroutine allocate_data

    end subroutine init


    ! This subroutine transpose all physical quantities. Theses quantities must been known for all sub-domains
    ! orientations.
    ! At entry:
    !   - U, Pr     must be in      Z-configuration
    !   - V         must be in      Y-configuration
    !   - W         must be in      X-configuration


    subroutine update_velocity(ntime, ns)

        use numerical_methods_settings
        use following_velocity_bc_controller
        use following_velocity_operations

        use following_fringe_data
        use following_fringe_solver, only:    &
        FRINGE_set_inflow=>set_inflow

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

        !   Time integration implicit viscous
        !   if nb_substep=1 Adams-Bashfort if nb_substep=3 Runge-Kutta
        call set_time_coeffs
        call set_NS_coeffs

        ! ****************  RHS COMPUTATION ************************
        call perform_velocity_products
        call perform_RHS

        ! Define the velocity field at the inlet/outlet for an open boundary problem
        ! MUST BE AFTER the RHS calculation: the inflow and outflow obtained are for the next time step (n+1)
        ! and must not been applied to the current step velocity field
        ! MUST BE BEFORE calculation of u* : the outflow condition is performed from velocity field u at the current
        ! time step

        ! Reinject the flow at the outlet at the inlet for fringe
        ! if (use_fringe) then
        !     call FRINGE_set_inflow
        ! endif

        call perform_mean_gradP1(mean_grad_P1)
        mean_grad_P3=0.d0

        if (nrank==0) write(*,*) '### Following channel mean_grad_P1', mean_grad_P1

        RHS1_x =  RHS1_x - mean_grad_P1
        RHS3_x =  RHS3_x - mean_grad_P3

        ! ******************  UPDATE VELOCITY **********************
        ! RHS=(Non linear terms + Diffusive terms)
        ! PreviousRHS = RHS(n-1)
        ! Obtention of û*, v*, w*

        ! If MHD not used, external forces are added after computation of pressure gradient
        RHS1_x = RHS1_x + Fext1
        RHS2_x = RHS2_x + Fext2
        RHS3_x = RHS3_x + Fext3

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

        q1_save=q1_x
        q2_save=q2_x
        q3_save=q3_x

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
            use following_fringe_data, only: use_fringe
            use following_mesh

            implicit none

            real*8, intent(out) :: mean_grad_P1
            real*8,dimension (xstart_f(1):xend_f(1),xstart_f(2):xend_f(2),xstart_f(3):xend_f(3))           :: tmp1, tmp2_x, tmp3_x
            real*8,dimension (ystart_f(1):yend_f(1),ystart_f(2):yend_f(2),ystart_f(3):yend_f(3))           :: tmp2, tmp3_y
            real*8,dimension (zstart_f(1):zend_f(1),zstart_f(2):zend_f(2),zstart_f(3):zend_f(3))           :: tmp3

            real*8 s1tot, s1tot_glob, total_volume
            integer k,i,j, mpi_err, n1e, n2e, n3e, n1s, n2s, n3s

            s1tot=0.d0
            s1tot_glob=0.d0
            tmp1=0.d0
            tmp2=0.d0
            tmp3=0.d0

            ! Z ------------------------------------------------------------------------

            n1e=(min(n1m, zend_f(1))-zstart_f(1))+1
            n2e=(min(n2m, zend_f(2))-zstart_f(2))+1

            call D2c_3Dz(q1_z, tmp3, zsize_f(1),n1e,zsize_f(2),n2e,n3, dx3, .true., NS_Q1_BC3)

            ! Y ------------------------------------------------------------------------

            n1e=(min(n1m, yend_f(1))-ystart_f(1))+1
            n3e=(min(n3m, yend_f(3))-ystart_f(3))+1

            call D1c_MULT_3Dy(q1_y, tmp2, ysize_f(1),n1e,n2,ysize_f(3),n3e, dx2, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,1))
            call D2c_MULTACC_3Dy(q1_y, tmp2, ysize_f(1),n1e,n2,ysize_f(3),n3e, dx2, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,2))

            do k=ystart_f(3), min(n3m, yend_f(3))   !do k=1,n3m
                do i=ystart_f(1), min(n1m, yend_f(1))   !do i=1,n1m

                    tmp2(i,1,k)= ( q1_y(i,2,k)*a3_d + q1_y(i,1,k)*a2_d + q1_wall20(i,k)*a1_d )/dx2**2
                    tmp2(i,n2m,k)= ( q1_y(i,n2m-1,k)*a1_u + q1_y(i,n2m,k) *a2_u + q1_wall21(i,k)*a3_u)/dx2**2

                enddo
            enddo

            ! X ------------------------------------------------------------------------

            n2e=(min(n2m, xend_f(2))-xstart_f(2))+1
            n3e=(min(n3m, xend_f(3))-xstart_f(3))+1

            call D2c_3Dx(q1_x, tmp1, n1,xsize_f(2),n2e,xsize_f(3),n3e, dx1, .false., NS_Q1_BC1)

            call transpose_y_to_x(tmp2, tmp2_x, decomp_following)

            call transpose_z_to_y(tmp3, tmp3_y, decomp_following)
            call transpose_y_to_x(tmp3_y, tmp3_x, decomp_following)

            do k=xstart_f(3), min(n3m, xend_f(3))   !do k=1,n3m
                do j=xstart_f(2), min(n2m, xend_f(2))   !do j=1,n2m

                    ! only in the fringe region
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
            use IBM_settings, only: ren_tau_ibm 
            use following_fringe_data, only: use_fringe
            use DNS_settings, only:dt

            implicit none

            real*8, intent(out) :: mean_grad_P1

            real*8 :: ut1, ut1_glob,total_volume
            integer k,i,j, mpi_err, n1e, n2e, n3e, n1s, n2s, n3s

            ! total_volume = 2.d0*L1*L3
            ! if (use_fringe) total_volume = total_volume/(1+fringe_length)

            ut1=0.d0
            ut1_glob=0.d0
            do k=xstart_f(3),min(xend_f(3),n3-1)
                do j=xstart_f(2),min(xend_f(2),n2-1)
                    ut1=ut1+q1_x(1,j,k)*(Y(j+1)-Y(j))*dx3
                enddo
            enddo

            call MPI_ALLREDUCE(ut1, ut1_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
            ut1_glob=ut1_glob/(L3*L2)

            ! if (nrank.eq.0) then
            !     write(*,*) 'instantaneous flowrate', ut1_glob
            !     write(*,*) 'initial flowrate', initial_flowrate
            !     ! write(*,*) 'added dp/dx', (1/total_volume)*(abs(ut1_glob-initial_flowrate))/(dt)
            ! endif

            ! mean_grad_P1 = mean_grad_P1 - (1/total_volume)*(abs(ut1_glob-initial_flowrate))/(dt)
            mean_grad_P1 = - (ren_tau_ibm / ren )**2

            return
        end subroutine

        subroutine perform_NSterms()
            implicit none
            integer             :: i, j, k
            integer             :: n1s, n1e, n2s,n2e, n3s,n3e


            ! Transpositions of grad(P) in X configuration for the time advancement
            call transpose_z_to_y(gradP3_z, gradP3_y, decomp_following)
            call transpose_y_to_x(gradP3_y, gradP3_x, decomp_following)
            call transpose_y_to_x(gradP2_y, gradP2_x, decomp_following)

            n1e=min(n1-1, xend_f(1))
            n2e=min(n2-1, xend_f(2))
            n3e=min(n3-1, xend_f(3))

            if ((BC1==UNBOUNDED).or.(BC1==FRINGE)) n1s=xstart_f(1)
            if ((BC1/=UNBOUNDED).and.(BC1/=FRINGE)) n1s=max(2,xstart_f(1))
            n2s=xstart_f(2)
            n3s=xstart_f(3)
            do k = n3s, n3e
                do j = n2s, n2e
                    do i = n1s, n1e
!                     NSTERMS1(i,j,k)=(RHS1_x(i,j,k)*gam + previousRHS1_x(i,j,k)*rom - al*(gradP1_x(i,j,k) + mean_grad_P1) + al*(It_Factor_1*PreviousFext1(i,j,k) + It_Factor_2*Fext1(i,j,k)) )*dt
                      NSTERMS1(i,j,k)=(RHS1_x(i,j,k)*gam + previousRHS1_x(i,j,k)*rom - al*gradP1_x(i,j,k) )*dt
                    end do
                end do
            end do


            if (BC2==UNBOUNDED) n2s=xstart_f(2)
            if (BC2/=UNBOUNDED) n2s=max(2,xstart_f(2))
            ! On ne touche que q2[2..n1-2]
            if (BC1==OPEN) then
                n1s=2
                n1e=n1-2
            endif

            n3s=xstart_f(3)

            do k = n3s, n3e
                do j = n2s, n2e
                    do i = n1s, n1e
                    !WARNING : Integral of Fext2 MUST be 0
                        NSTERMS2(i,j,k)=(RHS2_x(i,j,k)*gam + previousRHS2_x(i,j,k)*rom - al*gradP2_x(i,j,k) )*dt
                         !NSTERMS2(i,j,k)= (RHS2_x(i,j,k)*gam + previousRHS2_x(i,j,k)*rom + gradP2_x(i,j,k) )*dt
                    end do
                end do
            end do

            if (BC3==UNBOUNDED) n3s=xstart_f(3)
            if (BC3/=UNBOUNDED) n3s=max(2,xstart_f(3))

            ! On ne touche que q3[2..n1-2]
            if (BC1==OPEN) then
                n1s=2
                n1e=n1-2
            endif

            n2s=xstart_f(2)

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

    subroutine final_velocity(ntime, ns)

        use following_physical_fields
        use DNS_settings
        use following_velocity_operations
        use numerical_methods_settings
        use following_velocity_bc_controller

        use following_fringe_data

        use mpi

        implicit none

        integer, intent(in)                                     :: ns, ntime

        real*8, dimension(xstart_f(2):xend_f(2),xstart_f(3):xend_f(3)) :: dpdy_old1
        real*8, dimension(xstart_f(1):xend_f(1),xstart_f(2):xend_f(2),xstart_f(3):xend_f(3)) :: dpdx_old, dpdy_old, dpdz_old
        real*8, dimension(zstart_f(1):zend_f(1),zstart_f(2):zend_f(2)) :: dpdy_old3
        real*8                                                  :: errorsum, errorsum_glob, divy_mean
        integer                                                 :: mpi_err, s
        integer                                                 :: i,j,k
        real*8                                                  :: max_q1, min_q1, max_q2, min_q2, max_q3, min_q3
        real*8                                                  :: max_glob_q1, min_glob_q1, max_glob_q2, min_glob_q2, max_glob_q3, min_glob_q3
        ! real*8                                                  :: flowrate_in, flowrate_in_glob, flowrate_fringe_in, flowrate_fringe_in_glob, flowrate_out, flowrate_out_glob

        real*8 :: t_bef, t_aft

        character(200)                :: snaps_dir, snap_dir
        character(200)    :: file_path, snap_path

        ! BC and IBM are applied BEFORE correction step: => The velocity at BC and body must be adapted to be
        ! correct AFTER correction step.
        ! Use of an iterative process

        do s = 1, 1         ! Loop for iterative process

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

            ! if (BC1==FRINGE) call apply_gradp_fringe_velocity_1

            call transpose_x_to_y(q2_x, q2_y, decomp_following)
            call transpose_x_to_y(q3_x, q3_y, decomp_following)
            call transpose_y_to_z(q3_y, q3_z, decomp_following)

            ! Applying the boundary condition
            call apply_BC3
            call apply_BC2
            call apply_BC1

            ! use the coarsest grid
            call perform_divergence(divu_z, q3_z, q2_y, q1_x)

            call solve_Poisson(divu_z/(al*dt), dp_z)

            ! Correct velocity and update pressure
            call perform_gradp(dp_z)

        end do

        call correct_velocity

        ! ******************  PERFORM RESIDUAL DIVERGENCE **********
        call spread_to_all_pencil(q3_z, q2_y, q1_x, dp_z)

        ! RHS3_z=P*aldt
        call perform_divergence(divu_z, q3_z, q2_y, q1_x, divy_mean)
        if (nrank==0) write(*,*) '### Following channel divy_mean', divy_mean
        call transpose_z_to_y(divu_z, divu_y, decomp_following)
        call transpose_y_to_x(divu_y, divu_x, decomp_following)

        contains

            subroutine apply_gradp_fringe_velocity_1
                use following_fringe_data

                implicit none
                integer     :: i, j, k

                do k=xstart_f(3),min(xend_f(3),n3-1)
                    do j=xstart_f(2),min(xend_f(2),n2-1)

                        ! take the velocity after the inflow is set
                        q1_x(1,j,k)=q1_save_x(1,j,k)+dphidx1_x(1,j,k)*al*dt
                        q2_x(1,j,k)=q2_save_x(1,j,k)+dphidx2_x(1,j,k)*al*dt
                        q3_x(1,j,k)=q3_save_x(1,j,k)+dphidx3_x(1,j,k)*al*dt

                        ! do i=n_fringe_start,n1-1

                        !     q1_x(i,j,k)=q1_save(i,j,k)+dphidx1_x(i,j,k)*al*dt
                        !     q2_x(i,j,k)=q2_save(i,j,k)+dphidx2_x(i,j,k)*al*dt
                        !     q3_x(i,j,k)=q3_save(i,j,k)+dphidx3_x(i,j,k)*al*dt

                        ! enddo

                    enddo
                enddo

            end subroutine apply_gradp_fringe_velocity_1

            subroutine perform_gradp(dph_z)

                use numerical_methods_settings
                use schemes_interface

                use numerical_methods_settings, only: schemes_configuration

                implicit none

                real*8, dimension(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3))      :: dph_z
                real*8, dimension(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3))      :: dph_x
                real*8, dimension(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3))      :: dph_y
                integer ::k,j,i,n, n1e, n2e, n3e

                ! Gradx(P) -----------------------------------------------------------------
                n1e=zsize_f(1)
                n2e=zsize_f(2)

                call D1ssh_3Dz(dph_z, dphidx3_z, zsize_f(1),n1e,zsize_f(2),n2e,n3, dx3, .true., POISSON_PR_BC3)

                call transpose_z_to_y(dphidx3_z, dph_y, decomp_following)
                call transpose_y_to_x(dph_y, dphidx3_x, decomp_following)         ! dph_y is used as temporary array to switch from z to x pencil

                ! Grady(P) -----------------------------------------------------------------
                call transpose_z_to_y(dph_z, dph_y, decomp_following)

                n1e=ysize_f(1)
                n3e=ysize_f(3)

                call D1ssh_MULT_3Dy(dph_y, dphidx2_y, ysize_f(1),n1e,n2,ysize_f(3),n3e, dx2, .true., POISSON_PR_BC2, Y_to_YTr_for_D1)

                call transpose_y_to_x(dphidx2_y, dphidx2_x, decomp_following)

                ! Gradz(P) -----------------------------------------------------------------
                call transpose_y_to_x(dph_y, dph_x, decomp_following)

                n2e=xsize_f(2)
                n3e=xsize_f(3)

                call D1ssh_3Dx(dph_x, dphidx1_x, n1,xsize_f(2),n2e,xsize_f(3),n3e, dx1, .true., POISSON_PR_BC1)

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
                do j = zstart_f(2), min(n2m, zend_f(2))
                    do i=zstart_f(1), min(n1m, zend_f(1))
                        q3_z(i,j,:)=q3_z(i,j,:)-dphidx3_z(i,j,:)*al*dt
                    enddo
                enddo


                ! V correction

                !do k=1,n3m
                do k=ystart_f(3), min(n3m, yend_f(3))
                    !do i=1,n1m
                    do i=ystart_f(1), min(n1m, yend_f(1))
                        q2_y(i,:,k)=q2_y(i,:,k)-dphidx2_y(i,:,k)*al*dt
                    enddo
                enddo

                ! X ------------------------------------------------------------------------

                do k=xstart_f(3), min(n3m, xend_f(3))

                    !do j=1,n2m
                    do j=xstart_f(2), min(n2m, xend_f(2))
                        q1_x(:,j,k)=q1_x(:,j,k)-dphidx1_x(:,j,k)*al*dt

!                        pr_x(:,j,k)=0.d0!pr_x(:,j,k) + dph_x(:,j,k)/aldt

                    enddo
                enddo

                return
            end subroutine correct_velocity

    end subroutine final_velocity


    subroutine add_action_x(fx1, fx2, fx3)
        implicit none

        real*8, dimension(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)), intent(in)  :: fx1, fx2, fx3

        Fext1 = fx1
        Fext2 = fx2
        Fext3 = fx3

    end subroutine add_action_x

    subroutine substract_action(fx1, fx2, fx3)
        implicit none

        real*8, dimension(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3)), intent(in)  :: fx1, fx2, fx3

        Fext1 = Fext1 - fx1
        Fext2 = Fext2 - fx2
        Fext3 = Fext3 - fx3

    end subroutine substract_action

    !GENERIC
    subroutine perform_RHS()

        use boundary_scheme
        use following_mesh

        implicit none

        !real*8,dimension(n3)           :: uu
        real*8,dimension(xstart_f(1):xend_f(1), xstart_f(2):xend_f(2), xstart_f(3):xend_f(3))           :: q1q1_x, tmp1_x, tmp2_x, tmp3_x
        real*8,dimension(ystart_f(1):yend_f(1), ystart_f(2):yend_f(2), ystart_f(3):yend_f(3))           :: q2q2_y, tmp1_y, tmp2_y, tmp3_y
        real*8,dimension(zstart_f(1):zend_f(1), zstart_f(2):zend_f(2), zstart_f(3):zend_f(3))           :: q3q3_z, tmp1_z, tmp2_z, tmp3_z

        integer :: k,i,j,n, mpi_err, n1e, n2e, n3e, n1s, n2s, n3s

        ! VALGRIND
        q1q1_x=0.d0
        q2q2_y=0.d0
        q3q3_z=0.d0

        ! _______________________________________________________________________________________________________________________
        ! ***********************************************************************************************************************
        ! ****************************** Derivative along the Z direction *******************************************************
        ! ***********************************************************************************************************************

        n1e= zsize_f(1)!(min(n1m, zend(1))-zstart(1))+1
        n2e=zsize_f(2)!(min(n2m, zend(2))-zstart(2))+1

        ! ******** CORR_JO ********* !
        ! REinitialisation des RHS_z !
        ! ************************** !

        RHS1_z=0.d0
        RHS2_z=0.d0
        RHS3_z=0.d0

        tmp1_z = 0.d0
        tmp2_z = 0.d0
        tmp3_z = 0.d0

        ! For RHS1
        call D2c_3Dz(q1_z, tmp1_z, zsize_f(1),n1e,zsize_f(2),n2e,n3, diff3_coef, .true., NS_Q1_BC3)
        call D1s_ACC_3Dz(p13_z, tmp1_z, zsize_f(1),n1e,zsize_f(2),n2e,n3, conv3_coef, .false., NS_P13_BC3)
        ! *********************************** CORR_JO ******************************* !
        ! La D2c_3Dz sur q2 est shifted --> argument .true. (ancienne valeur : false) !
        ! *************************************************************************** !

        ! For RHS2
        call D2c_3Dz(q2_z, tmp2_z, zsize_f(1),n1e,zsize_f(2),n2e,n3, diff3_coef, .true., NS_Q2_BC3)
        call D1s_ACC_3Dz(p23_z, tmp2_z, zsize_f(1),n1e,zsize_f(2),n2e,n3, conv3_coef, .false., NS_P23_BC3)

        ! For RHS3
        call D2c_3Dz(q3_z, tmp3_z, zsize_f(1),n1e,zsize_f(2),n2e,n3, diff3_coef, .false., NS_Q3_BC3)
        call D0s_3Dz(q3_z, q3q3_z, zsize_f(1),n1e,zsize_f(2),n2e,n3, NS_Q3_BC3)
        q3q3_z=q3q3_z**2
        call D1ssh_ACC_3Dz(q3q3_z, tmp3_z, zsize_f(1),n1e,zsize_f(2),n2e,n3, conv3_coef, .true., NS_P33_BC3)

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

        if (BC3==NOSLIP) call diff_at_wall_3   ! ATTENTION

        RHS1_z = RHS1_z + tmp1_z
        RHS2_z = RHS2_z + tmp2_z
        RHS3_z = RHS3_z + tmp3_z

        ! _______________________________________________________________________________________________________________________
        ! ***********************************************************************************************************************
        ! ****************************** Derivative along the Y direction *******************************************************
        ! ***********************************************************************************************************************

        call transpose_z_to_y(p23_z, p23_y, decomp_following)

        call transpose_z_to_y(RHS1_z, RHS1_y, decomp_following)
        call transpose_z_to_y(RHS2_z, RHS2_y, decomp_following)
        call transpose_z_to_y(RHS3_z, RHS3_y, decomp_following)

        tmp1_y = 0.d0
        tmp2_y = 0.d0
        tmp3_y = 0.d0

        n1e=ysize_f(1)!(min(n1m, yend(1))-ystart(1))+1
        n3e=ysize_f(3)!(min(n3m, yend(3))-ystart(3))+1

        ! For RHS1
        call D1c_MULTACC_3Dy(q1_y, tmp1_y, ysize_f(1),n1e,n2,ysize_f(3),n3e, diff21_coef, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,1))
        call D2c_MULTACC_3Dy(q1_y, tmp1_y, ysize_f(1),n1e,n2,ysize_f(3),n3e, diff22_coef, .true., NS_Q1_BC2, Yc_to_YcTr_for_D2(:,2))

        ! For RHS2
        call D1c_MULTACC_3Dy(q2_y, tmp2_y, ysize_f(1),n1e,n2,ysize_f(3),n3e, diff21_coef, .false., NS_Q2_BC2, Y_to_YTr_for_D2(:,1))
        call D2c_MULTACC_3Dy(q2_y, tmp2_y, ysize_f(1),n1e,n2,ysize_f(3),n3e, diff22_coef, .false., NS_Q2_BC2, Y_to_YTr_for_D2(:,2))

        ! For RHS3
        call D1c_MULTACC_3Dy(q3_y, tmp3_y, ysize_f(1),n1e,n2,ysize_f(3),n3e, diff21_coef, .true., NS_Q3_BC2, Yc_to_YcTr_for_D2(:,1))
        call D2c_MULTACC_3Dy(q3_y, tmp3_y, ysize_f(1),n1e,n2,ysize_f(3),n3e, diff22_coef, .true., NS_Q3_BC2, Yc_to_YcTr_for_D2(:,2))

        call D1s_MULTACC_3Dy(p12_y, tmp1_y, ysize_f(1),n1e,n2,ysize_f(3),n3e, conv2_coef, .false., NS_P12_BC2, Yc_to_YcTr_for_D1)
        call D0s_3Dy(q2_y, q2q2_y, ysize_f(1),n1e,n2,ysize_f(3),n3e, NS_Q2_BC2)
        q2q2_y=q2q2_y**2
        call D1ssh_MULTACC_3Dy(q2q2_y, tmp2_y, ysize_f(1),n1e,n2,ysize_f(3),n3e, conv2_coef, .true., NS_P22_BC2, Y_to_YTr_for_D1)

        call D1s_MULTACC_3Dy(p23_y, tmp3_y, ysize_f(1),n1e,n2,ysize_f(3),n3e, conv2_coef, .false., NS_P23_BC2, Yc_to_YcTr_for_D1)

        if(BC2==NOSLIP) call diff_at_wall_2 ! implicit: NS_Q3_BC2=Dirichlet

        RHS1_y = RHS1_y + tmp1_y
        RHS2_y = RHS2_y + tmp2_y
        RHS3_y = RHS3_y + tmp3_y

        ! write(*,*) 'RHS1_y = ', RHS1_y(:,20,32)

        ! _______________________________________________________________________________________________________________________
        ! ***********************************************************************************************************************
        ! ****************************** Derivative along the X direction *******************************************************
        ! ***********************************************************************************************************************

        call transpose_z_to_y(p13_z, p13_y, decomp_following)
        call transpose_y_to_x(p13_y, p13_x, decomp_following)

        call transpose_y_to_x(p12_y, p12_x, decomp_following)

        call transpose_y_to_x(RHS1_y, RHS1_x, decomp_following)
        call transpose_y_to_x(RHS2_y, RHS2_x, decomp_following)
        call transpose_y_to_x(RHS3_y, RHS3_x, decomp_following)

        tmp1_x = 0.d0
        tmp2_x = 0.d0
        tmp3_x = 0.d0

        n2e=xsize_f(2)!(min(n2m, xend(2))-xstart(2))+1
        n3e=xsize_f(3)!(min(n3m, xend(3))-xstart(3))+1

        ! For RHS1
        call D0s_3Dx(q1_x, q1q1_x, n1,xsize_f(2),n2e,xsize_f(3),n3e, NS_Q1_BC1)
        q1q1_x=q1q1_x**2
        call D1ssh_ACC_3Dx(q1q1_x, tmp1_x, n1,xsize_f(2),n2e,xsize_f(3),n3e, conv1_coef, .true., NS_P11_BC1)
        call D2c_ACC_3Dx(q1_x, tmp1_x, n1,xsize_f(2),n2e,xsize_f(3),n3e, diff1_coef, .false., NS_Q1_BC1)

        ! For RHS2
        call D1s_ACC_3Dx(p12_x, tmp2_x, n1,xsize_f(2),n2e,xsize_f(3),n3e, conv1_coef, .false., NS_P12_BC1)
        call D2c_ACC_3Dx(q2_x, tmp2_x, n1,xsize_f(2),n2e,xsize_f(3),n3e, diff1_coef, .true., NS_Q2_BC1)

        ! For RHS3
        call D1s_ACC_3Dx(p13_x, tmp3_x, n1,xsize_f(2),n2e,xsize_f(3),n3e, conv1_coef, .false., NS_P13_BC1)
        call D2c_ACC_3Dx(q3_x, tmp3_x, n1,xsize_f(2),n2e,xsize_f(3),n3e, diff1_coef, .true., NS_Q3_BC1)

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
            do k=xstart_f(3), min(n3m, xend_f(3))       !do k=1,n3m
                do j=xstart_f(2), min(n2m, xend_f(2))       !do j=1,n2m
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

            do k=ystart_f(3), min(n3m, yend_f(3))       !do k=1,n3m
                do i=ystart_f(1), min(n1m, yend_f(1))       !do i=1,n1m

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

            do j = zstart_f(2), min(n2m, zend_f(2))     !do j=1,n2m
                do i=zstart_f(1), min(n1m, zend_f(1))       !do i=1,n1m
                    tmp1_z(i,j,1)   =tmp1_z(i,j,1)  +( q1_z(i,j,2)       /3.d0 - q1_z(i,j,1)      + q1_wall30(i,j)*2.d0/3.d0 ) /diff3_coef**2
                    tmp1_z(i,j,n3m) =tmp1_z(i,j,n3m)+( q1_z(i,j,n3m-1)   /3.d0 - q1_z(i,j,n3m)    + q1_wall31(i,j)*2.d0/3.d0)  /diff3_coef**2

                    tmp2_z(i,j,1)   =tmp2_z(i,j,1)  +( q2_z(i,j,2)       /3.d0 - q2_z(i,j,1)      + q2_wall30(i,j)*2.d0/3.d0)  /diff3_coef**2
                    tmp2_z(i,j,n3m) =tmp2_z(i,j,n3m)+( q2_z(i,j,n3m-1)   /3.d0 - q2_z(i,j,n3m)    + q2_wall31(i,j)*2.d0/3.d0 ) /diff3_coef**2
                enddo
            enddo

        end subroutine diff_at_wall_3

    end subroutine


    subroutine perform_intermediate_velocity()

        implicit none

        real*8                                                                                  :: RHS1i, RHS2i, RHS3i
        integer                                                                                 :: i,j,k
        integer                                                                                 :: n1s, n1e, n2s,n2e, n3s,n3e

        n1e=min(n1-1, xend_f(1))
        n2e=min(n2-1, xend_f(2))
        n3e=min(n3-1, xend_f(3))

        if ((BC1==UNBOUNDED).or.(BC1==FRINGE)) n1s=xstart_f(1)
        if ((BC1/=UNBOUNDED).and.(BC1/=FRINGE)) n1s=max(2,xstart_f(1))
        n2s=xstart_f(2)
        n3s=xstart_f(3)
        do k = n3s, n3e
            do j = n2s, n2e
                do i = n1s, n1e
                    q1_x(i,j,k) = q1_save_x(i,j,k) + NSTERMS1(i,j,k)
                end do
            end do
        end do


        if (BC2==UNBOUNDED) n2s=xstart_f(2)
        if (BC2/=UNBOUNDED) n2s=max(2,xstart_f(2))
        ! On ne touche que q2[2..n1-2]
        if (BC1==OPEN) then
            n1s=2
            n1e=n1-2
        endif

        n3s=xstart_f(3)

        do k = n3s, n3e
            do j = n2s, n2e
                do i = n1s, n1e
                    q2_x(i,j,k) = q2_save_x(i,j,k) + NSTERMS2(i,j,k)
                end do
            end do
        end do

        if (BC3==UNBOUNDED) n3s=xstart_f(3)
        if (BC3/=UNBOUNDED) n3s=max(2,xstart_f(3))

        ! On ne touche que q3[2..n1-2]
        if (BC1==OPEN) then
            n1s=2
            n1e=n1-2
        endif

        n2s=xstart_f(2)

        do k = n3s, n3e
            do j = n2s, n2e
                do i = n1s, n1e
                    q3_x(i,j,k) = q3_save_x(i,j,k) + NSTERMS3(i,j,k)
                end do
            end do
        end do

        ! q3_x, q2_x and w_c contains u*,v* and w*

    end subroutine perform_intermediate_velocity

    ! Perform UV, UW, VW at the center of edges
    subroutine perform_velocity_products()
        use following_mesh

        implicit none

        integer k,j,i
        real*8, dimension(ystart_f(1):yend_f(1),ystart_f(2):yend_f(2),ystart_f(3):yend_f(3))   :: tmp12_y
        real*8, dimension(zstart_f(1):zend_f(1),zstart_f(2):zend_f(2),zstart_f(3):zend_f(3))   :: tmp23_z, tmp13_z

        integer :: n1e, n2e, n3e, n1s, n2s, n3s

        ! X ------------------------------------------------------------------------

        call D0ssh_3Dx(q3_x, p13_x, n1,xsize_f(2),xsize_f(2),xsize_f(3),xsize_f(3), NS_Q3_BC1)
        call D0ssh_3Dx(q2_x, p12_x, n1,xsize_f(2),xsize_f(2),xsize_f(3),xsize_f(3), NS_Q2_BC1)

        if ((BC1==NOSLIP).or.(BC1==FREESLIP).or.(BC1==PSEUDO_PERIODIC).or.(BC1==OPEN)) then

            do k=xstart_f(3), xend_f(3)       !do k=1,n3
                do j=xstart_f(2), xend_f(2)       !do j=1,n2
                    p13_x(1, j,  k)   =q3_wall10(j, k)
                    p13_x(n1, j, k)   =q3_wall11(j, k)

                    p12_x(1, j,  k)   =q2_wall10(j, k)
                    p12_x(n1, j, k)   =q2_wall11(j, k)
                enddo
            enddo

        end if

        ! Y ------------------------------------------------------------------------
        call transpose_x_to_y(p12_x, p12_y, decomp_following)

        call D0ssh_3Dy(q3_y, p23_y, ysize_f(1),ysize_f(1),n2,ysize_f(3),ysize_f(3), NS_Q3_BC2)
        ! call D0ssh_MULT_3Dy(q1_y, p12_y, ysize(1),ysize(1),n2,ysize(3),ysize(3), NS_Q1_BC2)
        call D0ssh_3Dy(q1_y, tmp12_y, ysize_f(1),ysize_f(1),n2,ysize_f(3),ysize_f(3), NS_Q1_BC2)

        if ((BC2==NOSLIP).or.(BC2==FREESLIP)) then

            do k=ystart_f(3), yend_f(3)      !do k=1,n3
                do i=ystart_f(1),yend_f(1)       !do i=1,n1
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
        call transpose_y_to_z(p23_y, p23_z, decomp_following)

        call transpose_x_to_y(p13_x, p13_y, decomp_following)
        call transpose_y_to_z(p13_y, p13_z, decomp_following)

        ! call D0ssh_MULT_3Dz(q2_z, p23_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_Q2_BC3)
        ! call D0ssh_MULT_3Dz(q1_z, p13_z, zsize(1),zsize(1),zsize(2),zsize(2),n3, NS_Q1_BC3)

        call D0ssh_3Dz(q2_z, tmp23_z, zsize_f(1),zsize_f(1),zsize_f(2),zsize_f(2),n3, NS_Q2_BC3)
        call D0ssh_3Dz(q1_z, tmp13_z, zsize_f(1),zsize_f(1),zsize_f(2),zsize_f(2),n3, NS_Q1_BC3)

        ! ATTENTION
        if ((BC3==NOSLIP).or.(BC3==FREESLIP)) then

            do i=zstart_f(1), zend_f(1)
                do j=zstart_f(2), zend_f(2)
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

            deallocate(p13_z)
            deallocate(p23_z)

            deallocate(p23_y)
            deallocate(p13_y)
            deallocate(p12_y)

            deallocate(p13_x)
            deallocate(p12_x)

        end subroutine deallocate_data

    end subroutine finalize

end module following_velocity_solver