module embedded_scalar_solver

    use decomp_2d
    use embedded_mesh
    use boundaries
    use schemes3D_interface
    use embedded_irregular_derivative_coefficients
    use DNS_settings
    use embedded_scalar_data
    use embedded_physical_fields
    use embedded_data

    implicit none

    real*8, dimension(:,:,:), allocatable       :: previousRHS1, previousRHS2, previoussource_term, Fextsca, previousFextsca

    real*8  :: gam,rom

contains

    subroutine init()
        implicit none


        allocate(previousRHS1(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)))
        previousRHS1=0.d0
        allocate(previousRHS2(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)))
        previousRHS2=0.d0
        allocate(previoussource_term(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)))
        previoussource_term=0.d0

        allocate(Fextsca(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)))
        Fextsca=0.d0
        allocate(previousFextsca(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)))
        previousFextsca=0.d0

    end subroutine init

    subroutine solve_scalar(q1_x, q2_y, q3_z, ns)

        use time_schemes
        use run_ctxt_data
        use numerical_methods_settings

        use embedded_scalar_data

        implicit none

        integer, intent(in) :: ns

        real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)) :: q1_x
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)) :: q2_y, sca_y0, q_y
        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3)) :: q3_z

        real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)) :: q1S   ! scalar * Q1
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)) :: q2S   ! scalar * Q2
        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3)) :: q3S   ! scalar * Q3
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)) :: RHS
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)) :: RHS1
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)) :: RHS2, RHS2_1
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)) :: source_term

        call set_time_coeffs
        call perform_velocityscalar_products(q1_x, q2_y, q3_z, q1S, q2S, q3S)
        call perform_RHS1(q1S, q2S, q3S, RHS1)
        call perform_RHS2(RHS2)

        ! if (init_type==KAWAMURA_INIT) then

        !     if (streamwise==1) then

        !         call transpose_x_to_y(q1_x, q_y)

        !     elseif (streamwise==3) then

        !         call transpose_x_to_y(q3_x, q_y)

        !     endif

        !     call perform_source_term(q_y, source_term)

        ! else

        !     source_term = 0.d0

        ! endif

        ! for now, the set of equations solved are the ones of LYONS (1991)
        source_term = 0.d0

        ! RHS=RHS1+RHS2

        call next_scalar

        ! if (init_type==CONSTANT_HEAT_FLUX) call next_wall_scalar

        previousRHS1=RHS1
        previousRHS2=RHS2
        previoussource_term=source_term
        previousFextsca=Fextsca

        call transpose_y_to_x(sca_y, sca_x, decomp_embedded)
        call transpose_y_to_z(sca_y, sca_z, decomp_embedded)

    contains

        subroutine set_time_coeffs()

            logical, save :: first_time_advancement=.true.

!            previousRHS_are_available=.false.

            if (first_time_advancement) then

                gam=1.d0
                rom=0.d0

            else
                gam=ga(ns)
                rom=ro(ns)

            end if

            first_time_advancement=.false.

        end subroutine set_time_coeffs

        subroutine next_scalar

            implicit none
            integer :: i,j,k
            integer :: n1s, n1e, n2s,n2e, n3s,n3e
            real*8, dimension(xstart_e(1):xend_e(1),xstart_e(2):xend_e(2),xstart_e(3):xend_e(3)) :: sca_term_x
            real*8, dimension(ystart_e(1):yend_e(1),ystart_e(2):yend_e(2),ystart_e(3):yend_e(3)) :: sca_term

            n1s=ystart_e(1)
            n1e=min(n1-1, yend_e(1))
            n2s=1
            n2e=n2-1

            n3s=ystart_e(3)
            n3e=min(n3-1, yend_e(3))

            do k = n3s, n3e
                do j = n2s, n2e
                    do i = n1s, n1e
                        sca_y(i,j,k)=sca_y(i,j,k)+dt*gam*(RHS1(i,j,k)+RHS2(i,j,k)+source_term(i,j,k)+Fextsca(i,j,k))+dt*rom*(previousRHS1(i,j,k)+previousRHS2(i,j,k)+previoussource_term(i,j,k)+previousFextsca(i,j,k))

                    end do
                end do
            end do

        end subroutine next_scalar

    end subroutine solve_scalar

    subroutine perform_RHS1(q1S, q2S, q3S, RHS1)

        use embedded_scalar_data

        implicit none

        real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)), intent (in) :: q1S   ! scalar * Q1
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)), intent (in) :: q2S   ! scalar * Q2
        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3)), intent (in) :: q3S   ! scalar * Q3

        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)), intent (out)    :: RHS1

        real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3))  :: RHS1x   ! scalar * Q1
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3))  :: RHS1y   ! scalar * Q2
        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3))  :: RHS1z   ! scalar * Q3

        integer :: n1e, n2e, n3e

        RHS1x=0.d0
        RHS1y=0.d0
        RHS1z=0.d0
        RHS1=0.d0

        ! CONSERVATIVE


        ! Q3*scalar ************************************************
        ! velocity interpolation: j (1:n2-1), k(1,n3-1))
        n1e=(min(n1-1, zend_e(1))-zstart_e(1))+1
        n2e=(min(n2-1, zend_e(2))-zstart_e(2))+1
        call D1c_3Dz(q3S, RHS1z, zsize_e(1),n1e,zsize_e(2),n2e,n3, -dx3, .true., TRANSPORT_Q3S_BC3)


        ! Q2*scalar ************************************************
        ! velocity interpolation: i (1:n1-1), k(1,n3-1)
        n1e=(min(n1-1, yend_e(1))-ystart_e(1))+1
        n3e=(min(n3-1, yend_e(3))-ystart_e(3))+1

        call D1c_MULTACC_3Dy(q2S, RHS1y, ysize_e(1),n1e,n2,ysize_e(3),n3e, -dx2,.true., TRANSPORT_Q2S_BC2, Yc_to_YcTr_for_D1)

        !if (SCA_BC2==FIXED_VALUE) call d1c_wall2
        call d1c_wall2

        ! Q1*scalar ************************************************
        ! velocity interpolation: i (1:n1-1), j(1,n2-1)
        n2e=(min(n2m, xend_e(2))-xstart_e(2))+1
        n3e=(min(n3m, xend_e(3))-xstart_e(3))+1
        call D1c_3Dx(q1S, RHS1x, n1,xsize_e(2),n2e,xsize_e(3),n3e, -dx1, .true., TRANSPORT_Q1S_BC1)

        ! RHS1 et RHS2 Table compilation
        RHS1=RHS1y
        call transpose_x_to_y(RHS1x, RHS1y, decomp_embedded)
        RHS1=RHS1+RHS1y
        call transpose_z_to_y(RHS1z, RHS1y, decomp_embedded)
        RHS1=RHS1+RHS1y

    contains

        subroutine d1c_wall2()
            use boundary_scheme
            implicit none
            integer ::i,j,k

            if (SCA_BC2==FIXED_VALUE) then

                do k=ystart_e(3), min(n3m, yend_e(3))       !do k=1,n3m
                    do i=ystart_e(1), min(n1m, yend_e(1))       !do i=1,n1m
                        ! RHS1y(i,1,k)=-(q2S(i,1,k)-q2_wall20(i,k)*sca_wall20(i,k))/(Yc(1)-Y(1))
                        ! RHS1y(i,n2m,k)=-(q2_wall21(i,k)*sca_wall21(i,k)-q2S(i,n2m,k))/(Y(n2)-Yc(n2-1))

                        RHS1y(i,1,k)= - ( ( 0.5d0 * ( q2S(i,1,k) + q2S(i,2,k) ) - q2_wall20(i,k)*sca_wall20(i,k) )/dx2 ) * Yc_to_YcTr_for_D1(1)
                        RHS1y(i,n2m,k)= - ( ( q2_wall21(i,k)*sca_wall21(i,k) - 0.5d0 * ( q2S(i,n2-1,k) + q2S(i,n2-2,k) ) )/dx2 ) * Yc_to_YcTr_for_D1(n2m)

                    enddo
                enddo


            end if
            
            if (SCA_BC2==FIXED_FLUX) then
                do k=ystart_e(3), min(n3m, yend_e(3))       !do k=1,n3m
                    do i=ystart_e(1), min(n1m, yend_e(1))       !do i=1,n1m
                        RHS1y(i,1,k)=-(q2S(i,1,k)-q2_wall20(i,k)*(heat_flux*(Yc(1)-Y(1))+sca_y(i,1,k)))/(Yc(1)-Y(1))
                        !RHS1y(i,n2m,k)=-(q2_wall21(i,k)*(heat_flux*(Yc(1)-Y(1))+sca_y(i,1,k))-q2S(i,n2m,k))/(Y(n2)-Yc(n2-1))
                        RHS1y(i,n2m,k)=-q2_wall21(i,k)*heat_flux
                    enddo
                enddo
            end if

        end subroutine d1c_wall2

    end subroutine perform_RHS1

    subroutine perform_RHS2(RHS2)
        implicit none

        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)), intent (out)    :: RHS2

        real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3))  :: RHS2x   ! scalar * Q1
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3))  :: RHS2y   ! scalar * Q2
        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3))  :: RHS2z   ! scalar * Q3

        integer :: n1e, n2e, n3e

        RHS2x=0.d0
        RHS2y=0.d0
        RHS2z=0.d0
        RHS2=0.d0

        ! Diffusive scalar - Z direction ************************************************
        n1e=(min(n1-1, zend_e(1))-zstart_e(1))+1
        n2e=(min(n2-1, zend_e(2))-zstart_e(2))+1
        call D2c_3Dz(sca_z, RHS2z, zsize_e(1),n1e,zsize_e(2),n2e,n3, dx3*dsqrt(ren*prandtl), .true., TRANSPORT_SCA_BC3)


        ! Diffusive scalar - Y direction ************************************************
        n1e=(min(n1-1, yend_e(1))-ystart_e(1))+1
        n3e=(min(n3-1, yend_e(3))-ystart_e(3))+1


        call D1c_MULTACC_3Dy(sca_y, RHS2y, ysize_e(1),n1e,n2,ysize_e(3),n3e, dx2*ren*prandtl,.true., TRANSPORT_SCA_BC2, Yc_to_YcTr_for_D2(:,1))
        call D2c_MULTACC_3Dy(sca_y, RHS2y, ysize_e(1),n1e,n2,ysize_e(3),n3e, dx2*dsqrt(ren*prandtl),.true., TRANSPORT_SCA_BC2, Yc_to_YcTr_for_D2(:,2))

        !if (SCA_BC2==FIXED_VALUE) call d2c_wall2
        call d2c_wall2

        ! Diffusive scalar - X direction ************************************************
        n2e=(min(n2m, xend_e(2))-xstart_e(2))+1
        n3e=(min(n3m, xend_e(3))-xstart_e(3))+1
        call D2c_3Dx(sca_x, RHS2x, n1,xsize_e(2),n2e,xsize_e(3),n3e, dx1*dsqrt(ren*prandtl), .true., TRANSPORT_SCA_BC1)

        ! RHS2 Table compilation

        RHS2=RHS2y
        call transpose_x_to_y(RHS2x, RHS2y, decomp_embedded)
        RHS2=RHS2+RHS2y
        call transpose_z_to_y(RHS2z, RHS2y, decomp_embedded)
        RHS2=RHS2+RHS2y

    contains

        subroutine d2c_wall2()
            use boundary_scheme
            implicit none
            integer ::i,j,k

            if (SCA_BC2==FIXED_VALUE) then
                do k=ystart_e(3), min(n3m, yend_e(3))       !do k=1,n3m
                    do i=ystart_e(1), min(n1m, yend_e(1))       !do i=1,n1m
                        RHS2y(i,1,k)=RHS2y(i,1,k)+( sca_y(i,2,k)*a3_d + sca_y(i,1,k)*a2_d + sca_wall20(i,k)*a1_d)            /(dx2*dsqrt(ren*prandtl))**2
                        RHS2y(i,n2m,k)=RHS2y(i,n2m,k)+( sca_y(i,n2m-1,k)*a1_u + sca_y(i,n2m,k)*a2_u + sca_wall21(i,k)*a3_u ) /(dx2*dsqrt(ren*prandtl))**2
                    enddo
                enddo
            end if

            if (SCA_BC2==FIXED_FLUX) then
                do k=ystart_e(3), min(n3m, yend_e(3))       !do k=1,n3m
                    do i=ystart_e(1), min(n1m, yend_e(1))       !do i=1,n1m
                        RHS2y(i,1,k)=RHS2y(i,1,k)+( sca_y(i,2,k)*a3_d + sca_y(i,1,k)*(a2_d+a1_d) +heat_flux*(Yc(1)-Y(1))*a1_d)            /(dx2*dsqrt(ren*prandtl))**2
                        RHS2y(i,n2m,k)=RHS2y(i,n2m,k)+( sca_y(i,n2m-1,k)*a1_u + sca_y(i,n2m,k)*(a2_u+a3_u) + heat_flux*(Y(n2)-Yc(n2-1))*a3_u ) /(dx2*dsqrt(ren*prandtl))**2
                    enddo
                enddo
            end if

        end subroutine d2c_wall2

    end subroutine perform_RHS2

    ! subroutine perform_source_term(q_y, source_term)

    !     use mpi
    !     use COMMON_workspace_view, only: COMMON_log_path
    !     use run_ctxt_data, only: ntime

    !     implicit none

    !     real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent (in) :: q_y
    !     real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: q_x
    !     real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent (out) :: source_term
    !     real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: source_term_x
    !     ! real*8 :: ut1, ut1_global
    !     real*8, dimension(n1) :: ut1, ut1_global
    !     real*8 :: umean1, umean_glob1, ren_tau
    !     integer:: i,j,k
    !     integer     :: mpi_err
    !     character*200       :: temperature_history_file

    !     call transpose_y_to_x(q_y,q_x)

    !     ! Perform the integral over y-direction
    !     ut1=0.d0
    !     do i = xstart(1),min(xend(1),n1-1)
    !         do j=xstart(2),min(xend(2),n2-1)
    !             do k=xstart(3),min(xend(3),n3-1)
    !                 ut1(i)=ut1(i)+q_x(i,j,k)*(Y(j+1)-Y(j))
    !             enddo
    !         enddo
    !     enddo

    !     call MPI_ALLREDUCE(ut1, ut1_global, n1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
    !     ut1_global=ut1_global/n3m
    !     ut1_global=ut1_global/L2

    !     do i = xstart(1),min(xend(1),n1-1)
    !         do j=xstart(2),min(xend(2),n2-1)
    !             do k=xstart(3),min(xend(3),n3-1)

    !                 source_term_x(i,j,k) = q_x(i,j,k) / ut1_global(i)

    !             enddo
    !         enddo
    !     enddo

    !     call transpose_x_to_y(source_term_x,source_term)

    !     ! compute u_tau
    !     umean1 = 0.d0

    !     do k=ystart(3), min(n3m, yend(3))
    !         do i=ystart(1), min(n1m, yend(1))
    !             umean1=umean1+q_y(i,1,k)
    !             umean1=umean1+q_y(i,n2-1,k)
    !         enddo
    !     enddo

    !     call MPI_ALLREDUCE (umean1, umean_glob1, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

    !     umean_glob1=umean_glob1/dfloat(2*n1m*n3m)
    !     ren_tau=dsqrt(ren*(umean_glob1/Yc(1)))

    !     source_term = source_term * ren_tau/ren

    ! end subroutine perform_source_term

    subroutine perform_velocityscalar_products(q1_x, q2_y, q3_z, q1S, q2S, q3S)

        use embedded_scalar_data

        implicit none
        real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)), intent (in) :: q1_x
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)), intent (in) :: q2_y
        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3)), intent (in) :: q3_z

        real*8, dimension(xstart_e(1):xend_e(1), xstart_e(2):xend_e(2), xstart_e(3):xend_e(3)), intent (out)    :: q1S   ! scalar * Q1
        real*8, dimension(ystart_e(1):yend_e(1), ystart_e(2):yend_e(2), ystart_e(3):yend_e(3)), intent (out)    :: q2S   ! scalar * Q2
        real*8, dimension(zstart_e(1):zend_e(1), zstart_e(2):zend_e(2), zstart_e(3):zend_e(3)), intent (out)    :: q3S   ! scalar * Q3

        integer :: n1e, n2e, n3e
        integer :: i,j,k

        ! NON CONSERVATIVE


        ! Q1*scalar ************************************************
        ! velocity interpolation: j (1:n2-1), k(1,n3-1)
        n2e=(min(n2-1, xend_e(2))-xstart_e(2))+1
        n3e=(min(n3-1, xend_e(3))-xstart_e(3))+1
        q1S=0.d0
        call D0s_3Dx(q1_x, q1S, n1,xsize_e(2),n2e,xsize_e(3),n3e, NS_Q1_BC1)
        q1S=q1S*sca_x

        ! Q2*scalar ************************************************
        ! velocity interpolation: i (1:n1-1), k(1,n3-1)
        n1e=(min(n1-1, yend_e(1))-ystart_e(1))+1
        n3e=(min(n3-1, yend_e(3))-ystart_e(3))+1
        q2S=0.d0
        call D0s_3Dy(q2_y, q2S, ysize_e(1),n1e,n2,ysize_e(3),n3e, NS_Q2_BC2)
        q2S=q2S*sca_y

        ! Q3*scalar ************************************************
        ! velocity interpolation: i (1:n1-1), j(1,n2-1)
        n1e=(min(n1-1, zend_e(1))-zstart_e(1))+1
        n2e=(min(n2-1, zend_e(2))-zstart_e(2))+1
        q3S=0.d0
        call D0s_3Dz(q3_z, q3S, zsize_e(1),n1e,zsize_e(2),n2e,n3, NS_Q2_BC3)
        q3S=q3S*sca_z

    end subroutine perform_velocityscalar_products!
!
end module embedded_scalar_solver
