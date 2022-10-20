
! module VELOCITY_properties
!     use mpi
!     use decomp_2d

!     use mesh
!     use boundaries
!     use irregular_derivative_coefficients
!     use DNS_settings
!     use schemes3D_interface

!     use VELOCITY_operations

!     implicit none

! contains

!     subroutine perform_kinetic(q3_z, q2_y, q1_x, flow_rate, spanwise_flow_rate, kin_energy, enstrophy, streamwise)

!         use physical_fields, only:q3_y
!         implicit none

!         real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in)    :: q3_z
!         real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)    :: q2_y
!         real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: q1_x
!         real*8, intent(out)     :: kin_energy(3), enstrophy, flow_rate, spanwise_flow_rate
!         integer, intent(in)     :: streamwise

!         real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))      :: wc_x, vc_x, uc_x
!         real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))      :: uc_y, vc_y, wc_y
!         real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))      :: uc_z, vc_z, wc_z

!         real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))      :: vortX_y, vortY_y, vortZ_y


!         integer k, j, i, mpi_err
!         real*8  :: cflmax_glob=0.d0, kin_energy_glob(3)=0.d0, enstrophy_glob=0.d0, flow_rate_glob=0.d0, spanwise_flow_rate_glob=0.d0

!         call perform_velocity_at_center(q3_z, q2_y, q1_x, uc_z, vc_y, wc_x)

!         call transpose_x_to_y(wc_x, wc_y)
!         call transpose_z_to_y(uc_z, uc_y)


!         call perform_vorticity(wc_y, vc_y, uc_y, vortX_y, vortY_y, vortZ_y)
! !        call test_vorticity

!         kin_energy=0.d0
!         enstrophy=0.d0
!         flow_rate=0.d0
!         spanwise_flow_rate=0.d0

!         do k=ystart(3), min(n3m, yend(3))
!             !do i=1,n1m
!             do j=1, n2m
!                 do i=ystart(1), min(n1m, yend(1))

!                     if (streamwise==1) then
!                         flow_rate= flow_rate+(wc_y(i,j,k)) * cell_size_Y(j) * dx1 * dx3
!                         spanwise_flow_rate=spanwise_flow_rate+(uc_y(i,j,k)) * cell_size_Y(j) * dx1 * dx3
!                     elseif (streamwise==3) then
!                         flow_rate= flow_rate+(uc_y(i,j,k)) * cell_size_Y(j) * dx1 * dx3
!                         spanwise_flow_rate=spanwise_flow_rate+(wc_y(i,j,k)) * cell_size_Y(j) * dx1 * dx3
!                     endif

!                     kin_energy(1)= kin_energy(1)+(wc_y(i,j,k)**2) * cell_size_Y(j) * dx1 * dx3
!                     kin_energy(2)= kin_energy(2)+(vc_y(i,j,k)**2) * cell_size_Y(j) * dx1 * dx3
!                     kin_energy(3)= kin_energy(3)+(uc_y(i,j,k)**2) * cell_size_Y(j) * dx1 * dx3
!                     enstrophy= enstrophy+(vortX_y(i,j,k)**2 + vortY_y(i,j,k)**2 + vortZ_y(i,j,k)**2) * cell_size_Y(j) * dx1 *dx3
!                 enddo
!             enddo
!         enddo


!         call MPI_ALLREDUCE (enstrophy, enstrophy_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
!         call MPI_ALLREDUCE (kin_energy, kin_energy_glob, 3, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
!         call MPI_ALLREDUCE (flow_rate, flow_rate_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
!         call MPI_ALLREDUCE (spanwise_flow_rate, spanwise_flow_rate_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

!         kin_energy=kin_energy_glob/(L1*L2*L3)
!         enstrophy=enstrophy_glob/(L1*L2*L3)
!         flow_rate=flow_rate_glob/(L1*L2*L3)
!         spanwise_flow_rate=spanwise_flow_rate_glob/(L1*L2*L3)

!         return

!     contains


!         subroutine test_vorticity()

!             use mathematical_constants
!             use snapshot_writer
!             use COMMON_workspace_view

!             implicit none

!             real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: q1c_x, q2c_x, q3c_x
!             real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))  :: q2c_y, q1c_y, q3c_y
!             real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))  :: q3c_z

!             real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: q2c_2Dx, q1c_2Dy
!             real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: q1c_2Dz, q3c_2Dx
!             real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: q3c_2Dy, q2c_2Dz

!             real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: wc_x, vc_x, uc_x
!             real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))  :: wc_y, vc_y, uc_y
!             real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))  :: wc_z, vc_z, uc_z


!             real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: vortX_x, vortY_x, vortZ_x
!             real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: vortX_x_exp, vortY_x_exp, vortZ_x_exp
!             real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))  :: errX_x, errY_x, errZ_x

!             integer                                                                     :: i, j, k

!             real*8          ::  vortX_x_glob=0.d0, vortY_x_glob=0.d0, vortZ_x_glob=0.d0
!             real*8          ::  vortX_x_exp_glob=0.d0, vortY_x_exp_glob=0.d0, vortZ_x_exp_glob=0.d0
!             integer         ::  mpi_err

!             q1c_x=0.d0

!             q1c_y=0.d0
!             q2c_y=0.d0
!             q3c_y=0.d0

!             q3c_z=0.d0

!             wc_y=0.d0
!             vc_y=0.d0
!             uc_y=0.d0

!             vortX_x=0.d0
!             vortY_x=0.d0
!             vortZ_x=0.d0

!             q1c_2Dy=0.d0
!             q1c_2Dz=0.d0
!             q2c_2Dx=0.d0
!             q2c_2Dz=0.d0
!             q3c_2Dx=0.d0
!             q3c_2Dy=0.d0

!             vortX_x_exp=0.d0
!             vortY_x_exp=0.d0
!             vortZ_x_exp=0.d0

!             write(*,*)'n3m', n1m,n2m,n3m

!         do i = 1, n1m
!             do j = xstart(2), min(xend(2),n2m)
!                 do k = xstart(3), min(xend(3),n3m)
!                     q1c_x(i, j, k)=dcos(Yc(j)*2*PI/L2)*dcos(Xc(k)*2.d0*PI/L3)
!                     q2c_x(i, j, k)=dcos(Zc(i)*2*PI/L1)*dcos(Xc(k)*2.d0*PI/L3)
!                     q3c_x(i, j, k)=dcos(Zc(i)*2*PI/L1)*dcos(Yc(j)*2.d0*PI/L2)

!                 end do
!             end do
!         end do

!         wc_x=q1c_x
!         vc_x=q2c_x
!         uc_x=q3c_x

!         call perform_vorticity_x(wc_x, vc_x, uc_x, vortZ_x, vortY_x, vortX_x)

!         do i = 1, n1m
!             do j = xstart(2), min(xend(2),n2m)
!                 do k = xstart(3), min(xend(3),n3m)

!                     q2c_2Dx(i, j, k)=-(2*PI/L1)*dsin(Zc(i)*2*PI/L1)*dcos(Xc(k)*2.d0*PI/L3)
!                     q1c_2Dy(i, j, k)=-(2*PI/L2)*dsin(Yc(j)*2*PI/L2)*dcos(Xc(k)*2.d0*PI/L3)

!                     q1c_2Dz(i, j, k)=dcos(Yc(j)*2.d0*PI/L2)*(-(2*PI/L3)*dsin(Xc(k)*2*PI/L3))
!                     q3c_2Dx(i, j, k)=-(2*PI/L1)*dsin(Zc(i)*2*PI/L1)*dcos(Yc(j)*2.d0*PI/L2)

!                     q2c_2Dz(i, j, k)=dcos(Zc(i)*2.d0*PI/L1)*(-(2*PI/L3)*dsin(Xc(k)*2*PI/L3))
!                     q3c_2Dy(i, j, k)=dcos(Zc(i)*2.d0*PI/L1)*(-(2*PI/L2)*dsin(Yc(j)*2*PI/L2))
!                 end do
!             end do
!         end do

!         vortZ_x_exp = q3c_2Dy - q2c_2Dz
!         vortY_x_exp = q1c_2Dz - q3c_2Dx
!         vortX_x_exp = q2c_2Dx - q1c_2Dy

!         errZ_x = vortZ_x_exp-vortZ_x
!         errY_x = vortY_x_exp-vortY_x
!         errX_x = vortX_x_exp-vortX_x

!         call create_snapshot(COMMON_snapshot_path, "IBM", q1c_x(:,:,:), "q1c_x", 1)
!         call create_snapshot(COMMON_snapshot_path, "IBM", q2c_x(:,:,:), "q2c_x", 1)
!         call create_snapshot(COMMON_snapshot_path, "IBM", q3c_x(:,:,:), "q3c_x", 1)

!         call create_snapshot(COMMON_snapshot_path, "IBM", q1c_2Dy(:,:,:), "q1c_2Dy", 1)
!         call create_snapshot(COMMON_snapshot_path, "IBM", q1c_2Dz(:,:,:), "q1c_2Dz", 1)
!         call create_snapshot(COMMON_snapshot_path, "IBM", q2c_2Dx(:,:,:), "q2c_2Dx", 1)
!         call create_snapshot(COMMON_snapshot_path, "IBM", q2c_2Dz(:,:,:), "q2c_2Dz", 1)
!         call create_snapshot(COMMON_snapshot_path, "IBM", q3c_2Dx(:,:,:), "q3c_2Dx", 1)
!         call create_snapshot(COMMON_snapshot_path, "IBM", q3c_2Dy(:,:,:), "q3c_2Dy", 1)

!         call create_snapshot(COMMON_snapshot_path, "IBM", errZ_x(:,:,:), "errZ_x", 1)
!         call create_snapshot(COMMON_snapshot_path, "IBM", errY_x(:,:,:), "errY_x", 1)
!         call create_snapshot(COMMON_snapshot_path, "IBM", errX_x(:,:,:), "errX_x", 1)

!         end subroutine test_vorticity

!     end subroutine



!     subroutine perform_kinetic_IBM(q3c_z, q2c_y, q1c_x, flow_rate, k1, k2, k3)
!         implicit none
!         real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)      :: q1c_x
!         real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)      :: q2c_y
!         real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in)      :: q3c_z


!         real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))      :: q2c_x, q3c_x
!         real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))      :: q1c_y, q3c_y
!         real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))      :: q1c_z, q2c_z

!         real*8, dimension(xsize(1))                                                     :: flow_rate, k1,k2,k3

!         real*8  :: debit, debit_loc, k1_loc, k2_loc, k3_loc
!         integer k,i,j, mpi_err


!         call transpose_z_to_y(q3c_z, q3c_y)
!         call transpose_y_to_x(q3c_y, q3c_x)

!         call transpose_y_to_x(q2c_y, q2c_x)

!         do i=xstart(1), min(n1m, xend(1))

!             ! Flow rate calculation
!             debit_loc=0.d0
!             k1_loc=0.d0
!             k2_loc=0.d0
!             k3_loc=0.d0


!             do k=xstart(3), min(n3m, xend(3))
!                 do j=xstart(2), min(n2m, xend(2))
!                     debit_loc=debit_loc +q1c_x(i,j,k) * (Y(j+1)-Y(j))*dx3
!                     k1_loc=k1_loc       +(q1c_x(i,j,k)**2) * (Y(j+1)-Y(j)) * dx3
!                     k2_loc=k2_loc       +(q2c_x(i,j,k)**2) * (Y(j+1)-Y(j)) * dx3
!                     k3_loc=k3_loc       +(q3c_x(i,j,k)**2) * (Y(j+1)-Y(j)) * dx3
!                 enddo
!             enddo

!             call MPI_ALLREDUCE (debit_loc, flow_rate(i), 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
!             call MPI_ALLREDUCE (k1_loc, k1(i), 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
!             call MPI_ALLREDUCE (k2_loc, k2(i), 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
!             call MPI_ALLREDUCE (k3_loc, k3(i), 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

!         enddo

!         flow_rate=flow_rate/(L2*L3)
!         k1=k1/(L2*L3)
!         k2=k2/(L2*L3)
!         k3=k3/(L2*L3)

!     end subroutine perform_kinetic_IBM

!     subroutine perform_maxvel_IBM(q3_x, q2_x, q1_x, dp_x, q1_max, q2_max, q3_max, pr_max, q1_min, q2_min, q3_min, pr_min)
!         implicit none
!         real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)      :: q1_x, q2_x, q3_x, dp_x


!         real*8, dimension(xsize(1))                                                     :: q1_max, q2_max, q3_max, pr_max
!         real*8, dimension(xsize(1))                                                     :: q1_min, q2_min, q3_min, pr_min

!         real*8  :: debit, debit_loc, q1_max_loc, q2_max_loc, q3_max_loc, q1_min_loc, q2_min_loc, q3_min_loc, pr_min_loc, pr_max_loc
!         integer k,i,j, mpi_err


!         do i=xstart(1), min(n1m, xend(1))

!             ! Flow rate calculation
!             debit_loc=0.d0
!             q1_max_loc=-10000.d0
!             q2_max_loc=-10000.d0
!             q3_max_loc=-10000.d0
!             pr_max_loc=-10000.d0
!             q1_min_loc=10000.d0
!             q2_min_loc=10000.d0
!             q3_min_loc=10000.d0
!             pr_min_loc=10000.d0


!             do k=xstart(3), min(n3m, xend(3))
!                 do j=xstart(2), min(n2m, xend(2))

!                     q1_max_loc=max(q1_max_loc, q1_x(i,j,k))
!                     q2_max_loc=max(q2_max_loc, q2_x(i,j,k))
!                     q3_max_loc=max(q3_max_loc, q3_x(i,j,k))
!                     pr_max_loc=max(pr_max_loc, dp_x(i,j,k))

!                     q1_min_loc=min(q1_min_loc, q1_x(i,j,k))
!                     q2_min_loc=min(q2_min_loc, q2_x(i,j,k))
!                     q3_min_loc=min(q3_min_loc, q3_x(i,j,k))
!                     pr_min_loc=min(pr_min_loc, dp_x(i,j,k))

!                 enddo
!             enddo

!             call MPI_ALLREDUCE (q1_max_loc, q1_max(i), 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
!             call MPI_ALLREDUCE (q2_max_loc, q2_max(i), 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
!             call MPI_ALLREDUCE (q3_max_loc, q3_max(i), 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)
!             call MPI_ALLREDUCE (pr_max_loc, pr_max(i), 1, MPI_DOUBLE_PRECISION , MPI_MAX , MPI_COMM_WORLD , mpi_err)

!             call MPI_ALLREDUCE (q1_min_loc, q1_min(i), 1, MPI_DOUBLE_PRECISION , MPI_MIN , MPI_COMM_WORLD , mpi_err)
!             call MPI_ALLREDUCE (q2_min_loc, q2_min(i), 1, MPI_DOUBLE_PRECISION , MPI_MIN , MPI_COMM_WORLD , mpi_err)
!             call MPI_ALLREDUCE (q3_min_loc, q3_min(i), 1, MPI_DOUBLE_PRECISION , MPI_MIN , MPI_COMM_WORLD , mpi_err)
!             call MPI_ALLREDUCE (pr_min_loc, pr_min(i), 1, MPI_DOUBLE_PRECISION , MPI_MIN , MPI_COMM_WORLD , mpi_err)

!         enddo

!     end subroutine perform_maxvel_IBM

! end module VELOCITY_properties
