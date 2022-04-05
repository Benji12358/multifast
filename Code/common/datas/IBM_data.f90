module IBM_data
    implicit none
    real*8, dimension(:,:,:), allocatable   :: IBM_mask1_x, IBM_mask2_x, IBM_mask3_x
    ! real*8, dimension(:,:,:), allocatable   :: IBM_mask1_z, IBM_mask2_z, IBM_mask3_z
    real*8, dimension(:,:,:), allocatable   :: vel_term1, vel_term2, vel_term3
    real*8                                  :: q_bound(3)
    real*8, dimension(:), allocatable       :: flow_rate_IBM, kin1_IBM,kin2_IBM,kin3_IBM
    real*8, dimension(:), allocatable       :: q1max_IBM,q2max_IBM,q3max_IBM
    real*8, dimension(:), allocatable       :: q1min_IBM,q2min_IBM,q3min_IBM, prmin_IBM, prmax_IBM

    ! 2nd order box
    real*8,dimension (:,:,:), allocatable       :: q3_x_ibm, q3_y_ibm, q3_z_ibm
    real*8,dimension (:,:,:), allocatable       :: q2_x_ibm, q2_y_ibm, q2_z_ibm
    real*8,dimension (:,:,:), allocatable       :: q1_x_ibm, q1_y_ibm, q1_z_ibm
    real*8,dimension (:,:,:), allocatable       :: dp_x_ibm, dp_y_ibm, dp_z_ibm

    real*8,dimension (:,:), allocatable         :: q3_wall10_ibm, q2_wall10_ibm, q1_wall10_ibm
    real*8,dimension (:,:), allocatable         :: q3_wall11_ibm, q2_wall11_ibm, q1_wall11_ibm

    real*8,dimension (:,:), allocatable         :: q3_wall20_ibm, q2_wall20_ibm, q1_wall20_ibm
    real*8,dimension (:,:), allocatable         :: q3_wall21_ibm, q2_wall21_ibm, q1_wall21_ibm

    real*8,dimension (:,:), allocatable         :: q3_wall30_ibm, q2_wall30_ibm, q1_wall30_ibm
    real*8,dimension (:,:), allocatable         :: q3_wall31_ibm, q2_wall31_ibm, q1_wall31_ibm

    real*8,dimension (:,:,:), allocatable       :: divu_z_ibm
    real*8,dimension (:,:,:), allocatable       :: dphidx1_x_ibm, dphidx2_x_ibm, dphidx3_x_ibm
    real*8,dimension (:,:,:), allocatable       :: dphidx1_y_ibm, dphidx2_y_ibm, dphidx3_y_ibm
    real*8,dimension (:,:,:), allocatable       :: dphidx1_z_ibm, dphidx2_z_ibm, dphidx3_z_ibm
end module IBM_data
