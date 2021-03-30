module IBM_data
    implicit none
    real*8, dimension(:,:,:), allocatable   :: IBM_mask1_x, IBM_mask2_x, IBM_mask3_x
    real*8, dimension(:,:,:), allocatable   :: IBM_mask1_z, IBM_mask2_z, IBM_mask3_z
    real*8, dimension(:,:,:), allocatable   :: vel_term1, vel_term2, vel_term3
    real*8                                  :: q_bound(3)
    real*8, dimension(:), allocatable       :: flow_rate_IBM, kin1_IBM,kin2_IBM,kin3_IBM
    real*8, dimension(:), allocatable       :: q1max_IBM,q2max_IBM,q3max_IBM
    real*8, dimension(:), allocatable       :: q1min_IBM,q2min_IBM,q3min_IBM, prmin_IBM, prmax_IBM
end module IBM_data
