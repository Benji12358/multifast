module IBM_data
    use decomp_2d, only:DECOMP_INFO

    implicit none
    real*8, dimension(:,:,:), allocatable   :: IBM_mask1_x, IBM_mask2_x, IBM_mask3_x
    real*8, dimension(:,:,:), allocatable   :: IBM_modulation_x, IBM_modulation_y, IBM_modulation_z
    real*8, dimension(:,:,:), allocatable   :: vel_term1, vel_term2, vel_term3
    real*8                                  :: q_bound(3)
    real*8									:: xs, xe, ys, ye, zs, ze
    real*8                                  :: ibm_volume
    real*8, dimension(:), allocatable       :: flow_rate_IBM, kin1_IBM,kin2_IBM,kin3_IBM
    real*8, dimension(:), allocatable       :: q1max_IBM,q2max_IBM,q3max_IBM
    real*8, dimension(:), allocatable       :: q1min_IBM,q2min_IBM,q3min_IBM, prmin_IBM, prmax_IBM

    ! 2nd order box
    TYPE(DECOMP_INFO)                       :: decomp_fine
    TYPE(DECOMP_INFO)                       :: subdecomp_fine_computational, subdecomp_fine_physical, subdecomp_coarse_computational
    integer                                 :: i_start_ibm_2nd_cp, i_end_ibm_2nd_cp, j_start_ibm_2nd_cp, j_end_ibm_2nd_cp, k_start_ibm_2nd_cp, k_end_ibm_2nd_cp
    integer                                 :: i_start_ibm_2nd_cc, i_end_ibm_2nd_cc, j_start_ibm_2nd_cc, j_end_ibm_2nd_cc, k_start_ibm_2nd_cc, k_end_ibm_2nd_cc
    integer                                 :: i_start_ibm_2nd_fc, i_end_ibm_2nd_fc, j_start_ibm_2nd_fc, j_end_ibm_2nd_fc, k_start_ibm_2nd_fc, k_end_ibm_2nd_fc
    integer                                 :: i_start_ibm_2nd_fp, i_end_ibm_2nd_fp, j_start_ibm_2nd_fp, j_end_ibm_2nd_fp, k_start_ibm_2nd_fp, k_end_ibm_2nd_fp
    integer                                 :: n1_ibm, n2_ibm, n3_ibm
    real*8, dimension(:,:,:), allocatable   :: vel_term1_ibm, vel_term2_ibm, vel_term3_ibm

    real*8, dimension(:,:,:), allocatable   :: IBM_mask1_x_ibm, IBM_mask2_x_ibm, IBM_mask3_x_ibm

    real*8,dimension (:,:,:), allocatable   :: q3_x_ibm, q3_y_ibm, q3_z_ibm
    real*8,dimension (:,:,:), allocatable   :: q2_x_ibm, q2_y_ibm, q2_z_ibm
    real*8,dimension (:,:,:), allocatable   :: q1_x_ibm, q1_y_ibm, q1_z_ibm
    real*8,dimension (:,:,:), allocatable   :: dp_x_ibm, dp_y_ibm, dp_z_ibm

    real*8,dimension (:,:), allocatable     :: q3_wall10_ibm, q2_wall10_ibm, q1_wall10_ibm
    real*8,dimension (:,:), allocatable     :: q3_wall11_ibm, q2_wall11_ibm, q1_wall11_ibm

    real*8,dimension (:,:), allocatable     :: q3_wall20_ibm, q2_wall20_ibm, q1_wall20_ibm
    real*8,dimension (:,:), allocatable     :: q3_wall21_ibm, q2_wall21_ibm, q1_wall21_ibm

    real*8,dimension (:,:), allocatable     :: q3_wall30_ibm, q2_wall30_ibm, q1_wall30_ibm
    real*8,dimension (:,:), allocatable     :: q3_wall31_ibm, q2_wall31_ibm, q1_wall31_ibm

    ! Mesh info global fine mesh
    integer, dimension(3)                   :: xstart_fine, xend_fine, xsize_fine  ! x-pencil
    integer, dimension(3)                   :: ystart_fine, yend_fine, ysize_fine  ! y-pencil
    integer, dimension(3)                   :: zstart_fine, zend_fine, zsize_fine  ! z-pencil

    ! Mesh info local coarse computational mesh
    integer, dimension(3)                   :: xstart_ibm_ccm, xend_ibm_ccm, xsize_ibm_ccm  ! x-pencil
    integer, dimension(3)                   :: ystart_ibm_ccm, yend_ibm_ccm, ysize_ibm_ccm  ! y-pencil
    integer, dimension(3)                   :: zstart_ibm_ccm, zend_ibm_ccm, zsize_ibm_ccm  ! z-pencil

    ! Mesh info local coarse physical mesh
    integer, dimension(3)                   :: xstart_ibm_cpm, xend_ibm_cpm, xsize_ibm_cpm  ! x-pencil
    integer, dimension(3)                   :: ystart_ibm_cpm, yend_ibm_cpm, ysize_ibm_cpm  ! y-pencil
    integer, dimension(3)                   :: zstart_ibm_cpm, zend_ibm_cpm, zsize_ibm_cpm  ! z-pencil

    ! Mesh info local fine computational mesh
    integer, dimension(3)                   :: xstart_ibm_fcm, xend_ibm_fcm, xsize_ibm_fcm  ! x-pencil
    integer, dimension(3)                   :: ystart_ibm_fcm, yend_ibm_fcm, ysize_ibm_fcm  ! y-pencil
    integer, dimension(3)                   :: zstart_ibm_fcm, zend_ibm_fcm, zsize_ibm_fcm  ! z-pencil

    ! Mesh info local fine physical mesh
    integer, dimension(3)                   :: xstart_ibm_fpm, xend_ibm_fpm, xsize_ibm_fpm  ! x-pencil
    integer, dimension(3)                   :: ystart_ibm_fpm, yend_ibm_fpm, ysize_ibm_fpm  ! y-pencil
    integer, dimension(3)                   :: zstart_ibm_fpm, zend_ibm_fpm, zsize_ibm_fpm  ! z-pencil

    logical                                 :: object_in_current_proc_ccm_x, object_in_current_proc_ccm_y, object_in_current_proc_ccm_z
    logical                                 :: object_in_current_proc_cpm_x, object_in_current_proc_cpm_y, object_in_current_proc_cpm_z
    logical                                 :: object_in_current_proc_fcm_x, object_in_current_proc_fcm_y, object_in_current_proc_fcm_z
    logical                                 :: object_in_current_proc_fpm_x, object_in_current_proc_fpm_y, object_in_current_proc_fpm_z

    real*8                                  :: dx1_ibm, dx2_ibm, dx3_ibm

    real*8, dimension(:), allocatable       :: X_ibm, Xc_ibm
    real*8, dimension(:), allocatable       :: Y_ibm, Yc_ibm
    real*8, dimension(:), allocatable       :: Z_ibm, Zc_ibm
    real*8, dimension(:), allocatable       :: cell_size_Y_ibm

    ! Irregular coef info
    real*8,dimension(:), allocatable        :: Y_to_YTr_for_D1_ibm
    real*8,dimension(:), allocatable        :: Yc_to_YcTr_for_D1_ibm

    real*8,dimension(:, :), allocatable     :: Y_to_YTr_for_D2_ibm
    real*8,dimension(:, :), allocatable     :: Yc_to_YcTr_for_D2_ibm

end module IBM_data
