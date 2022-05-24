
module IBM_settings
    implicit none
    integer, parameter                  :: IBM_INTERPOL_NONE=0, IBM_FADLUN=1, ANTISYMMETRIC_INTERPOL=2, SECOND_ORDER_INTERPOL=3
    real*8, dimension(:), allocatable   :: body_x1, body_x2, body_x3
    real*8, dimension(:), allocatable   :: body_scale_x1, body_scale_x2, body_scale_x3
    integer, dimension(:), allocatable   :: i_start_obj_q1, i_end_obj_q1, j_start_obj_q1, j_end_obj_q1, k_start_obj_q1, k_end_obj_q1
    integer, dimension(:), allocatable   :: i_start_obj_q2, i_end_obj_q2, j_start_obj_q2, j_end_obj_q2, k_start_obj_q2, k_end_obj_q2
    integer, dimension(:), allocatable   :: i_start_obj_q3, i_end_obj_q3, j_start_obj_q3, j_end_obj_q3, k_start_obj_q3, k_end_obj_q3
    integer, dimension(:), allocatable   :: i_start_obj_sca, i_end_obj_sca, j_start_obj_sca, j_end_obj_sca, k_start_obj_sca, k_end_obj_sca
    ! Mesh info local ibm objects
    integer, dimension(:,:), allocatable    :: xstart_ibm_q1, xend_ibm_q1, xsize_ibm_q1  ! x-pencil
    integer, dimension(:,:), allocatable    :: ystart_ibm_q1, yend_ibm_q1, ysize_ibm_q1  ! y-pencil
    integer, dimension(:,:), allocatable    :: zstart_ibm_q1, zend_ibm_q1, zsize_ibm_q1  ! z-pencil
    logical, dimension(:), allocatable      :: object_in_current_proc_x_q1, object_in_current_proc_y_q1, object_in_current_proc_z_q1
    ! Mesh info local ibm objects
    integer, dimension(:,:), allocatable    :: xstart_ibm_q2, xend_ibm_q2, xsize_ibm_q2  ! x-pencil
    integer, dimension(:,:), allocatable    :: ystart_ibm_q2, yend_ibm_q2, ysize_ibm_q2  ! y-pencil
    integer, dimension(:,:), allocatable    :: zstart_ibm_q2, zend_ibm_q2, zsize_ibm_q2  ! z-pencil
    logical, dimension(:), allocatable      :: object_in_current_proc_x_q2, object_in_current_proc_y_q2, object_in_current_proc_z_q2
    ! Mesh info local ibm objects
    integer, dimension(:,:), allocatable    :: xstart_ibm_q3, xend_ibm_q3, xsize_ibm_q3  ! x-pencil
    integer, dimension(:,:), allocatable    :: ystart_ibm_q3, yend_ibm_q3, ysize_ibm_q3  ! y-pencil
    integer, dimension(:,:), allocatable    :: zstart_ibm_q3, zend_ibm_q3, zsize_ibm_q3  ! z-pencil
    logical, dimension(:), allocatable      :: object_in_current_proc_x_q3, object_in_current_proc_y_q3, object_in_current_proc_z_q3
    ! Mesh info local ibm objects
    integer, dimension(:,:), allocatable    :: xstart_ibm_sca, xend_ibm_sca, xsize_ibm_sca  ! x-pencil
    integer, dimension(:,:), allocatable    :: ystart_ibm_sca, yend_ibm_sca, ysize_ibm_sca  ! y-pencil
    integer, dimension(:,:), allocatable    :: zstart_ibm_sca, zend_ibm_sca, zsize_ibm_sca  ! z-pencil
    logical, dimension(:), allocatable      :: object_in_current_proc_x_sca, object_in_current_proc_y_sca, object_in_current_proc_z_sca

    !!!!!!! IN OBJECT
    ! Mesh info local ibm objects
    integer, dimension(:,:), allocatable    :: xstart_ibm_inobj_q1, xend_ibm_inobj_q1, xsize_ibm_inobj_q1  ! x-pencil
    integer, dimension(:,:), allocatable    :: ystart_ibm_inobj_q1, yend_ibm_inobj_q1, ysize_ibm_inobj_q1  ! y-pencil
    integer, dimension(:,:), allocatable    :: zstart_ibm_inobj_q1, zend_ibm_inobj_q1, zsize_ibm_inobj_q1  ! z-pencil
    logical, dimension(:), allocatable      :: object_in_current_proc_x_inobj_q1, object_in_current_proc_y_inobj_q1, object_in_current_proc_z_inobj_q1
    ! Mesh info local ibm objects
    integer, dimension(:,:), allocatable    :: xstart_ibm_inobj_q2, xend_ibm_inobj_q2, xsize_ibm_inobj_q2  ! x-pencil
    integer, dimension(:,:), allocatable    :: ystart_ibm_inobj_q2, yend_ibm_inobj_q2, ysize_ibm_inobj_q2  ! y-pencil
    integer, dimension(:,:), allocatable    :: zstart_ibm_inobj_q2, zend_ibm_inobj_q2, zsize_ibm_inobj_q2  ! z-pencil
    logical, dimension(:), allocatable      :: object_in_current_proc_x_inobj_q2, object_in_current_proc_y_inobj_q2, object_in_current_proc_z_inobj_q2
    ! Mesh info local ibm objects
    integer, dimension(:,:), allocatable    :: xstart_ibm_inobj_q3, xend_ibm_inobj_q3, xsize_ibm_inobj_q3  ! x-pencil
    integer, dimension(:,:), allocatable    :: ystart_ibm_inobj_q3, yend_ibm_inobj_q3, ysize_ibm_inobj_q3  ! y-pencil
    integer, dimension(:,:), allocatable    :: zstart_ibm_inobj_q3, zend_ibm_inobj_q3, zsize_ibm_inobj_q3  ! z-pencil
    logical, dimension(:), allocatable      :: object_in_current_proc_x_inobj_q3, object_in_current_proc_y_inobj_q3, object_in_current_proc_z_inobj_q3
    ! Mesh info local ibm objects
    integer, dimension(:,:), allocatable    :: xstart_ibm_inobj_sca, xend_ibm_inobj_sca, xsize_ibm_inobj_sca  ! x-pencil
    integer, dimension(:,:), allocatable    :: ystart_ibm_inobj_sca, yend_ibm_inobj_sca, ysize_ibm_inobj_sca  ! y-pencil
    integer, dimension(:,:), allocatable    :: zstart_ibm_inobj_sca, zend_ibm_inobj_sca, zsize_ibm_inobj_sca  ! z-pencil
    logical, dimension(:), allocatable      :: object_in_current_proc_x_inobj_sca, object_in_current_proc_y_inobj_sca, object_in_current_proc_z_inobj_sca
    real*8,dimension (3)                :: qbound
    integer                             :: number_of_objects
    integer                             :: interpol_type
    logical                             :: IBM_activated
    character(200)                      :: obj_file_path
    real*8                              :: initial_flowrate
    real*8                              :: ren_tau_ibm
    integer                             :: nb_flow_rate=200

end module IBM_settings

