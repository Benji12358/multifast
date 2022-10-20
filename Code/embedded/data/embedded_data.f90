module embedded_data
    use decomp_2d, only:DECOMP_INFO 
    implicit none

    TYPE(DECOMP_INFO) :: decomp_embedded

    ! Mesh info embedded mesh
    integer, dimension(3)                   :: xstart_e, xend_e, xsize_e  ! x-pencil
    integer, dimension(3)                   :: ystart_e, yend_e, ysize_e  ! y-pencil
    integer, dimension(3)                   :: zstart_e, zend_e, zsize_e  ! z-pencil

    real*8                                  :: u_bulk

    logical                                 :: use_embedded

end module embedded_data

module embedded_physical_fields
    implicit none

    real*8,dimension (:,:,:), allocatable       :: q3_x, q3_y, q3_z
    real*8,dimension (:,:,:), allocatable       :: q2_x, q2_y, q2_z
    real*8,dimension (:,:,:), allocatable       :: q1_x, q1_y, q1_z
    real*8,dimension (:,:,:), allocatable       :: dp_x, dp_y, dp_z

    real*8      :: dP_streamwise, dP_spanwise

    real*8,dimension (:,:), allocatable         :: q3_wall10, q2_wall10, q1_wall10
    real*8,dimension (:,:), allocatable         :: q3_wall11, q2_wall11, q1_wall11

    real*8,dimension (:,:), allocatable         :: q3_wall20, q2_wall20, q1_wall20
    real*8,dimension (:,:), allocatable         :: q3_wall21, q2_wall21, q1_wall21

    real*8,dimension (:,:), allocatable         :: q3_wall30, q2_wall30, q1_wall30
    real*8,dimension (:,:), allocatable         :: q3_wall31, q2_wall31, q1_wall31

    real*8,dimension (:,:,:), allocatable       :: divu_x, divu_y, divu_z, source_term
    real*8,dimension (:,:,:), allocatable       :: dphidx1_x, dphidx2_x, dphidx3_x
    real*8,dimension (:,:,:), allocatable       :: dphidx1_y, dphidx2_y, dphidx3_y
    real*8,dimension (:,:,:), allocatable       :: dphidx1_z, dphidx2_z, dphidx3_z

end module embedded_physical_fields

module embedded_start_settings
    implicit none

    integer, parameter  :: NO_SOURCE=0, HDF5_FILE=1
    integer             :: start_it
    integer             :: index_for_output

    character*200       :: external_fields_path
    integer             :: start_source_type
    real*8              :: vper

end module embedded_start_settings

module embedded_fringe_data
    implicit none

    integer, parameter   :: POISEUILLE_INFLOW=0, SQUARE_INFLOW=1, INFLOW_FROM_FILE=2, INFLOW_FROM_PREV_RUN=3, BOUNDARY_LAYER_INFLOW=4

    logical                                   :: use_fringe=.false.
    real*8                                    :: max_strength_damping, delta_rise, delta_fall, fringe_length
    integer                                   :: n_fringe_start, n_fringe_end, n_delta_fringe, n_fringe_region, n_interest_region
    integer                                   :: inflow_buff, inflow_type
    integer                                   :: n_delta_rise, n_delta_fall, delta_activation

    real*8, dimension(:,:,:), allocatable     :: f1_fringe_x, f2_fringe_x, f3_fringe_x
    real*8, dimension(:,:), allocatable       :: sca_inflow
    real*8, dimension(:), allocatable         :: lambda_x

end module embedded_fringe_data

module embedded_scalar_data
    implicit none

    integer                                     :: SCA_state
    real*8,dimension (:,:,:), allocatable       :: sca_x, sca_y, sca_z
    real*8                                      :: prandtl, delta_T, heat_flux

    real*8,dimension (:,:), allocatable         :: sca_wall10, sca_wall20, sca_wall30
    real*8,dimension (:,:), allocatable         :: sca_wall11, sca_wall21, sca_wall31

end module embedded_scalar_data

