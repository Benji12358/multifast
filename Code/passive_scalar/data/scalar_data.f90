module SCALAR_data
    implicit none

    integer, parameter  :: CLASSIC_INIT=0, INIT_FROM_FILE=1, KAWAMURA_INIT=2, CONSTANT_HEAT_FLUX=3

    real*8,dimension (:,:,:), allocatable       :: sca_x, sca_y, sca_z
    real*8  :: renprandtl, delta_T, heat_flux
    integer :: SCA_state
    integer :: init_type
    logical :: reset_scalar_field

    real*8,dimension (:,:), allocatable         :: sca_wall10, sca_wall20, sca_wall30
    real*8,dimension (:,:), allocatable         :: sca_wall11, sca_wall21, sca_wall31

    real*8,dimension (:,:,:,:), allocatable       :: sca_outflow
    real*8,dimension (:,:,:,:), allocatable       :: sca_inflow

    real*8, dimension(:,:,:), allocatable       :: fb1_scalar_x, fb2_scalar_x, fb3_scalar_x
    real*8, dimension(:,:,:), allocatable       :: fb1_scalar_y, fb2_scalar_y, fb3_scalar_y
    real*8, dimension(:,:,:), allocatable       :: fb1_scalar_z, fb2_scalar_z, fb3_scalar_z

    real*8                                      :: vol_exp_coef = 0.000207d0

end module SCALAR_data
