module scalar_datas
    implicit none

    real*8,dimension (:,:,:), allocatable       :: sca_x, sca_y, sca_z
    real*8  :: renprandtl, delta_T, heat_flux
    integer :: sca_activated

    real*8,dimension (:,:), allocatable         :: sca_wall10, sca_wall20, sca_wall30
    real*8,dimension (:,:), allocatable         :: sca_wall11, sca_wall21, sca_wall31

    real*8,dimension (:,:,:,:), allocatable       :: sca_outflow
    real*8,dimension (:,:,:,:), allocatable       :: sca_inflow

    real*8, dimension(:,:,:), allocatable       :: fb1_scalar_x, fb2_scalar_x, fb3_scalar_x
    real*8, dimension(:,:,:), allocatable       :: fb1_scalar_y, fb2_scalar_y, fb3_scalar_y
    real*8, dimension(:,:,:), allocatable       :: fb1_scalar_z, fb2_scalar_z, fb3_scalar_z

    real*8                                      :: vol_exp_coef = 0.000207d0

end module scalar_datas
