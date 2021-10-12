
module IBM_settings
    implicit none
    integer, parameter      :: IBM_INTERPOL_NONE=0, IBM_FADLUN=1, ANTISYMMETRIC_INTERPOL=2, SECOND_ORDER_INTERPOL=3
    real*8                  :: body_x1, body_x2, body_x3
    real*8                  :: body_scale_x1, body_scale_x2, body_scale_x3
    real*8,dimension (3)    :: qbound
    integer                 :: interpol_type
    logical                 :: IBM_activated
    character(200)          :: obj_file_path
    integer                 :: nb_flow_rate=200

end module IBM_settings

