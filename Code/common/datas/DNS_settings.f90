
module DNS_settings
    implicit none

    ! FLOW_FROM_INFLOW: In turbulence generator, the flow will be build from the first 2D fields of inflow,
    ! obtained from an external simulation. This non-physical first solution garantees a minimal consistancy between
    ! inflow and inner field and avoid divergence of code
    integer, parameter  :: CONSTANT_FLOW=0, CHANNEL_FLOW=1, FLOW_FROM_INFLOW=2

    real*8         :: dt,ren, Uc, g, h_height
    integer         :: save_gradP_frequency ! mean_gradP export frequency
    integer         :: first_it=0
    integer         :: last_it
    integer         :: flow_type
    real*8          :: inflow_int=0.4d0
    integer         :: streamwise
    real*8          :: q1_x_av, q2_x_av, q3_x_av

end module DNS_settings

module inflow_settings
    implicit none

    integer, parameter   :: INFLOW_NONE=0, INFLOW_SQUARE=1

    integer         :: inflow_mode

    integer         :: inflow_sq_k1=22, inflow_sq_k2=42, inflow_sq_j1=40, inflow_sq_j2=70
    integer         :: outflow_buff, inflow_buff

contains

end module inflow_settings

module start_settings
    implicit none

    integer, parameter  :: NO_SOURCE=0, HDF5_FILE=1
    integer             :: start_it

    character*200       :: external_fields_path
    integer             :: start_source_type
    logical             :: start_from_coarse_file
    integer             :: n1c, n2c, n3c
    real*8              :: vper

end module start_settings

module blow_settings
    implicit none
    integer :: blow_end

end module blow_settings

module twave_settings
    implicit none
    integer :: twave_on
	integer :: inner_units
	real*8	:: Re_tau

end module twave_settings

module numerical_methods_settings
    implicit none
    integer         :: schemes_configuration
    integer         :: time_advancement
    logical         :: use_generic_poisson=.false.

    logical         :: inout_newversion=.true.

end module numerical_methods_settings
