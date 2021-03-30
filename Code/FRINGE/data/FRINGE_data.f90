module fringe_data
    implicit none

    integer, parameter   :: POISEUILLE_INFLOW=0, SQUARE_INFLOW=1, INFLOW_FROM_FILE=2, BOUNDARY_LAYER_INFLOW=4

    logical                                   :: use_fringe=.false.
    real*8                                    :: max_strength_damping, delta_rise, delta_fall, fringe_length
    integer                                   :: n_fringe_start, n_fringe_end, n_delta_fringe, n_fringe_region
    integer                                   :: inflow_buff, inflow_type
    integer                                   :: n_delta_rise, n_delta_fall, delta_activation

    real*8, dimension(:,:,:), allocatable     :: f1_fringe_x, f2_fringe_x, f3_fringe_x
    real*8, dimension(:,:,:), allocatable     :: f1_fringe_z, f2_fringe_z, f3_fringe_z
    real*8, dimension(:,:), allocatable       :: q1_inflow, q2_inflow, q3_inflow
    real*8, dimension(:), allocatable         :: lambda_x, lambda_z

end module fringe_data
