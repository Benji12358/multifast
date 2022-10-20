module embedded_mesh
    implicit none
    integer     :: n1, n2, n3
    integer     :: n1m, n2m, n3m
    integer     :: n1_fringe, n3_fringe

    real*8  :: L1, L2, L3
    real*8  :: stretch_Y

    integer, parameter  :: ORLANDI_MESH=1, LAMBALLAIS_MESH=2, NO_TRANSFO_MESH=3
    integer :: mesh_type

    real*8  :: dx1, dx2, dx3

    real*8, dimension(:), allocatable :: X, Xc
    real*8, dimension(:), allocatable :: Y, Yc
    real*8, dimension(:), allocatable :: Z, Zc
    real*8, dimension(:), allocatable :: cell_size_Y
    real*8, dimension(:,:,:), allocatable :: Y_field

end module embedded_mesh



module embedded_irregular_derivative_coefficients

    implicit none
    real*8,dimension(:), allocatable    :: Y_to_YTr_for_D1
    real*8,dimension(:), allocatable    :: Yc_to_YcTr_for_D1

    real*8,dimension(:, :), allocatable :: Y_to_YTr_for_D2
    real*8,dimension(:, :), allocatable :: Yc_to_YcTr_for_D2

end module embedded_irregular_derivative_coefficients
