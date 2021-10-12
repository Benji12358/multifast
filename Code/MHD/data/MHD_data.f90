module MHD_data
    implicit none

    integer                                   :: MHD_state ! 0: do not solve MHD ; 1: active MHD (RK3 ou AB2) ; 2: active MHD (Euler1) ; 3: passive MHD (Fext=0)
    integer                                   :: Monodirec_Magnetic_field, magnet_number_pairs
    real*8                                    :: Magnetic_field_unit_vector_x, Magnetic_field_unit_vector_y , Magnetic_field_unit_vector_z
    real*8,dimension (:,:,:), allocatable     :: phi_x, phi_y, phi_z
    real*8                                    :: Stuart_number,Hartmann_number, delta_B,smooth, num_period
    integer                                   :: nexport_MHD, layer_type, MHD_export_3D

    real*8, dimension(:,:,:), allocatable     :: fb1_MHD_x, fb2_MHD_x, fb3_MHD_x
    real*8,dimension (:,:,:), allocatable     :: gradphi1_x, gradphi2_y, gradphi3_z
    real*8,dimension (:,:,:), allocatable     :: ucrossB1_x, ucrossB2_y, ucrossB3_z

    real*8,dimension (:,:,:), allocatable     :: B01_x, B02_y, B03_z
    real*8,dimension (:,:,:), allocatable     :: B01c_y, B02c_y, B03c_y
    real*8,dimension (:,:,:), allocatable     :: A1_x, A1_y, A2_y, A3_y, A3_z
    real*8,dimension (:,:,:), allocatable     :: RHS_z,RHS_z2

    type magnet_array
        integer  :: magnet_wall_location
        integer  :: sigma
        real*8  :: magnet_size_x
        real*8  :: magnet_size_z
        real*8  :: magnet_center_x
        real*8  :: magnet_center_z
        real*8  :: magnet_center_y
    end type magnet_array

    integer, parameter                              :: MAX_MAGNET=50
    type(magnet_array), dimension(MAX_MAGNET)       :: mag_array

end module MHD_data
