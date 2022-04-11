

module time_schemes
    implicit none
    real*8, dimension(3)  ::ga,ro
    integer, parameter  :: EULER=0, ADAM_BASHFORTH=1, RK3=3
    integer             :: nb_substep

end module time_schemes

module schemes_interface
    use d1c_schemes
    use d2c_schemes
    use d0s_schemes
    use d1s_schemes
    use d1s_schemes_ibm
    implicit none

    procedure(D1c_ExpCtr_O2Fp0), pointer  ::  D1c, D2c
    procedure(D1c_ExpCtr_O2Fp0_ACC), pointer  ::  D1c_ACC, D2c_ACC
    procedure(D1c_ExpCtr_O2Fp0_MULT), pointer  ::  D1c_MULT, D2c_MULT
    procedure(D1c_ExpCtr_O2Fp0_MULT_ACC), pointer  ::  D1c_MULT_ACC, D2c_MULT_ACC

    procedure(D0s_ExpCtr_O2Fp0), pointer  ::  D0s, D1s
    procedure(D1s_ExpCtr_O2Fp0_ibm), pointer  ::  D1s_ibm

    procedure(D0s_ExpCtr_O2Fp0_MULT), pointer  :: D0s_MULTbyHimself
    procedure(D1s_ExpCtr_O2Fp0_MULT), pointer  ::  D1s_MULT
    procedure(D1s_ExpCtr_O2Fp0_ACC), pointer  ::  D1s_ACC
    procedure(D1s_ExpCtr_O2Fp0_MULT_ACC), pointer  :: D1s_MULT_ACC

end module schemes_interface


module schemes3D_interface
    use DRP_D1c
    use DRP_D2c
    use DRP_D1s
    use DRP_D1s_sh
    use DRP_D0s
    use DRP_D0s_sh

    use CPT_D1c
    use CPT_D2c
    use CPT_D1s
    use CPT_D1s_sh
    use CPT_D0s
    use CPT_D0s_sh

    use O2_D1c
    use O2_D2c
    use O2_D1s
    use O2_D1s_sh
    use O2_D0s
    use O2_D0s_sh
    implicit none

        procedure(D0s_DRP5_3Dx), pointer            :: D0s_3Dx
        procedure(D0ssh_DRP5_3Dx), pointer          :: D0ssh_3Dx
        procedure(D0s_DRP5_3Dy), pointer            :: D0s_3Dy
        procedure(D0ssh_DRP5_3Dy), pointer          :: D0ssh_3Dy
        procedure(D0ssh_DRP5_MULT_3Dy), pointer     :: D0ssh_MULT_3Dy
        procedure(D0s_DRP5_3Dz), pointer            :: D0s_3Dz
        procedure(D0ssh_DRP5_3Dz), pointer          :: D0ssh_3Dz
        procedure(D0ssh_DRP5_MULT_3Dz), pointer     :: D0ssh_MULT_3Dz

        procedure(D1s_DRP5_3Dx), pointer            :: D1s_3Dx
        procedure(D1s_DRP5_ACC_3Dx), pointer        :: D1s_ACC_3Dx
        procedure(D1ssh_DRP5_3Dx), pointer          :: D1ssh_3Dx
        procedure(D1ssh_DRP5_ACC_3Dx), pointer      :: D1ssh_ACC_3Dx
        procedure(D1s_DRP5_3Dy), pointer            :: D1s_3Dy
        procedure(D1ssh_DRP5_3Dy), pointer          :: D1ssh_3Dy
        procedure(D1s_DRP5_MULT_3Dy), pointer       :: D1s_MULT_3Dy
        procedure(D1s_DRP5_MULTACC_3Dy), pointer    :: D1s_MULTACC_3Dy
        procedure(D1ssh_DRP5_MULT_3Dy), pointer     :: D1ssh_MULT_3Dy
        procedure(D1ssh_DRP5_MULTACC_3Dy), pointer  :: D1ssh_MULTACC_3Dy
        procedure(D1s_DRP5_3Dz), pointer            :: D1s_3Dz
        procedure(D1ssh_DRP5_3Dz), pointer          :: D1ssh_3Dz
        procedure(D1s_DRP5_ACC_3Dz), pointer        :: D1s_ACC_3Dz
        procedure(D1ssh_DRP5_ACC_3Dz), pointer      :: D1ssh_ACC_3Dz

        procedure(D1c_DRP6_3Dx), pointer            :: D1c_3Dx
        procedure(D1c_DRP6_3Dy), pointer            :: D1c_3Dy
        procedure(D1c_DRP6_MULT_3Dy), pointer       :: D1c_MULT_3Dy
        procedure(D1c_DRP6_MULTACC_3Dy), pointer    :: D1c_MULTACC_3Dy
        procedure(D1c_DRP6_3Dz), pointer            :: D1c_3Dz

        procedure(D2c_DRP6_3Dx), pointer            :: D2c_3Dx
        procedure(D2c_DRP6_ACC_3Dx), pointer        :: D2c_ACC_3Dx
        procedure(D2c_DRP6_3Dy), pointer            :: D2c_3Dy
        procedure(D2c_DRP6_MULT_3Dy), pointer       :: D2c_MULT_3Dy
        procedure(D2c_DRP6_MULTACC_3Dy), pointer    :: D2c_MULTACC_3Dy
        procedure(D2c_DRP6_3Dz), pointer            :: D2c_3Dz


end module schemes3D_interface

module boundary_scheme
    implicit none
    real*8    :: a1_d,a2_d,a3_d,a1_u,a2_u,a3_u

contains

end module boundary_scheme
