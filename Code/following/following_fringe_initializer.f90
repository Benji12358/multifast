module following_fringe_initializer
    use decomp_2d
    use mpi
    use following_mesh
    use following_fringe_data
    use DNS_settings
    use following_physical_fields
    use following_data

    use following_scalar_data

!
    implicit none
!
contains
!
!
    subroutine get_inflow()

        implicit none
        integer :: j, k, mpi_err
        logical :: pair_n2
        real*8 :: delta, ut1
        real*8, dimension(n2) :: tmp_profile

        ! Assume streamwise direction is 1
        ! so streamwise velocity is q1

        sca_inflow = 0.d0 
        q1_inflow = 0.d0 
        q2_inflow = 0.d0 
        q3_inflow = 0.d0

        return

    end subroutine get_inflow

end module following_fringe_initializer
