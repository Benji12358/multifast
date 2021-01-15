

!module MHD_bc_controller
!    implicit none
!
!contains
!
!    subroutine set_numerical_boundaries()
!
!        use boundaries
!        implicit none
!
!        select case (MHD_BC1)
!
!            case (FIXED_VALUE, OPEN)
!
!            TRANSPORT_Q1S_BC1=Dirichlet
!            TRANSPORT_MHD_BC1=Dirichlet
!
!            case (UNBOUNDED)
!
!            TRANSPORT_Q1S_BC1=periodic
!            TRANSPORT_MHD_BC1=periodic
!
!        end select
!
!        select case (SCA_BC2)
!            ! TODO : Add open condition
!            case (FIXED_VALUE)
!
!            TRANSPORT_Q2S_BC2=Dirichlet
!            TRANSPORT_MHD_BC2=Dirichlet
!
!            case (FIXED_FLUX)
!                ! TODO
!
!            case (UNBOUNDED)
!
!            TRANSPORT_Q2S_BC2=periodic
!            TRANSPORT_MHD_BC2=periodic
!
!        end select
!
!        select case (MHD_BC3)
!
!            ! TODO : Add open condition
!            case (FIXED_VALUE)
!
!            TRANSPORT_Q3S_BC3=Dirichlet
!            TRANSPORT_MHD_BC3=Dirichlet
!
!            case (FIXED_FLUX)
!                ! TODO
!
!            case (UNBOUNDED)
!
!            TRANSPORT_Q3S_BC3=periodic
!            TRANSPORT_MHD_BC3=periodic
!
!        end select
!
!        !write(*,*)'TRANSPORT_Q2S_BC2', TRANSPORT_Q2S_BC2, Dirichlet
!
!    end subroutine set_numerical_boundaries
!
!end module MHD_bc_controller
