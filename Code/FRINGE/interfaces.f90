
! ***********************************************************************************
! ***********************************************************************************
! ******************************** IO subroutines ***********************************
! ***********************************************************************************
! ***********************************************************************************


! -----------------------------------------------------------------------------------
!                          NAMESPACES FOR USE IN CORE MODULE
! -----------------------------------------------------------------------------------
module IFRINGE_IO_settings

    use FRINGE_settings, only: FRINGE_IO_read_settings => read_settings

end module IFRINGE_IO_settings

! module IFRINGE_IO_dao
! !
!     use FRINGE_dao, only:       &
!     FRINGE_IO_write_fields => write_fields
! !
! end module IFRINGE_IO_dao

! module IFRINGE_IO_results_writer
! !
!     use FRINGE_results_writer, only:       &
!     FRINGE_IO_write_fields => write_fields
! !
! end module IFRINGE_IO_results_writer

!module IFRINGE_IO_loader
!
!    use FRINGE_loader, only:   &
!    FRINGE_IO_load_fields => load_fields
!
!end module IFRINGE_IO_loader

! module IFRINGE_IO
!     use IFRINGE_IO_settings
!     use IFRINGE_IO_results_writer
! !    use IFRINGE_IO_loader
! contains

! end module IFRINGE_IO





module IFRINGE_LIFE

    use FRINGE_init, only: FRINGE_initialize => initialize

!    use FRINGE_recovery, only:    &
!    FRINGE_save_state             => save_state

end module IFRINGE_LIFE

!module IFRINGE_OPEN
!    use FRINGE_inout_flow_old, only: &
!    FRINGE_OPEN_init=>inoutflow_init,     &
!    FRINGE_OPEN_add_outflow=>add_outflow,     &
!    FRINGE_OPEN_get_inflow=>get_inflow, &
!    FRINGE_OPEN_newversion=>inout_newversion
!
!end module IFRINGE_OPEN
