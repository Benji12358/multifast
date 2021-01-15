
! ***********************************************************************************
! ***********************************************************************************
! ******************************** IO subroutines ***********************************
! ***********************************************************************************
! ***********************************************************************************


! -----------------------------------------------------------------------------------
!                          NAMESPACES FOR USE IN CORE MODULE
! -----------------------------------------------------------------------------------
module IMHD_IO_settings

    use MHD_settings, only: MHD_IO_read_settings => read_settings

end module IMHD_IO_settings

module IMHD_IO_dao
!
    use MHD_dao, only:       &
    MHD_IO_write_fields => write_fields
!
end module IMHD_IO_dao

module IMHD_IO_results_writer
!
    use MHD_results_writer, only:       &
    MHD_IO_write_fields => write_fields
!
end module IMHD_IO_results_writer

!module IMHD_IO_loader
!
!    use MHD_loader, only:   &
!    MHD_IO_load_fields => load_fields
!
!end module IMHD_IO_loader

module IMHD_IO
    use IMHD_IO_settings
    use IMHD_IO_results_writer
!    use IMHD_IO_loader
contains

end module IMHD_IO





module IMHD_LIFE

    use MHD_init, only: MHD_initialize => initialize

!    use MHD_recovery, only:    &
!    MHD_save_state             => save_state

end module IMHD_LIFE

!module IMHD_OPEN
!    use MHD_inout_flow_old, only: &
!    MHD_OPEN_init=>inoutflow_init,     &
!    MHD_OPEN_add_outflow=>add_outflow,     &
!    MHD_OPEN_get_inflow=>get_inflow, &
!    MHD_OPEN_newversion=>inout_newversion
!
!end module IMHD_OPEN
