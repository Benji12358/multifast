
! ***********************************************************************************
! ***********************************************************************************
! ******************************** IO subroutines ***********************************
! ***********************************************************************************
! ***********************************************************************************


! -----------------------------------------------------------------------------------
!                          NAMESPACES FOR USE IN CORE MODULE
! -----------------------------------------------------------------------------------
module ISCALAR_IO_settings

    use SCALAR_settings, only:       &
    SCALAR_IO_read_settings => read_settings

end module ISCALAR_IO_settings

module ISCALAR_IO_dao

    use SCALAR_dao, only:       &
    SCALAR_IO_write_fields => write_fields

end module ISCALAR_IO_dao

module ISCALAR_IO_results_writer

    use SCALAR_results_writer, only:       &
    SCALAR_IO_write_fields => write_fields

end module ISCALAR_IO_results_writer

module ISCALAR_IO_loader

    use SCALAR_loader, only:   &
    SCALAR_IO_load_fields => load_fields

end module ISCALAR_IO_loader

module ISCALAR_IO
    use ISCALAR_IO_settings
    use ISCALAR_IO_results_writer
    use ISCALAR_IO_loader
contains

end module ISCALAR_IO





module ISCALAR_LIFE

    use SCALAR_init, only:    &
    SCALAR_initialize              => initialize

    use SCALAR_recovery, only:    &
    SCALAR_save_state             => save_state

end module ISCALAR_LIFE

module ISCALAR_OPEN
    use SCALAR_inout_flow_old, only: &
    SCALAR_OPEN_init=>inoutflow_init,     &
    SCALAR_OPEN_add_outflow=>add_outflow,     &
    SCALAR_OPEN_get_inflow=>get_inflow, &
    SCALAR_OPEN_newversion=>inout_newversion

end module ISCALAR_OPEN
