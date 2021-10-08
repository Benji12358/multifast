
! ***********************************************************************************
! ***********************************************************************************
! ******************************** IO subroutines ***********************************
! ***********************************************************************************
! ***********************************************************************************


! -----------------------------------------------------------------------------------
!                          NAMESPACES FOR USE IN CORE MODULE
! -----------------------------------------------------------------------------------
module IVELOCITY_IO_settings

    use VELOCITY_settings, only:       &
    VELOCITY_IO_read_blowing_settings => read_blowing_settings

end module IVELOCITY_IO_settings

module IVELOCITY_IO_dao

    use VELOCITY_dao, only:       &
    VELOCITY_IO_write_fields => write_fields

end module IVELOCITY_IO_dao

module IVELOCITY_IO_results_writer

    use VELOCITY_results_writer, only:       &
    VELOCITY_IO_write_fields    => write_fields

    use VELOCITY_dao, only:       &
    CELL_CENTERED_IO_WRITE_FIELDS    => hdf5output

end module IVELOCITY_IO_results_writer

module IVELOCITY_IO_loader

    use VELOCITY_loader, only:   &
    VELOCITY_IO_load_fields => load_fields

end module IVELOCITY_IO_loader

module IVELOCITY_IO_log_writer

    use VELOCITY_log_writer, only:  &
    VELOCITY_IO_write_divergence   => export_divergence_file,      &
    VELOCITY_IO_write_kinetic      => export_kinetic,              &
    VELOCITY_IO_write_kinetic_IBM  =>export_kinetic_IBM,           &
    VELOCITY_IO_write_velmax_IBM   =>export_velmax_IBM,           &
    VELOCITY_IO_write_properties   =>export_properties,         &
    VELOCITY_IO_write_gradP   =>export_gradP


end module IVELOCITY_IO_log_writer

module IVELOCITY_IO
    use IVELOCITY_IO_results_writer
    use IVELOCITY_IO_loader
    use IVELOCITY_IO_log_writer

end module IVELOCITY_IO



module IVELOCITY_LIFE
    use VELOCITY_init, only:    &
    VELOCITY_initialize              => initialize

    use VELOCITY_recovery, only:    &
    VELOCITY_save_state             => save_state

end module IVELOCITY_LIFE


module IVELOCITY_OPEN
    use VELOCITY_inout_flow_old, only: &
    VELOCITY_OPEN_init=>inoutflow_init,     &
    VELOCITY_OPEN_add_outflow=>add_outflow,     &
    VELOCITY_OPEN_get_inflow=>get_inflow

end module IVELOCITY_OPEN
