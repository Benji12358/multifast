
! ***********************************************************************************
! ***********************************************************************************
! ******************************** IO subroutines ***********************************
! ***********************************************************************************
! ***********************************************************************************


! -----------------------------------------------------------------------------------
!                          NAMESPACES FOR USE IN CORE MODULE
! -----------------------------------------------------------------------------------








module IEMBEDDED_LIFE
    use embedded_initializer, only:    &
    EMBEDDED_INITIALIZE              => initialize

end module IEMBEDDED_LIFE

module IEMBEDDED_VELOCITY_IO_results_writer

    use embedded_velocity_results_writer, only:       &
    EMBEDDED_VELOCITY_IO_write_fields    => write_fields

end module IEMBEDDED_VELOCITY_IO_results_writer

module IEMBEDDED_SCALAR_IO_results_writer

    use embedded_scalar_results_writer, only:       &
    EMBEDDED_SCALAR_IO_write_fields    => write_fields

end module IEMBEDDED_SCALAR_IO_results_writer

module IEMBEDDED_VELOCITY_IO_loader

    use embedded_velocity_loader, only:   &
    EMBEDDED_VELOCITY_IO_load_fields => load_fields

end module IEMBEDDED_VELOCITY_IO_loader




module IEMBEDDED_VELOCITY_LIFE
    use embedded_velocity_init, only:    &
    embedded_velocity_initialize              => initialize

end module IEMBEDDED_VELOCITY_LIFE









module IEMBEDDED_IO
    use IEMBEDDED_VELOCITY_IO_results_writer
    use IEMBEDDED_VELOCITY_IO_loader

    use IEMBEDDED_SCALAR_IO_results_writer

end module IEMBEDDED_IO