
! ***********************************************************************************
! ***********************************************************************************
! ******************************** IO subroutines ***********************************
! ***********************************************************************************
! ***********************************************************************************


! -----------------------------------------------------------------------------------
!                          NAMESPACES FOR USE IN CORE MODULE
! -----------------------------------------------------------------------------------








module Ifollowing_LIFE
    use following_initializer, only:    &
    following_INITIALIZE              => initialize

end module Ifollowing_LIFE

module Ifollowing_VELOCITY_IO_results_writer

    use following_velocity_results_writer, only:       &
    following_VELOCITY_IO_write_fields    => write_fields

end module Ifollowing_VELOCITY_IO_results_writer

module Ifollowing_SCALAR_IO_results_writer

    use following_scalar_results_writer, only:       &
    following_SCALAR_IO_write_fields    => write_fields

end module Ifollowing_SCALAR_IO_results_writer

module Ifollowing_VELOCITY_IO_loader

    use following_velocity_loader, only:   &
    following_VELOCITY_IO_load_fields => load_fields

end module Ifollowing_VELOCITY_IO_loader




module Ifollowing_VELOCITY_LIFE
    use following_velocity_init, only:    &
    following_velocity_initialize              => initialize

end module Ifollowing_VELOCITY_LIFE









module Ifollowing_IO
    use Ifollowing_VELOCITY_IO_results_writer
    use Ifollowing_VELOCITY_IO_loader

    use Ifollowing_SCALAR_IO_results_writer

end module Ifollowing_IO