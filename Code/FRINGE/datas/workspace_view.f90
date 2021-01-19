module FRINGE_workspace_view
    implicit none
!
    character(200)   :: results_path
    character(200)   :: results3D_path
!
!    character(200)   :: log_path
!    character(200)   :: sensor_file
!
!
!    character(200)   :: recovery_RHS_dir
!    character(200)   :: recovery_fields_dir
!    character(200)   :: recovery_sensor_file
!
!
end module FRINGE_workspace_view
!
module IFRINGE_workspace_view
!
    use FRINGE_workspace_view, only:                              &
!
    FRINGE_results_path             =>results_path,                 &
    FRINGE_results3D_path           =>results3D_path
!
!    FRINGE_sensor_file              =>sensor_file,                  &
!    FRINGE_recovery_sensor_file     =>recovery_sensor_file,         &
!    FRINGE_recovery_fields_dir      =>recovery_fields_dir,          &
!    FRINGE_recovery_RHS_dir         =>recovery_RHS_dir
!
end module IFRINGE_workspace_view
