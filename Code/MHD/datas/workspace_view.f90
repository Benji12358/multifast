module MHD_workspace_view
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
end module MHD_workspace_view
!
module IMHD_workspace_view
!
    use MHD_workspace_view, only:                              &
!
    MHD_results_path             =>results_path,                 &
    MHD_results3D_path           =>results3D_path
!
!    MHD_sensor_file              =>sensor_file,                  &
!    MHD_recovery_sensor_file     =>recovery_sensor_file,         &
!    MHD_recovery_fields_dir      =>recovery_fields_dir,          &
!    MHD_recovery_RHS_dir         =>recovery_RHS_dir
!
end module IMHD_workspace_view
