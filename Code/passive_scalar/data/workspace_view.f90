module SCALAR_workspace_view
    implicit none

    character(200)   :: results_path
    character(200)   :: results3D_path

    character(200)   :: log_path
    character(200)   :: sensor_file


    character(200)   :: recovery_RHS_dir
    character(200)   :: recovery_fields_dir
    character(200)   :: recovery_sensor_file


end module SCALAR_workspace_view

module ISCALAR_workspace_view

    use SCALAR_workspace_view, only:                              &

    SCALAR_results_path             =>results_path,                 &
    SCALAR_results3D_path           =>results3D_path,               &

    SCALAR_sensor_file              =>sensor_file,                  &
    SCALAR_recovery_sensor_file     =>recovery_sensor_file,         &
    SCALAR_recovery_fields_dir      =>recovery_fields_dir,          &
    SCALAR_recovery_RHS_dir         =>recovery_RHS_dir

end module ISCALAR_workspace_view
