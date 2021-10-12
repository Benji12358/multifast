module VELOCITY_workspace_view
    implicit none

    character(200)   :: lapY_file

    character(200)   :: results_path
    character(200)   :: results3D_path

    character(200)   :: log_path
    character(200)   :: divergence_history_file
    character(200)   :: kinetic_history_file
    character(200)   :: sensor_file


    character(200)   :: recovery_RHS_dir
    character(200)   :: recovery_fields_dir
    character(200)   :: recovery_run_ctxt_file
    character(200)   :: recovery_sensor_file
    character(200)   :: recovery_divergence_history_file
    character(200)   :: recovery_kinetic_history_file


end module VELOCITY_workspace_view

module IVELOCITY_workspace_view

    use VELOCITY_workspace_view, only:                              &

    VELOCITY_results_path               =>results_path,             &
    VELOCITY_results3D_path             =>results3D_path,           &

    VELOCITY_sensor_file              =>sensor_file,                  &
    VELOCITY_lapY_file                =>lapY_file,                    &
    VELOCITY_divergence_history_file  =>divergence_history_file,      &
    VELOCITY_kinetic_history_file     =>kinetic_history_file,         &
    VELOCITY_recovery_sensor_file                 =>recovery_sensor_file,                     &
    VELOCITY_recovery_divergence_history_file     =>recovery_divergence_history_file,         &
    VELOCITY_recovery_kinetic_history_file        =>recovery_kinetic_history_file,            &
    VELOCITY_recovery_fields_dir                =>recovery_fields_dir,                    &
    VELOCITY_recovery_RHS_dir                     =>recovery_RHS_dir

end module IVELOCITY_workspace_view
