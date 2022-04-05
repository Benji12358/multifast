module workspace_view
    implicit none

    character(50)   :: simulation_name
    character(200)   :: simulation_path

    character(200)   :: settings_path
    character(200)   :: lapY_file

    character(200)   :: log_path
    character(200)   :: divergence_history_file
    character(200)   :: kinetic_history_file

    character(200)   :: results_path
    character(200)   :: results3D_path


    character(200)   :: recovery_path
    character(200)   :: recovery_outflow_dir
    character(200)   :: recovery_RHS_dir
    character(200)   :: recovery_run_ctxt_file
    character(200)   :: recovery_sensor_file
    character(200)   :: recovery_divergence_history_file
    character(200)   :: recovery_kinetic_history_file


    character(200)   :: snapshot_path
    character(200)   :: anim2D_path
    character(200)   :: anim2D_bubble_path


    character(200)   :: outflow_path
    character(200)   :: inflow_path

end module workspace_view

module COMMON_workspace_view
    use workspace_view, only:                                       &
    COMMON_simulation_name          =>simulation_name,              &
    COMMON_simulation_path          =>simulation_path,              &
    COMMON_settings_path            =>settings_path,                &
    COMMON_lapY_file                =>lapY_file,                    &
    COMMON_log_path                 =>log_path,                     &
    COMMON_results_path             =>results_path,                 &
    COMMON_recovery_path                        =>recovery_path,                            &
    COMMON_recovery_outflow_dir                 =>recovery_outflow_dir,                            &
    COMMON_recovery_run_ctxt_file               =>recovery_run_ctxt_file,                   &
    COMMON_snapshot_path                        =>snapshot_path,                            &
    COMMON_anim2D_path                          =>anim2D_path,                              &
    COMMON_anim2D_bubble_path                   =>anim2D_bubble_path,                       &
    COMMON_outflow_path                         =>outflow_path,                             &
    COMMON_inflow_path                          =>inflow_path,                              &
    COMMON_results3D_path                       =>results3D_path

end module COMMON_workspace_view
