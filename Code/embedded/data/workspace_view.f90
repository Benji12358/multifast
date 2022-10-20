module embedded_workspace_view
    implicit none

    character(200)   :: results3D_path

    character(200)   :: lapY_file

    character(200)   :: recovery_path
    character(200)   :: snapshot_path

end module embedded_workspace_view

module embedded_common_workspace_view

    use workspace_view, only:                           &
    embedded_common_results3D_path      =>results3D_path,   &
    embedded_common_lapY_file           =>lapY_file,        &
    embedded_common_recovery_path       =>recovery_path,    &
    embedded_common_snapshot_path       =>snapshot_path

end module embedded_common_workspace_view
