module following_workspace_view
    implicit none

    character(200)   :: results3D_path

    character(200)   :: lapY_file

    character(200)   :: recovery_path
    character(200)   :: snapshot_path

end module following_workspace_view

module following_common_workspace_view

    use following_workspace_view, only:                           &
    following_common_results3D_path      =>results3D_path,   &
    following_common_lapY_file           =>lapY_file,        &
    following_common_recovery_path       =>recovery_path,    &
    following_common_snapshot_path       =>snapshot_path

end module following_common_workspace_view
