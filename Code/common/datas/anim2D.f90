module anim2D
    implicit none

    type anim_parameters
        integer                             :: nb_slices
        integer, dimension(10)              :: slices
        integer                             :: nb_steps
        integer                             :: step_size
        logical                             :: export_q1, export_q2, export_q3, export_pr
        logical                             :: export_om1, export_om2, export_om3
        logical                             :: export_fb1, export_fb2, export_fb3 !A ENLEVER
    end type anim_parameters

    type(anim_parameters)                   :: param_anim2D_1, param_anim2D_2, param_anim2D_3
    logical                                 :: anim2D_1, anim2D_2, anim2D_3

end module anim2D

module animBubble
    implicit none

    type slice
        integer                             :: p1, p2, p3    ! Central position (index in 1, 2 & 3)
        integer                             :: s1, s2, s3     ! "+-" dimension from central position into direction 1, 2, 3
        integer                             :: step_size
    end type slice

    integer                                 :: slice_cpt=0
    integer, parameter                      :: MAX_SLICES=100
    integer, parameter                      :: MAX_EXPORT=100
    type(slice), dimension(MAX_SLICES)      :: slice_tbl

end module animBubble
