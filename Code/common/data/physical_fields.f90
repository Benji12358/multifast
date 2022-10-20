module physical_fields
    implicit none

    real*8,dimension (:,:,:), allocatable       :: q3_x, q3_y, q3_z
    real*8,dimension (:,:,:), allocatable       :: q2_x, q2_y, q2_z
    real*8,dimension (:,:,:), allocatable       :: q1_x, q1_y, q1_z
    real*8,dimension (:,:,:), allocatable       :: dp_x, dp_y, dp_z
    !real*8,dimension (:,:,:), allocatable       :: pr_x, pr_y, pr_z

    real*8      :: dP_streamwise, dP_spanwise

    real*8,dimension (:,:), allocatable         :: q3_wall10, q2_wall10, q1_wall10
    real*8,dimension (:,:), allocatable         :: q3_wall11, q2_wall11, q1_wall11

    real*8,dimension (:,:), allocatable         :: q3_wall20, q2_wall20, q1_wall20
    real*8,dimension (:,:), allocatable         :: q3_wall21, q2_wall21, q1_wall21

    real*8,dimension (:,:), allocatable         :: q3_wall30, q2_wall30, q1_wall30
    real*8,dimension (:,:), allocatable         :: q3_wall31, q2_wall31, q1_wall31

    real*8,dimension (:,:,:), allocatable       :: divu_x, divu_y, divu_z, source_term
    real*8,dimension (:,:,:), allocatable       :: dphidx1_x, dphidx2_x, dphidx3_x
    real*8,dimension (:,:,:), allocatable       :: dphidx1_y, dphidx2_y, dphidx3_y
    real*8,dimension (:,:,:), allocatable       :: dphidx1_z, dphidx2_z, dphidx3_z

    real*8,dimension (:,:,:), allocatable       :: om1_x, om2_x, om3_x
    real*8,dimension (:,:,:), allocatable       :: void_x
    real*8,dimension (:,:,:), allocatable       :: gradtau3D_x


end module physical_fields


