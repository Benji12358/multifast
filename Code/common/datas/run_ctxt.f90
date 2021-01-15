module run_ctxt_data
    implicit none

    integer, parameter  :: FINISHED=1, IN_PROGRESS=2
    integer, parameter  :: NEW_SIMULATION=0, CONTINUE_FROM_PREVIOUS_RUN=FINISHED, RECOVERY_A_RUN=IN_PROGRESS
    logical             :: test_mode
    integer             :: run_ctxt
    integer             :: run_status
    integer             :: ntime
    real*8              :: t=0.d0, tinit=0.d0

end module run_ctxt_data
