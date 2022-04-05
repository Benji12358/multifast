module schemes_settings
    implicit none
    logical :: treat_boundary_by_cpt_scheme=.true.
    integer, parameter  :: CPT_MIN=7

contains


    subroutine configure_explicit_schemes(treat_boundary_by_cpt_scheme_arg)
        implicit none
        logical :: treat_boundary_by_cpt_scheme_arg

        treat_boundary_by_cpt_scheme=treat_boundary_by_cpt_scheme_arg

    end subroutine configure_explicit_schemes

end module schemes_settings
