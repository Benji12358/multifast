module mod_params

	use COMMON_workspace_view   ! To get the path of the input file
 	use VELOCITY_workspace_view ! To get common results path
 	use mesh                    ! To get mesh parameters
	use decomp_2d               ! To get domain decomposition
 	use run_ctxt_data           ! To get time 't' and 'ntime'
 	use mathematical_constants  ! To get the value of pi
 	use DNS_settings            ! To get the streamwise direction
 	use boundaries              ! To get Travelling wave parameters
 	use twave_settings          ! To get request for inner to outer units conversion
 	use hdf5                    ! For collective IO
 	use start_settings          ! To get start settings
 	use velocity_operations     ! To get subroutines for interpolating velocities at
                              	! cell center and calculating vorticity
	
	integer :: d_freq,restart

	integer :: phase_averaging
	integer :: xbins
	integer :: xbins_old

	real*8  :: amp,kappa,omega,wave_speed
	real*8  :: amp_old,kappa_old,omega_old
	real*8  :: cf_amp,cf_kappa,cf_omega
	real*8  :: lambda,xmax,xmin,delX

	integer :: recover
	
end module mod_params