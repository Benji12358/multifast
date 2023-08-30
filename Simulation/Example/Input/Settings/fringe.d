1		fringe_state 0: no fringe ; 1: fringe activated (if fringe is activated, activate at least one BC to 5 in global.d)
0.3	fringe_length (in percentage of L_streamwise
0.5		delta_rise (in percentage of L_fringe)
0.15		delta_fall (in percentage of L_fringe)
0		delta_activation
0.05		max_strength_damping
------------- Inflow ---------------------------------------------
1		inflow_type 0 : poiseuille inflow, 1 : square inflow, 2 : inflow from file, 3: inlfow from previous run, 4 : boundary layer inflow
---------------- Export -----------------------------------------------------------
500		nexport_MHD : mean current export frequency (meanJ.dat)
1		MHD_export_3D, exporting 3D fields (Elec_curr), 0:OFF 1:ON
