1000  100  11000    		nb_iteration   save3D_frequency  checkpoint_frequency
5	2	0			Boundary conditions
1					streamwise direction
1			flow_type (CONSTANT_FLOW=0, CHANNEL_FLOW=1, FLOW_FROM_INFLOW=2, BOUNDARY_LAYER_FLOW=4)
0.25			delta_BL
0			outflow_type 														# Only for inflow/outflow
0			outflow_buff (0: not outflow export, >0 size of outflow buffer)		# Only for inflow/outflow
0			inflow_buff (0: not inflow import, >0 size of inflow buffer)		# Only for inflow/outflow
0			inflow_type (INFLOW_NONE=0, INFLOW_SQUARE=1)						# Only for inflow/outflow
/home/bureau/WORKSPACE_MAURICIO/MachineVirtuelleTrial/DNS_MAURICIO/Simulations/OPEN_VERSION0/Input/Outflow/
34		nx_start, i.e. which x is used as inflow
1.225		divro, Liquid_Density/Gas_Density
0.1			d, diameter_of_the_bubble (mm) - Mean Diameter
-9810.0		g, gravity (mm/(s*s))
1000.0		Reynolds_number (U_poiseuille*H/kinematic_viscosity)
1.0			kinematic_viscosity (mm*mm/s) ! characteristic viscosity
1.0			h_height, half_channel_height (mm) ! characteristic length
1000.0		Uc, velocity_center_Poiseuille (mm/s) ! characteristic velocity
500			save_gradP_frequency, mean streamwise pressure gradient export frequency (gradP.dat)
