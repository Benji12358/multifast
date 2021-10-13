1		MHD_state 0:MHD not resolved; 1: active MHD (RK3 ou AB2) ; 2: active MHD (Euler1) ; 3: passive MHD (Fext=0)
600.		Hartmann_number                         
0.		Magnetic_field_unit_vector_x 
0.		Magnetic_field_unit_vector_y 
0.		Magnetic_field_unit_vector_z
------------- Magnetic boundary layer ---------------------------------------------
0		layer_type ! 0 : uniform field, 1 : sinusoïdal, 2 : bipolar
0.3		delta_B, length (y-axis) of the layer of magnetic field ( /!\ applied on both walls /!\ )
0.1		smooth
42.		num_period ! Apply for sinusoïdal field => number of periods over channel length
-------------- Magnets ------------------------------------------------------------
0                   					Number of Magnet Pair
1 -1  -1  1  1 -1  -1  1  1 -1  -1  1  1 -1  -1  1     			Sigma                 ! -1:Magnetic Sink   1:Magnetic Source
40. 40. 40. 40.  40. 40. 40. 40. 40. 40. 40. 40.  40. 40. 40. 40.          	Magnet_size_x !! Direction 3 !! divided by pi
0.02 0.02 0.02 0.02  0.02 0.02 0.02 0.02   0.02 0.02 0.02 0.02  0.02 0.02 0.02 0.02	Magnet_size_z !! Direction 1 !! divided by pi
2. 2. 2. 2. 2. 2. 2. 2.  2. 2. 2. 2. 2. 2. 2. 2.          		Magnet_center_x !! Direction 3 !! divided by pi
0.083 0.25 0.083 0.25 0.417 0.583 0.417 0.583 0.75 0.917 0.75 0.917 1.083 1.25 1.083 1.25	Magnet_center_z !! Direction 1 !!divided by pi
-0.01 -0.01 2.01 2.01 -0.01 -0.01 2.01 2.01 -0.01 -0.01 2.01 2.01 -0.01 -0.01 2.01 2.01   		Magnet_center_y !! Direction 2 (<0 or >2)
-------------- Magnetic field from file -------------------------------------------
1					B_from_file 0: not from file, 1: from file
---------------- Export -----------------------------------------------------------
5000		nexport_MHD : mean current export frequency (meanJ.dat)
1		MHD_export_3D, exporting 3D fields (Elec_curr), 0:OFF 1:ON
