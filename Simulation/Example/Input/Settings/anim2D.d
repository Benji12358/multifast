0 0 0			! export X1, X2 and X3 slices. 0: No, 1: export in single hdf5 file, 2: export in multiple hdf5 files
X1 direction
1  	1		! slices
1000	200			! nb_steps step
1 	1	1	1	! q1 q2 q3 pr
1				! sca

X2 direction
3  	10 40 50 60		! slices
2000	200			! nb_steps step
0 	1	0	0	! q1 q2 q3 pr

X3 direction
1  	1			! slices
1000	200			! nb_steps step
1 	1	1	1	! q1 q2 q3 pr
1				! sca

Bubbles
0		# Number of exported slice/box. Each slice has a next dedicated column: slice1 --> 1st column, slice2 --> 2nd column, ...
1		# Bubble coordinate - Direction 1 : being set from 1 to 129 (use the current computational domain)	
1		# Bubble coordinate - Direction 2 : being set from 1 to 109 (use the current computational domain)
1		# Bubble coordinate - Direction 3 : being set from 1 to 65 (use the current computational domain)
129		# Span around the bubble coordinate - Direction 1 : being set from 1 to 129 (use the current computational domain)
109		# Span around the bubble coordinate - Direction 2 : being set from 1 to 109 (use the current computational domain)
65		# Span around the bubble coordinate - Direction 3 : being set from 1 to 65 (use the current computational domain)
0		# Save/Export frequency --> 1 = every dt, 2 = every 2 dt, etc.
