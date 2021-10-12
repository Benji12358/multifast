include $(ENVIRONMENT_FILE)

FFLAGS		= -c $(COMPILER_OPTIONS) $(DECOMP2D_INC)
SOURCES		= sources
INCS 		= $(FLAG_FFT) $(HDF5_INC) $(IO_INC) $(SOLVERS_INC) $(MISCELLANEOUS_INC) $(SCHEMES_INC) -I $(COMMON) -I $(CORE)/velocity

OBJ= mod_params.o mod_vars.o mod_myMPI.o mod_hdf5IO.o mod_inOut.o mod_restart.o mod_init.o mod_phase_average.o       

./.lib/libphaseaverage.a: $(OBJ)
	@echo;	\
	rm -f ./.lib/libphaseaverage.a;	\
	$(AR_CMD) -q ./.lib/libphaseaverage.a ./*.o ; 
	
	
%.o:	%.f90
	$(print_current_build);
	$(MUTE)$(FF) $(FFLAGS) $(INCS) $<;
	@if [ -z "$(MUTE)" ]; then echo; fi;
	
clean :
	$(MUTE)rm -f ./*.o; rm -f ./*.mod;		\
	rm -f ./.lib/libphaseaverage.a
