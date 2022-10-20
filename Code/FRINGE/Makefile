include $(ENVIRONMENT_FILE)

FFLAGS		= -c $(COMPILER_OPTIONS) $(DECOMP2D_INC)
SOURCES		= sources
INCS 		= $(FLAG_FFT) $(HDF5_INC) $(IO_INC) $(SOLVERS_INC) $(MISCELLANEOUS_INC) $(SCHEMES_INC) -I $(COMMON) -I $(CORE)/passive_scalar
DATA		=./data

OBJ=$(DATA)/FRINGE_data.o $(DATA)/workspace_view.o IO.o initializer.o solver.o Lifecycle.o interfaces.o


./.lib/libfringe.a: $(OBJ)
	@echo;	\
	rm -f ./.lib/libfringe.a;	\
	$(AR_CMD) -q ./.lib/libfringe.a ./*.o $(DATA)/*.o ; 
	
	
$(DATA)/%.o:	$(DATA)/%.f90
	$(print_current_build);
	$(MUTE)$(FF) -o $@ $(FFLAGS) $(INCS) $< ;
	@if [ -z "$(MUTE)" ]; then echo; fi;
	
%.o:	%.f90
	$(print_current_build);
	$(MUTE)$(FF) $(FFLAGS) $(INCS) $<;
	@if [ -z "$(MUTE)" ]; then echo; fi;
	
clean :
	$(MUTE)rm -f ./*.o; rm -f ./*.mod;		\
	rm -f $(DATA)/*.o; rm -f $(DATA)/*.mod;		\
	rm -f ./.lib/libfringe.a 
