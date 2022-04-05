include $(ENVIRONMENT_FILE)

FFLAGS		= -c $(COMPILER_OPTIONS) $(DECOMP2D_INC)
SOURCES		= sources
INCS 		= $(FLAG_FFT) $(HDF5_INC) $(IO_INC) $(SOLVERS_INC) $(MISCELLANEOUS_INC) $(SCHEMES_INC) -I $(COMMON)
DATA		=./data

OBJ=$(DATA)/scalar_data.o $(DATA)/workspace_view.o IO.o BC_controller.o initializer.o open.o solver.o Lifecycle.o interfaces.o 


./.lib/libscalar.a: $(OBJ)
	@echo;	\
	rm -f ./.lib/libscalar.a;	\
	$(AR_CMD) -q ./.lib/libscalar.a ./*.o $(DATA)/*.o ; 
	
	
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
	rm -f ./.lib/libscalar.a 