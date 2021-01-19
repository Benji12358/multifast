include $(ENVIRONMENT_FILE)

FFLAGS		= -c $(COMPILER_OPTIONS) $(DECOMP2D_INC)
SOURCES		= sources
INCS 		= $(FLAG_FFT) $(HDF5_INC) $(IO_INC) $(SOLVERS_INC) $(MISCELLANEOUS_INC) $(SCHEMES_INC) -I $(COMMON) -I $(CORE) -I $(CORE)/FRINGE
DATA		=./data
FIELDS		=./fields_utils
OBJ			=$(DATA)/workspace_view.o $(DATA)/settings.o $(FIELDS)/operations.o $(FIELDS)/properties.o $(FIELDS)/quality_criteriums.o initializer.o BC_controller.o \
			open.o solver.o IO.o Lifecycle.o interfaces.o

./.lib/libvelocity.a: $(OBJ)
	@echo;	\
	rm -f ./.lib/libvelocity.a;	\
	$(AR_CMD) -q ./.lib/libvelocity.a ./*.o  $(FIELDS)/*.o $(DATA)/*.o ; 
	
	
$(DATA)/%.o:	$(DATA)/%.f90
	$(print_current_build);
	$(MUTE)$(FF) -o $@ $(FFLAGS) $(INCS) $<	;
	@if [ -z "$(MUTE)" ]; then echo; fi;
	
$(FIELDS)/%.o:	$(FIELDS)/%.f90
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
	rm -f $(FIELDS)/*.o; rm -f $(FIELDS)/*.mod;		\
	rm -f ./.lib/libvelocity.a 
