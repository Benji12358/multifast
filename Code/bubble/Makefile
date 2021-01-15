include $(ENVIRONMENT_FILE)

FFLAGS		= -c $(COMPILER_OPTIONS) $(DECOMP2D_INC)
SOURCES		= sources
INCS 		= $(FLAG_FFT) $(HDF5_INC) $(IO_INC) $(SOLVERS_INC) $(MISCELLANEOUS_INC) $(SCHEMES_INC) -I $(COMMON)

OBJ=bubble_data.o fields_utils.o Bubble_generator.o 


./.lib/libbubble.a: $(OBJ)
	@echo;	\
	rm -f ./.lib/libbubble.a;	\
	$(AR_CMD) -q ./.lib/libbubble.a ./*.o ; 
	
%.o:	%.f90
	$(print_current_build);
	$(MUTE)$(FF) $(FFLAGS) $(INCS) $<;
	@if [ -z "$(MUTE)" ]; then echo; fi;
	
clean :
	$(MUTE)rm -f ./*.o; rm -f ./*.mod;		\
	rm -f ./.lib/libbubble.a 
