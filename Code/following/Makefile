include $(ENVIRONMENT_FILE)

FFLAGS		= -c $(COMPILER_OPTIONS) $(DECOMP2D_INC)
SOURCES		= sources
INCS 		= $(FLAG_FFT) $(HDF5_INC) $(IO_INC) $(SOLVERS_INC) $(MISCELLANEOUS_INC) $(SCHEMES_INC) -I $(CORE) -I $(COMMON)

DATA		=./data
MESH		=./mesh
FIELDS		=./fields_utils


OBJ = $(DATA)/following_data.o $(DATA)/workspace_view.o $(MESH)/mesh_settings.o $(MESH)/mesh_generator.o $(FIELDS)/operations.o $(FIELDS)/quality_criteriums.o \
		IO.o following_vel_BC_controller.o \
		poisson_solver.o passive_scalar_solver.o fringe_solver.o vel_solver.o \
		following_vel_initializer.o following_fringe_initializer.o following_sca_initializer.o \
		Lifecycle.o initializer.o interfaces.o 


./.lib/libfollowing.a: $(OBJ)
	@echo;	\
	rm -f ./.lib/libfollowing.a;	\
	$(AR_CMD) -q ./.lib/libfollowing.a ./*.o $(DATA)/*.o ; 
	
	
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
	rm -f ./.lib/libfollowing.a 
