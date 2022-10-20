include $(ENVIRONMENT_FILE)

# Flags _________________________________________________________________
ifeq ($(COMPILER), mpixlf95_r)
	LIBS =  $(IO_LIB) $(HDF5_LIB) $(SCHEMES_LIB) $(SOLVERS_LIB) $(MISCELLANEOUS_LIB) $(LIBS_FFT) $(DECOMP2D_LIB) $(LAPACK_LIB) 
endif

ifeq ($(COMPILER),ifort)
	LIBS =  $(SOLVERS_LIB) $(MISCELLANEOUS_LIB) $(LIBS_FFT) $(SCHEMES_LIB) $(LIBS_MKL) $(DECOMP2D_LIB)
endif

ifeq ($(COMPILER),gfortran)
	LIBS 		= $(IO_LIB) $(HDF5_LIB) $(SCHEMES_LIB) $(SOLVERS_LIB) $(MISCELLANEOUS_LIB) $(LIBS_FFT) $(DECOMP2D_LIB) $(LAPACK_LIB)
endif

FFLAGS		= -c $(COMPILER_OPTIONS) $(DECOMP2D_INC)
INCS 		= $(FLAG_FFT) $(HDF5_INC) $(IO_INC) $(SOLVERS_INC) $(MISCELLANEOUS_INC) $(SCHEMES_INC)

DATA=./data
MESH=./mesh
UTILS=./utils

OBJ=$(DATA)/run_ctxt.o $(DATA)/anim2D.o $(DATA)/workspace_view.o $(DATA)/boundaries.o $(DATA)/DNS_settings.o $(DATA)/IBM_data.o $(DATA)/IBM_settings.o $(MESH)/mesh_settings.o $(MESH)/mesh_generator.o $(DATA)/physical_fields.o $(DATA)/schemes_interface.o \
			Fields_tools.o $(UTILS)/buffers2d.o subdomains_view.o IO.o schemes_loader.o Poisson_generic_solver.o DRP_IBM.o DRP_IBM_fine.o IBM.o Poisson_solver.o


./.lib/libcommon.a 	: $(OBJ)
			@echo; \
			rm -f ./.lib/libcommon.a;	\
			$(AR_CMD) -q ./.lib/libcommon.a ./*.o $(DATA)/*.o $(MESH)/*.o $(UTILS)/*.o ; 	
			
$(DATA)/%.o:	$(DATA)/%.f90
	$(print_current_build);
	$(MUTE)$(FF) -o $@ $(FFLAGS) $(INCS) $< ;
	@if [ -z "$(MUTE)" ]; then echo; fi;
	
$(MESH)/%.o:	$(MESH)/%.f90
	$(print_current_build);
	$(MUTE)$(FF) -o $@ $(FFLAGS) $(INCS) $< ;
	@if [ -z "$(MUTE)" ]; then echo; fi;
	
$(UTILS)/%.o:	$(UTILS)/%.f90
	$(print_current_build);
	$(MUTE)$(FF) -o $@ $(FFLAGS) $(INCS) $< ;
	@if [ -z "$(MUTE)" ]; then echo; fi;
	
%.o:	%.f90
	$(print_current_build);
	$(MUTE)$(FF) $(FFLAGS) $(INCS) $< ;
	@if [ -z "$(MUTE)" ]; then echo; fi;
	
clean :
	@ rm -f $(DATA)/*.o; rm -f $(DATA)/*.mod ;	\
	rm -f $(MESH)/*.o;	rm -f $(MESH)/*.mod ;		\
	rm -f $(UTILS)/*.o;	rm -f $(UTILS)/*.mod ;		\
	rm -f ./*.o;  rm -f ./*.mod;						\
	rm -f ./.lib/libcommon.a 
