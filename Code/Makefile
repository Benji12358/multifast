include $(ENVIRONMENT_FILE)


CORE:=$(shell pwd)
COMMON				= $(CORE)/common
COMMON_LIB			=-L $(COMMON)/.lib -lcommon
SCALAR				= $(CORE)/passive_scalar
SCALAR_LIB			=-L $(SCALAR)/.lib -lscalar
MHD					= $(CORE)/MHD
MHD_LIB				=-L $(MHD)/.lib -lmhd
FRINGE				= $(CORE)/FRINGE
FRINGE_LIB			=-L $(FRINGE)/.lib -lfringe
BUBBLE				= $(CORE)/bubble
BUBBLE_LIB			=-L $(BUBBLE)/.lib -lbubble
VELOCITY			= $(CORE)/velocity
VELOCITY_LIB		=-L $(VELOCITY)/.lib -lvelocity
PHASE_AVG           = $(CORE)/phase_average
PHASE_AVG_LIB       =-L $(PHASE_AVG)/.lib -lphaseaverage
EMBEDDED            = $(CORE)/embedded
EMBEDDED_LIB        =-L $(EMBEDDED)/.lib -lembedded

# Flags _________________________________________________________________
ifeq ($(COMPILER), mpixlf95_r)
	LIBS =  $(IO_LIB) $(HDF5_LIB) $(SCHEMES_LIB) $(SOLVERS_LIB) $(MISCELLANEOUS_LIB) $(LIBS_FFT) $(DECOMP2D_LIB) $(LAPACK_LIB) 
endif

ifeq ($(COMPILER),ifort)
	LIBS =  $(SOLVERS_LIB) $(MISCELLANEOUS_LIB) $(LIBS_FFT) $(SCHEMES_LIB) $(LIBS_MKL) $(DECOMP2D_LIB)
endif

ifeq ($(COMPILER),gfortran)
	LIBS 		= $(SCHEMES_LIB) $(SOLVERS_LIB) $(MISCELLANEOUS_LIB) $(IO_LIB) $(HDF5_LIB) $(LIBS_FFT) $(DECOMP2D_LIB) $(LAPACK_LIB)
endif



FFLAGS		= -c $(COMPILER_OPTIONS) $(DECOMP2D_INC)
LDFLAG		= $(COMPILER_OPTIONS) -o
DATA		= $(CORE)//common/data

MESH		= $(CORE)/common/mesh

INCS 		= $(FLAG_FFT) $(HDF5_INC) $(IO_INC) $(SOLVERS_INC) $(MISCELLANEOUS_INC) $(SCHEMES_INC) -I $(COMMON) -I $(SCALAR) -I $(MHD) -I $(FRINGE) -I $(BUBBLE) -I $(VELOCITY) -I $(PHASE_AVG) -I $(EMBEDDED)

OBJ		= main.o diagnostic.o multiphysics.o IO_handler.o \
		   


# Dependencies __________________________________________________________

#$(SCALAR)/scalar_data.o 	: $(COMMON)/mesh_settings.o 


main.o			: build_common build_scalar build_mhd build_fringe build_bubble build_velocity build_embedded build_phaseaverage build_core   

# Rules _________________________________________________________________
$(EXEC_NAME): subsection_text=Linking objects files
$(EXEC_NAME):	main.o    
	@$(print_subsection);
	$(FF) $(LDFLAG) $(EXEC_NAME) $(OBJ)  $(BUBBLE_LIB) $(SCALAR_LIB) $(MHD_LIB) $(FRINGE_LIB) $(VELOCITY_LIB) $(EMBEDDED_LIB) $(PHASE_AVG_LIB) $(COMMON_LIB) $(LIBS);
	@echo; echo;
	
build: $(EXEC_NAME)

clean:
	$(MUTE)rm -f ./$(EXEC_NAME);							\
	rm -f ./*.o; rm -f ./*.mod; 						\
	make -C $(COMMON) clean ;					\
	make -C $(SCALAR) clean ;					\
	make -C $(MHD) clean ;					\
	make -C $(FRINGE) clean ;					\
	make -C $(PHASE_AVG) clean ; \
	make -C $(EMBEDDED) clean ; \
	make -C $(VELOCITY) clean ;					\
	make -C $(BUBBLE) clean ;
	
build_common: subsection_text=Building common module
build_common :	
	@$(print_subsection);
	$(MUTE)make -C $(COMMON) ENVIRONMENT_FILE=$(ENVIRONMENT_FILE) ;
	
build_core: core_section IO_handler.o diagnostic.o multiphysics.o

core_section :	subsection_text=Building core module
core_section : 
	@$(print_subsection);
	
build_bubble :	subsection_text=Building bubble module
build_bubble :	
	@$(print_subsection);
	$(MUTE)make -C $(BUBBLE) ENVIRONMENT_FILE=$(ENVIRONMENT_FILE) COMMON=$(COMMON) ;
	
build_scalar :	subsection_text=Building scalar module
build_scalar :	
	@$(print_subsection);
	$(MUTE)make -C $(SCALAR) ENVIRONMENT_FILE=$(ENVIRONMENT_FILE) COMMON=$(COMMON) ;
	
build_mhd :	subsection_text=Building MHD module
build_mhd :	
	@$(print_subsection);
	$(MUTE)make -C $(MHD) ENVIRONMENT_FILE=$(ENVIRONMENT_FILE) COMMON=$(COMMON) ;
	
build_fringe :	subsection_text=Building FRINGE module
build_fringe :	
	@$(print_subsection);
	$(MUTE)make -C $(FRINGE) ENVIRONMENT_FILE=$(ENVIRONMENT_FILE) COMMON=$(COMMON) CORE=$(CORE) ;
	
build_velocity :	subsection_text=Building velocity module
build_velocity:
	@$(print_subsection);
	$(MUTE)make -C $(VELOCITY) ENVIRONMENT_FILE=$(ENVIRONMENT_FILE) COMMON=$(COMMON) CORE=$(CORE) ;

build_phaseaverage :	subsection_text=Building phase average module
build_phaseaverage :	
	@$(print_subsection);
	$(MUTE)make -C $(PHASE_AVG) ENVIRONMENT_FILE=$(ENVIRONMENT_FILE) COMMON=$(COMMON) CORE=$(CORE);

build_embedded :	subsection_text=Building embedded module
build_embedded :	
	@$(print_subsection);
	$(MUTE)make -C $(EMBEDDED) ENVIRONMENT_FILE=$(ENVIRONMENT_FILE) COMMON=$(COMMON) CORE=$(CORE);

	
%.o:	%.f90
	$(print_current_build);
	$(MUTE)$(FF) $(FFLAGS) $(INCS) $<
