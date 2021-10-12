
####################################################################################################################################
####################################################################################################################################
####################################      BUILDING AND SIMULATIONS         #########################################################
####################################        ENVIRONMENT VARIABLES          ########################################################
####################################                                       #########################################################
####################################################################################################################################
####################################################################################################################################



####################################################################################################################################
############################### Set the environment informations HERE ##############################################################
####################################################################################################################################

# In most case, you only have to edit this section


ROOT				    =/home/users/arrondea7b/WORKSPACE
CODES_PATH			    =$(ROOT)/Codes
MULTIFAST_PATH			=$(CODES_PATH)/DNS/MULTIFAST_MHD_v2_IBM
ENVIRONMENT_FILE		=$(MULTIFAST_PATH)/Makefile.inc


COMPILER				=gfortran
COMPILER_OPTIONS_FILE	=$(MULTIFAST_PATH)/compiler_options/LEGI/$(COMPILER).conf
COMPIL_OPT_LEVEL		= HIGH
COMPIL_CHECK_LEVEL		= ZERO

# Default number of processes for a simulation NB_PROCS=PROW*PCOL (set PCOL=PROW=0 to let decomp2D to chose the decomposition)
NB_PROCS=2
PROW=0
PCOL=0

# Comment this line to see details at compilation time
# MUTE=@


####################################################################################################################################
############################### Definition of working path from above informations #################################################
##########################################  TO NOT CHANGE   ########################################################################


TMPDIR=$(MULTIFAST_PATH)/TMP
SIMULATION_DIR=$(MULTIFAST_PATH)/Simulations


SYSTEM_INFO_FILE=$(CODES_PATH)/system.info

EXTERN_LIBS=$(CODES_PATH)/libs.$(COMPILER)



# Utilities paths *************************************************************************************


UTILITIES=$(MULTIFAST_PATH)/Code/Utilities

LOC_SOLVERS= $(UTILITIES)/Matrix_solver
LOC_SCHEMES=$(UTILITIES)/Schemes
LOC_IO= $(UTILITIES)/IO
LOC_MISCELLANEOUS=$(UTILITIES)/Miscellaneous

SOLVERS_INC= -I $(UTILITIES)/Matrix_solver
SOLVERS_LIB= -L $(UTILITIES)/Matrix_solver/.lib -lsolvers

SCHEMES_INC=-I $(UTILITIES)/Schemes
SCHEMES_LIB=-L $(UTILITIES)/Schemes/.lib -lschemes

IO_INC= -I $(UTILITIES)/IO
IO_LIB= -L $(UTILITIES)/IO/.lib -lIO

MISCELLANEOUS_INC=-I $(UTILITIES)/Miscellaneous
MISCELLANEOUS_LIB=-L $(UTILITIES)/Miscellaneous/.lib -lmiscellaneous


include $(SYSTEM_INFO_FILE)
include $(COMPILER_OPTIONS_FILE)
include $(EXTERN_LIBS)

####################################################################################################################################
####################################################################################################################################
####################################                                       #########################################################
####################################        SHARED SCRIPT OR               ########################################################
####################################               MAKEFILE FUNCTIONS      #########################################################
####################################################################################################################################
####################################################################################################################################

SCRIPTS_DIR=$(MULTIFAST_PATH)/scripts

JUSTIFY=$(SCRIPTS_DIR)/justify.sh
TEST_HANDLER=$(SCRIPTS_DIR)/texthandler.sh

ECHO = @/bin/echo -e
IMPORTANT_TEXT=\e[1;34m
NORMAL_TEXT=\e[0m



define print_section
@$(TEST_HANDLER) $(section_text)
endef


define print_subsection
@$(ECHO) "\n\n$(IMPORTANT_TEXT)" `$(JUSTIFY) "$(subsection_text) " ':' 98 ":"` "$(NORMAL_TEXT)"
endef


define print_current_build
@$(JUSTIFY) "Build $@" "from $<" 50
endef
