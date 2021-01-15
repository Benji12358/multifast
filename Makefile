include ./Makefile.inc
include $(COMPILER_OPTIONS_FILE)


EXEC_NAME="multiFAST"
SIM_EXEC="DNS_EXEC"
TEST_EXEC="TEST_EXEC"
CODE_DIR:=./Code

KILL_SCRIPT=$(SCRIPTS_DIR)/killtest.sh
TEST_SCRIPT=$(SCRIPTS_DIR)/test.sh
RUN_SCRIPT=$(SCRIPTS_DIR)/run.sh

SMSTO:=""
MAILTO:=""

BM:=SINGLETON
PWD:=$(shell pwd)


####################################################################################################################################
####################################################################################################################################
####################################                                       #########################################################
####################################        RULES FOR BUILDING CODE 	   #########################################################
####################################                                       #########################################################
####################################################################################################################################
####################################################################################################################################
	

build: builds_tools build_multifast
	@cp ./Code/$(EXEC_NAME) ./$(SIM_EXEC) ;	\
	cp ./Code/$(EXEC_NAME) ./$(TEST_EXEC) ;


####################################################################################################################################
#################################### BUILD TOOLS RULES #############################################################################
####################################################################################################################################
builds_tools: builds_tools_header build_miscellaneous build_solvers build_schemes build_io 

builds_tools_header: section_text=BUILD TOOLS LIBRARIES
builds_tools_header:
	$(print_section)

build_miscellaneous :	subsection_text=Building miscellaneous library
build_miscellaneous:
	@$(print_subsection);																	\
	make -C $(LOC_MISCELLANEOUS) build ENVIRONMENT_FILE=$(ENVIRONMENT_FILE);

build_solvers :	subsection_text=Building solvers library
build_solvers:
	@$(print_subsection);													 				\
	make -C $(LOC_SOLVERS) build ENVIRONMENT_FILE=$(ENVIRONMENT_FILE);

build_schemes :	subsection_text=Building schemes library
build_schemes:
	@$(print_subsection);														 			\
	make -C $(LOC_SCHEMES) build ENVIRONMENT_FILE=$(ENVIRONMENT_FILE) ;

build_io :	subsection_text=Building IO library
build_io :
	@$(print_subsection);																	\
	make -C $(LOC_IO) build ENVIRONMENT_FILE=$(ENVIRONMENT_FILE) ;	


####################################################################################################################################
########################################### BUILD MULTIFAST RULE ###################################################################
####################################################################################################################################
build_multifast: section_text="BUILD MULTIFAST"
build_multifast:
	@$(print_section);	\
	make -C Code/ build EXEC_NAME=$(EXEC_NAME) ENVIRONMENT_FILE=$(ENVIRONMENT_FILE);



####################################################################################################################################
############################################# CLEAN RULES ##########################################################################
####################################################################################################################################


clean:
	make -C Code/ clean EXEC_NAME=$(EXEC_NAME);	

clean_all:
	make -C $(LOC_MISCELLANEOUS) clean;
	make -C $(LOC_SOLVERS) clean;
	make -C $(LOC_SCHEMES) clean;
	make -C $(LOC_IO) clean;
	make -C Code/ clean EXEC_NAME=$(EXEC_NAME);

