include $(ENVIRONMENT_FILE)


FFLAGS= -c $(COMPILER_OPTIONS) $(MISCELLANEOUS_INC)
INCS 		= $(SOLVERS_INC)

SOURCES=.
LIBDIR=.lib
INCDIR=./

DRP=$(SOURCES)/schemes3D/DRP
O2=$(SOURCES)/schemes3D/O2
CPT=$(SOURCES)/schemes3D/CPT

DRP_OBJ=DRP_D0s.o DRP_D1s.o DRP_D1c.o DRP_D2c.o
CPT_OBJ=CPT.o CPT_D0s.o CPT_D1s.o CPT_D1c.o CPT_D2c.o
O2_OBJ=O2_D0s.o O2_D1s.o O2_D1c.o O2_D2c.o

SCHEME1D_OBJ=schemes_settings.o staggered_schemes.o collocated_schemes.o interpol.o second_order_schemes.o

OBJ= $(SCHEME1D_OBJ) $(DRP_OBJ) $(CPT_OBJ) $(O2_OBJ) 

build: $(LIBDIR)/libschemes.a 

$(LIBDIR)/libschemes.a: $(OBJ)
	@echo;
	@rm -f $(LIBDIR)/libschemes.a;
	@$(AR_CMD) -q $(LIBDIR)/libschemes.a ./*.o;
	
	
%.o:	$(SOURCES)/%.f90
	$(print_current_build);
	$(MUTE)$(FF) $(FFLAGS) $(INCS) $<;
	@if [ -z "$(MUTE)" ]; then echo; fi;
	
	
%.o:	$(DRP)/%.f90
	$(print_current_build);
	$(MUTE)$(FF) $(FFLAGS) $(INCS) $<;
	@if [ -z "$(MUTE)" ]; then echo; fi;
	
	
%.o:	$(CPT)/%.f90
	$(print_current_build);
	$(MUTE)$(FF) $(FFLAGS) $(INCS) $<;
	@if [ -z "$(MUTE)" ]; then echo; fi;
	
	
%.o:	$(O2)/%.f90
	$(print_current_build);
	$(MUTE)$(FF) $(FFLAGS) $(INCS) $<;
	@if [ -z "$(MUTE)" ]; then echo; fi;
	
	

clean:
	$(MUTE)rm -f ./$(LIBDIR)/*.a; rm -f ./*.o; rm -f ./$(INCDIR)/*.mod ;


