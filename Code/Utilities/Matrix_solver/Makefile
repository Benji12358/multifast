include $(ENVIRONMENT_FILE)

FFLAGS= -c $(COMPILER_OPTIONS)

SOURCES=.
LIBDIR=.lib
INCDIR=./

build: $(LIBDIR)/libsolvers.a

OBJ=solvers.o

$(LIBDIR)/libsolvers.a: $(OBJ)
	@echo;
	@rm -f $(LIBDIR)/libsolvers.a;
	@$(AR_CMD) -q $(LIBDIR)/libsolvers.a ./*.o;

%.o:	$(SOURCES)/%.f90
	$(print_current_build);
	$(MUTE)$(FF) $(FFLAGS) $(INCS) $<;
	@if [ -z "$(MUTE)" ]; then echo; fi;
	
clean:
	$(MUTE)rm -f ./$(LIBDIR)/*.a ; rm -f ./*.o; rm -f ./$(INCDIR)/*.mod ;
