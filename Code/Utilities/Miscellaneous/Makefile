include $(ENVIRONMENT_FILE)

FFLAGS= -c $(COMPILER_OPTIONS)
INCS 		= $(DECOMP2D_INC)

SOURCES=.
LIBDIR=.lib
INCDIR=./

build: $(LIBDIR)/libmiscellaneous.a


OBJ=boundaries_types.o mathematical_constants.o file_copy.o formatter.o decomp2D_utils.o

$(LIBDIR)/libmiscellaneous.a: $(OBJ)
	@echo;
	@rm -f $(LIBDIR)/libmiscellaneous.a;
	@$(AR_CMD) -q $(LIBDIR)/libmiscellaneous.a boundaries_types.o  mathematical_constants.o file_copy.o decomp2D_utils.o formatter.o;

%.o:	$(SOURCES)/%.f90
	$(print_current_build);
	$(MUTE)$(FF) $(FFLAGS) $(INCS) $<;
	@if [ -z "$(MUTE)" ]; then echo; fi;

clean:
	$(MUTE)rm -f ./$(LIBDIR)/*.a; rm -f ./*.o; rm -f ./$(INCDIR)/*.mod ;
