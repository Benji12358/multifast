include $(ENVIRONMENT_FILE)

FFLAGS= -c $(COMPILER_OPTIONS)
INCS 		= $(HDF5_INC)

SOURCES=.
LIBDIR=.lib
INCDIR=./

OBJ=HDF5_IO.o OBJ_reader.o VTK_mod.o

build: $(LIBDIR)/libIO.a

$(LIBDIR)/libIO.a: $(OBJ)
	@echo;
	@rm -f $(LIBDIR)/libIO.a;
	@$(AR_CMD) -q $(LIBDIR)/libIO.a ./*.o;


%.o:	$(SOURCES)/%.f90
	$(print_current_build);
	$(MUTE)$(FF) $(FFLAGS) $(INCS) $<;
	@if [ -z "$(MUTE)" ]; then echo; fi;
	
clean:
	$(MUTE)rm -f ./$(LIBDIR)/*.a; rm -f ./*.o; rm -f ./$(INCDIR)/*.mod ;