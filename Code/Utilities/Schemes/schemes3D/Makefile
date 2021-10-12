FF=gfortran

DRP_SCHEMES=./DRP
O2_SCHEMES=./O2
CPT_SCHEMES=./CPT
TURING_ROOT=/data2/TURING
CODES_PATH=$(TURING_ROOT)/Codes
DNS_PATH=$(CODES_PATH)/DNS/DNS4.0
UTILITIES=$(DNS_PATH)/Code/Utilities
MISCELLANEOUS_INC=-I $(UTILITIES)/Miscellaneous/include

OBJ=check.o DRP_D0s.o DRP_D1s.o DRP_D1c.o DRP_D2c.o DRP_IBM.o O2_D0s.o O2_D1s.o O2_D1c.o O2_D2c.o \
CPT.o CPT_D0s.o CPT_D1s.o CPT_D1c.o CPT_D2c.o

all: check

check0: check.o $(DRP_SCHEMES)/DRP_D0s.o $(DRP_SCHEMES)/DRP_D1s.o
	$(FF) check.o $(DRP_SCHEMES)/DRP_D0s.o $(DRP_SCHEMES)/DRP_D1s.o

check.o: $(DRP_SCHEMES)/DRP_D0s.o $(DRP_SCHEMES)/DRP_D1s.o $(DRP_SCHEMES)/DRP_D1c.o $(DRP_SCHEMES)/DRP_D2c.o $(O2_SCHEMES)/O2_D0s.o $(O2_SCHEMES)/O2_D1s.o $(O2_SCHEMES)/O2_D1c.o $(O2_SCHEMES)/O2_D2c.o \
$(CPT_SCHEMES)/CPT.o $(CPT_SCHEMES)/CPT_D0s.o $(CPT_SCHEMES)/CPT_D1s.o $(CPT_SCHEMES)/CPT_D1c.o $(CPT_SCHEMES)/CPT_D2c.o


check:	check.o    
	$(FF) -o check $(OBJ) $(LIBS)		

	
$(DRP_SCHEMES)/%.o:	$(DRP_SCHEMES)/%.f90
	$(FF) -c -ffree-line-length-0 $(MISCELLANEOUS_INC)  $<	

	
$(O2_SCHEMES)/%.o:	$(O2_SCHEMES)/%.f90
	$(FF) -c -ffree-line-length-0 $(MISCELLANEOUS_INC)  $<	

	
$(CPT_SCHEMES)/%.o:	$(CPT_SCHEMES)/%.f90
	$(FF) -c -ffree-line-length-0 $(MISCELLANEOUS_INC)  $<


%.o:	%.f90
	$(FF) -c -ffree-line-length-0 $(INCS) $<


clean:
	rm ./*.o;		\
	rm ./*.mod
