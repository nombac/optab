### SYSTEM DEPENDENT SECTION ###
OPENMPI = /opt/local
HDF5 = /opt/local

FFLAGS  = -fallow-argument-mismatch -Wall -O3 #-fbounds-check

### SYSTEM INDEPENDENT SECTION
SHELL = /bin/sh

H5PFC = false
ifeq ($(H5PFC),true)
FC = $(HDF5)/bin/h5pfc
else
FC = $(OPENMPI)/bin/mpif90
FFLAGS += -I$(OPENMPI)/include -I$(HDF5)/include
LDFLAGS = -L$(HDF5)/lib -lhdf5 -lhdf5_hl -lhdf5_fortran -lhdf5_hl_fortran
endif

#LDFLAGS += -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/lib
#LDFLAGS += -mmacosx-version-min=11.2

#--------------------------------


### SYSTEM INDEPENDENT SECTION ###

# Check: REAL128 (quad precision) is not supported on arm64 macOS
ifneq ($(MAKECMDGOALS),clean)
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
ifeq ($(UNAME_S),Darwin)
ifeq ($(UNAME_M),arm64)
$(error ERROR: REAL128 (quad precision) is not supported on arm64 macOS. Use an x86_64 machine to compile and run both this program (prep_FastChem) and FastChem itself.)
endif
endif
endif

prep_FastChem: prep_FastChem.o const_module.o
	$(FC) -o $@ $^ $(LDFLAGS)

const_module.o:

prep_FastChem.o: const_module.o

%.o: %.F90
	$(FC) -c $(FFLAGS) -o $@ $<

clean:
	-\rm -rf *.o *~ *.mod prep_FastChem config.input* lnK_* tpgrid_*.dat
