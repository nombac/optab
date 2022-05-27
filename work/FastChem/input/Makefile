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

logK_FastChem: const_module.o logK_FastChem.o
	$(FC) -o $@ $^ $(LDFLAGS)

table: table.o
	$(FC) -o $@ $^ $(LDFLAGS)

const_module.o:

logK_FastChem.o: const_module.o

%.o: %.F90
	$(FC) -c $(FFLAGS) -o $@ $<

clean:
	-\rm *.o *~ *.mod logK_FastChem table
