### SYSTEM DEPENDENT SECTION ###
OPENMPI = /opt/local
HDF5 = /opt/local

FFLAGS  = -Wall -O3 #-fbounds-check

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


### OPTIONS FOR SPECIFIC PURPOSES ####
# output absorption coefficient per number density [cm^2] (rather than absorption coefficient [cm^-1]
SUFFIX :=

DEBUG = false
ifeq ($(DEBUG),true)
FFLAGS += -DDEBUG
SUFFIX := $(join $(SUFFIX),_debug)
endif

NLINES = false
ifeq ($(NLINES),true)
FFLAGS += -DLINE_STATISTICS
SUFFIX := $(join $(SUFFIX),_nlines)
endif

CROSS = false
ifeq ($(CROSS),true)
FFLAGS += -DCROSS_CHECK
SUFFIX := $(join $(SUFFIX),_cross)
endif

# MPIO
MPIO = false
ifeq ($(MPIO),true)
FFLAGS += -DMPIO
SUFFIX := $(join $(SUFFIX),_mpio)
endif

# use strength cutoff
STRENGTH = false
ifeq ($(STRENGTH),true)
FFLAGS += -DSTRENGTH_CUTOFF
SUFFIX := $(join $(SUFFIX),_strength)
endif

# use verner95 only for ground state photoionization (rather than hybrid of verner95 and verner96)
VERNER95 = false
ifeq ($(VERNER95),true)
FFLAGS += -DVERNER95_ONLY
SUFFIX := $(join $(SUFFIX),_verner95)
endif
# use verner96 only for ground state photoionization (rather than hybrid of verner95 and verner96)
VERNER96 = false
ifeq ($(VERNER96),true)
FFLAGS += -DVERNER96_ONLY
SUFFIX := $(join $(SUFFIX),_verner96)
endif

# use Phoenix grid
PHOENIXCOMP = false
ifeq ($(PHOENIXCOMP),true)
FFLAGS += -DPHOENIX_GRID
FFLAGS += -DPHOENIX_BREMS
FFLAGS += -DWAVELENGTH_BASE
FFLAGS += -DPHOENIX_BROADENING
FFLAGS += -DUNUSE_GFGAM
SUFFIX := $(join $(SUFFIX),_phoenixcomp)
endif
PHOENIXCOMPMIN = false
ifeq ($(PHOENIXCOMPMIN),true)
FFLAGS += -DPHOENIX_GRID
SUFFIX := $(join $(SUFFIX),_phoenixcomp_min)
endif

#--------------------------------


### SYSTEM INDEPENDENT SECTION ###
EXEC = a.out$(SUFFIX)

$(EXEC): code_module.o const_module.o continuum_module.o main.o h5lx_module.o input_module.o line_module.o mean_module.o mpi_module.o output_module.o voigt_module.o gffree_module.o
	$(FC) -o $@ $^ $(LDFLAGS)

code_module.o:

const_module.o:

continuum_module.o: code_module.o const_module.o mpi_module.o gffree_module.o

gffree_module.o:

h5lx_module.o:

input_module.o: mpi_module.o

line_module.o: code_module.o const_module.o h5lx_module.o mpi_module.o voigt_module.o

main.o: const_module.o continuum_module.o input_module.o line_module.o h5lx_module.o mean_module.o mpi_module.o output_module.o

mean_module.o: const_module.o mpi_module.o

mpi_module.o:

output_module.o: h5lx_module.o mpi_module.o

voigt_module.o:

%.o: %.F90
	$(FC) -c $(FFLAGS) -o $@ $<

clean:
	-\rm *.o *~ *.mod a.out*
