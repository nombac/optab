SHELL = /bin/sh

HDF5 = /opt/local
OPENMPI = /opt/local
FFLAGS = -Wall #-fbounds-check

FC = $(OPENMPI)/bin/mpif90
FFLAGS  += -O3 -I$(HDF5)/include
LDFLAGS  = -L$(HDF5)/lib -lhdf5 -lhdf5_hl -lhdf5_fortran -lhdf5_hl_fortran
BIN = .

all: convert_Phoenix convert_TEA convert_FastChem

convert_Phoenix: const_module.o convert_Phoenix.o string_module.o
	$(FC) -o $(BIN)/$@ $^ $(LDFLAGS)

convert_TEA: const_module.o convert_TEA.o string_module.o
	$(FC) -o $(BIN)/$@ $^ $(LDFLAGS)

convert_FastChem: const_module.o convert_FastChem.o string_module.o
	$(FC) -o $(BIN)/$@ $^ $(LDFLAGS)

list_id: list_id.o
	$(FC) -o $(BIN)/$@ $^ $(LDFLAGS)

const_module.o:

convert_Phoenix.o: const_module.o string_module.o

convert_TEA.o: const_module.o string_module.o

convert_FastChem.o: const_module.o string_module.o

list_id.o:

string_module.o:

%.o: %.F90
	$(FC) -c $(FFLAGS) -o $@ $<

clean:
	-\rm *.o a.out *~ *.mod list_id convert_Phoenix convert_TEA convert_FastChem
