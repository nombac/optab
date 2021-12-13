
SHELL = /bin/sh

OPENMPI = /opt/local
HDF5 = /opt/local

FC = $(OPENMPI)/bin/mpif90
FFLAGS  = -Wall -O3 -I$(HDF5)/include #-fbounds-check
LDFLAGS  = -L$(HDF5)/lib -lhdf5 -lhdf5_hl -lhdf5_fortran -lhdf5_hl_fortran

#FFLAGS += -DDIRNIST=\'$(DIRNIST)\'
#FFLAGS += -DDIRH5=\'$(DIRH5)\'
#FFLAGS += -DDIRTOPBASE=\'$(DIRTOPBASE)\'
#FFLAGS += -DDIRGFGAM=\'$(DIRGFGAM)\'

all: convert_lines_h5 convert_nist_h5 convert_cia_h5 preproc_hitran convert_topbase_h5 convert_gfgam_h5

convert_lines_h5: code_module.o const_module.o convert_lines_h5.o h5lx_module.o nist_module.o sort_module.o
	$(FC) -o $@ $^ $(LDFLAGS)
convert_nist_h5: convert_nist_h5.o nist_module.o
	$(FC) -o $@ $^ $(LDFLAGS)
convert_cia_h5: convert_cia_h5.o
	$(FC) -o $@ $^ $(LDFLAGS)
convert_topbase_h5: convert_topbase_h5.o const_module.o
	$(FC) -o $@ $^ $(LDFLAGS)
convert_gfgam_h5: convert_gfgam_h5.o const_module.o
	$(FC) -o $@ $^ $(LDFLAGS)
decompose_hitran: const_module.o decompose_hitran.o
	$(FC) -o $@ $^ $(LDFLAGS)
preproc_hitran: const_module.o preproc_hitran.o
	$(FC) -o $@ $^ $(LDFLAGS)
test: test.o h5lx_module.o
	$(FC) -o $@ $^ $(LDFLAGS)

test.o: h5lx_module.o
code_module.o:
convert_cia.o:
convert_lines_h5.o: code_module.o const_module.o h5lx_module.o nist_module.o sort_module.o
convert_nist_h5.o: nist_module.o
convert_topbase_h5.o: const_module.o
convert_gfgam_h5.o: const_module.o
decompose_hitran.o: const_module.o
preproc_hitran.o: const_module.o

code_module.o:
h5lx_module.o:
nist_module.o:
sort_module.o:

%.o: %.F90
	$(FC) -c $(FFLAGS) -o $@ $<

clean:
	-\rm *.o a.out *~ *.mod convert_lines_h5 convert_nist_h5 decompose_hitran convert_cia_h5 preproc_hitran convert_gfgam_h5 convert_topbase_h5 test
