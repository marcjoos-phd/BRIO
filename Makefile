INTERACTIVE=1
POSIX=1
POTOK=1
PNCDF=1
PHDF5=1
ADIOS=1
MPIIO=1

F90      = mpif90
FFLAGS   = -O3 -fbacktrace -g -ffree-line-length-0
CPPFLAGS = -x f95-cpp-input -DINTERACTIVE=$(INTERACTIVE) -DPOSIX=$(POSIX) -DPOTOK=$(POTOK) -DPNCDF=$(PNCDF) -DPHDF5=$(PHDF5) -DADIOS=$(ADIOS) -DMPIIO=$(MPIIO) 
HDFINC   = -I/local/home/mjoos/soft/hdf5-para/include/
HDFLIB   = -L/local/home/mjoos/soft/hdf5-para/lib -lhdf5_fortran -lhdf5 -lz
CDFINC   = -I/local/home/mjoos/soft/build/include/
CDFLIB   = -L/local/home/mjoos/soft/build/lib -lpnetcdf
ifeq ($(ADIOS),1)
ADIOSDIR = /local/home/mjoos/soft/build
ADIOSINC = $(shell ${ADIOSDIR}/bin/adios_config -c -f)
ADIOSLIB = $(shell ${ADIOSDIR}/bin/adios_config -l -f)
endif

VPATH  = src:src/modules
SRCOBJ = $(wildcard src/modules/*f90)
MODOBJ = $(SRCOBJ:.f90=.o)
BRIOBJ = BRIO.o
ALLOBJ = $(MODOBJ) $(BRIOBJ)

#============================================================================
$(info ======================== BRIO ==========================)
$(info ============= Benchmark for parallel I/O ===============)
$(info ========================================================)
$(info With: )
ifeq ($(POSIX),1)
$(info  - Sequential POSIX I/O)
endif
ifeq ($(POTOK),1)
$(info  - Sequential POSIX I/O with token management)
endif
ifeq ($(PNCDF),1)
$(info  - Parallel NetCDF I/O)
endif
ifeq ($(PHDF5),1)
$(info  - Parallel HDF5 I/O)
endif
ifeq ($(ADIOS),1)
$(info  - Adaptive I/O System)
endif
ifeq ($(MPIIO),1)
$(info  - MPI-IO)
endif
$(info ========================================================)
$(info  )
#============================================================================
%.o:%.f90
	$(F90) $(HDFINC) $(CDFINC) $(ADIOSINC) $(FLAGS) $(MPIINC) $(CPPFLAGS) -c $^ -o $@
#============================================================================
BRIO: adios $(ALLOBJ)
	$(F90) $(FLAGS) -o ./BRIO $(ALLOBJ) $(HDFLIB) $(CDFLIB) $(MPILIB) $(ADIOSLIB)
	mv *o *mod BRIO bin/
#============================================================================
ifeq ($(ADIOS),1)
adios:
	gpp.py input/adios_BRIO.xml
	mv *fh src/modules/
else
adios:
endif
#============================================================================
all: BRIO
#============================================================================
clean:
	rm -rf bin/*.o src/modules/*.o src/modules/*.fh bin/*.mod input/*.h5 input/*.nc input/*.bp input/*.mp input/sequentialio/

