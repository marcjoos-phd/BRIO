#===============================================================================
#======================== Benchmark for paRallel I/O ===========================
#===============================================================================
# Author: Marc B.R. Joos
#
# Created/last modified: jun 26, 2013/jul 9, 2013
#
# This file is distributed under GNU/GPL licence, 
# see <http://www.gnu.org/licenses/>.
#===============================================================================
# Define your environment & preprocessor options
# IOTYPE is in (POSIX, PNCDF, PHDF5, ALL)
INTERACTIVE=0
IOTYPE=PNCDF
ARCH=LOCAL

# IOTYPE definitions
#===============================================================================
ifeq ($(IOTYPE),POSIX)
POSIX=1
PNCDF=0
PHDF5=0
endif
ifeq ($(IOTYPE),PNCDF)
POSIX=0
PNCDF=1
PHDF5=0
endif
ifeq ($(IOTYPE),PHDF5)
POSIX=0
PNCDF=0
PHDF5=1
endif
ifeq ($(IOTYPE),ALL)
POSIX=1
PNCDF=1
PHDF5=1
endif

# On sappcm197
#===============================================================================
F90_LOCAL      = gfortran
FFLAGS_LOCAL   = -O3
CPPFLAGS_LOCAL = -x f95-cpp-input -DINTERACTIVE=$(INTERACTIVE) -DPOSIX=$(POSIX) -DPNCDF=$(PNCDF) -DPHDF5=$(PHDF5)
MPIINC_LOCAL   = -I/local/home/mjoos/soft/build/include/ 
MPILIB_LOCAL   = -L/local/home/mjoos/soft/build/lib -lmpi -lmpi_f77
HDFINC_LOCAL   = -I/local/home/mjoos/soft/hdf5-para/include/ 
HDFLIB_LOCAL   = -L/local/home/mjoos/soft/hdf5-para/lib -lhdf5_fortran -lhdf5 -lz
CDFINC_LOCAL   = -I/local/home/mjoos/soft/build/include/
CDFLIB_LOCAL   = -L/local/home/mjoos/soft/build/lib -lpnetcdf

# On Turing (BG/Q)
#============================================================================
F90_TURING      = mpixlf90_r
FFLAGS_TURING   = -O3
CPPFLAGS_TURING = -qsuffix=cpp=f90 -WF,-DINTERACTIVE=$(INTERACTIVE),-DPOSIX=$(POSIX),-DPNCDF=$(PNCDF),-DPHDF5=$(PHDF5)
MPIINC_TURING   =
MPILIB_TURING   =
HDFINC_TURING   = -I/bglocal/cn/pub/HDF5/1.8.9/par/include/
HDFLIB_TURING   = -L/bglocal/cn/pub/HDF5/1.8.9/par/lib -lhdf5_fortran -lhdf5
CDFINC_TURING   = -I/bglocal/cn/pub/Parallel-NetCDF/1.3.1/include/
CDFLIB_TURING   = -L/bglocal/cn/pub/Parallel-NetCDF/1.3.1/lib -lpnetcdf

# Gather all flags, libraries and objects; then compile
#===============================================================================
F90      = $(F90_$(ARCH))
FFLAGS   = $(FFLAGS_$(ARCH))
CPPFLAGS = $(CPPFLAGS_$(ARCH))
MPIINC   = $(MPIINC_$(ARCH))
MPILIB   = $(MPILIB_$(ARCH))
HDFINC   = $(HDFINC_$(ARCH))
HDFLIB   = $(HDFLIB_$(ARCH))
CDFINC   = $(CDFINC_$(ARCH))
CDFLIB   = $(CDFLIB_$(ARCH))

MODOBJ = commons.o
BRIOBJ = BRIO.o
ALLOBJ = $(MODOBJ) $(BRIOBJ)

#============================================================================
%.o:%.f90
	$(F90) $(HDFINC) $(CDFINC) $(FLAGS) $(MPIINC) $(CPPFLAGS) -c $^ -o $@
#============================================================================
BRIO: $(ALLOBJ)
	@echo "======================== BRIO =========================="
	@echo "============= Benchmark for parallel I/O ==============="
	@echo "======== Test parallel NetCDF and parallel HDF5 ========"
	@echo "========================================================"
	@echo "With: "
ifeq ($(IOTYPE),POSIX)
	@echo " - Sequential POSIX I/O"
endif
ifeq ($(IOTYPE),PNCDF)
	@echo " - Parallel NetCDF I/O"
endif
ifeq ($(IOTYPE),PHDF5)
	@echo " - Parallel HDF5 I/O"
endif
ifeq ($(IOTYPE),ALL)
	@echo " - Sequential POSIX I/O"
	@echo " - Parallel NetCDF I/O"
	@echo " - Parallel HDF5 I/O"
endif
ifeq ($(INTERACTIVE),1)
	@echo " - interactive interface"
else
	@echo " - namelist interface"
endif
	@echo "========================================================"
	@echo 
	$(F90) $(FLAGS) -o ./BRIO $(ALLOBJ) $(HDFLIB) $(CDFLIB) $(MPILIB)
all: BRIO
#============================================================================
clean:
	rm *.o *.mod *.h5 *.nc 
	rm -rf sequentialio/