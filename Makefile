#===============================================================================
# 1. Compiler Settings
#===============================================================================
#--------------------------------------------------------------
# Basic Fortran compiler arguments 
#--------------------------------------------------------------
FC = h5pfc -cpp -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -fallow-argument-mismatch
FC += -Ofast
# FC += -O0 -g -fbacktrace -Warray-bounds

#--------------------------------------------------------------
# Advanced Fortran compiler options 
#--------------------------------------------------------------
#FC += -fma -finline-functions
#FC += -align array64byte
#FC += -heap-arrays -ip -fno-alias -xHost
#FC += -xCORE-AVX512 -axCORE-AVX512 -mtune=skylake
#FC += -axCORE-AVX512 
#FC += -xAVX -axCORE-AVX2

#--------------------------------------------------------------
# System specific libaries
#--------------------------------------------------------------
# CARTESIUS #
#FFTW3_LIB = -L/nfs/admin/hpc/sw/fftw3-3.3.3-intel-impi/lib

# GALILEO   #
#FFTW3_LIB = -L/cineca/prod/opt/libraries/fftw/3.3.8/intelmpi--2018--binary/lib

# IRENE     #
#FFTW3_LIB = -L/ccc/products/python-2.7.14/intel--17.0.4.196__openmpi--2.0.2/default/lib
#BLAS_LIB = -L/ccc/products/mkl-19.0.5.281/intel--19.0.5.281__openmpi--4.0.1/default/19.0.5.281/mkl/lib/intel64

# MareNos4  #
#FFTW3_LIB = -L/apps/FFTW/3.3.6/INTEL/IMPI/lib
#BLAS_LIB = -L/apps/INTEL/2018.4.057/mkl/lib/intel64

# DISCOVERER
# FFTW3_LIB = -L/opt/software/fftw/3/3.3.10-gcc-openmpi/lib/
# BLAS_LIB = -L/opt/software/lapack/3/3.11.0-gcc/lib64/

# Local
FFTW3_LIB = -L/usr/local/lib/

# Common build flags
FFTW3_FLAGS = -lfftw3
BLAS_FLAGS = -llapack
# BLAS_FLAGS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread 
# HDF5_FLAGS = -lhdf5_fortran -lhdf5  -lsz -lz -ldl -lm
HDF5_FLAGS = -lhdf5_fortran -lhdf5 -lz -ldl -lm

# LDFLAGS = $(FFTW3_FLAGS) $(BLAS_FLAGS) $(HDF5_FLAGS)
# LDFLAGS = -I$(MPI_BASE)/Include $(FFTW_MPI_SHLIB) $(BLAS_FLAGS) $(HDF5_FLAGS) #-qmkl=sequential
#LDFLAGS = -L$(LD_LIBRARY_PATH) $(FFTW3_FLAGS) $(BLAS_FLAGS) $(HDF5_FLAGS)
#LDFLAGS = -L$(LD_LIBRARY_PATH) $(FFTW3_FLAGS) $(BLAS_FLAGS) $(HDF5_FLAGS) -mkl=parallel
LDFLAGS = $(FFTW3_LIB) $(FFTW3_FLAGS) $(BLAS_FLAGS) $(HDF5_FLAGS)
#LDFLAGS = $(FFTW3_LIB) $(FFTW3_FLAGS) $(BLAS_LIB) $(BLAS_FLAGS) $(HDF5_FLAGS)
#LDFLAGS = $(FFTW_SHLIB) $(MKL_SHLIB) $(HDF5_F90_SHLIB) $(HDF5_SHLIB) -lz #$(SZIP_LIB) -lz
#LDFLAGS = -lfftw3 -lmkl_core -lmkl_intel_lp64 -mkl=parallel


#===============================================================================
# 2. Make Rules
#===============================================================================
#-------------------------------------------------------------------------------
# Non-module Fortran files to be compiled:
#-------------------------------------------------------------------------------
EXTRA_DIST =    transpose_z_to_x.F90 transpose_x_to_z.F90 transpose_x_to_y.F90\
                transpose_y_to_x.F90 transpose_y_to_z.F90 transpose_z_to_y.F90\
                factor.F90 halo.F90 fft_common.F90 alloc.F90 halo_common.F90

FFILES +=   CalcMaxCFL.F90 CalcLocalDivergence.F90\
            CheckDivergence.F90 CorrectPressure.F90 CorrectVelocity.F90\
            CreateGrid.F90 CreateInitialConditions.F90 CreateResultDirectories.F90\
            DeallocateVariables.F90\
            ExplicitTermsCO2.F90 ExplicitTermsH2O.F90 ExplicitTermsTemp.F90\
            ExplicitTermsVX.F90 ExplicitTermsVY.F90 ExplicitTermsVZ.F90\
            geometry.F90 GlobalQuantities.F90 HdfRoutines.F90\
            IBMRoutines.F90\
            ImplicitAndUpdateCO2.F90 ImplicitAndUpdateH2O.F90 ImplicitAndUpdateTemp.F90\
            ImplicitAndUpdateVX.F90 ImplicitAndUpdateVY.F90 ImplicitAndUpdateVZ.F90\
            InitPressureSolver.F90 InitTimeMarchScheme.F90 InitVariables.F90\
            LocateLargeDivergence.F90 LocateLargeVelocity.F90\
            MovieRoutines.F90 MpiAuxRoutines.F90 obj_io.F90 QuitRoutine.F90\
            ReadFlowField.F90 ReadInputFile.F90 ReadOutlet.F90 SetWallBCs.F90\
            SolveImpEqnUpdate_CO2.F90 SolveImpEqnUpdate_H2O.F90 SolveImpEqnUpdate_Temp.F90\
            SolveImpEqnUpdate_X.F90 SolveImpEqnUpdate_Y.F90 SolveImpEqnUpdate_Z.F90\
            SolvePressureCorrection.F90 StatRoutines.F90 TimeMarcher.F90\
            WriteFlowField.F90 WriteFlowFieldSnapshot.F90 WriteGridInfo.F90 WriteOutlet.F90\
            factorize.F90 main.F90

#-------------------------------------------------------------------------------
# Files that create modules:
#-------------------------------------------------------------------------------
MFILES = param.F90 decomp_2d.F90 AuxiliaryRoutines.F90 decomp_2d_fft.F90 

# Object and module directory:
SRCDIR=src
OBJDIR=obj
OBJS  := $(FFILES:%.F90=$(OBJDIR)/%.o)
MOBJS := $(MFILES:%.F90=$(OBJDIR)/%.o)

FC += -J $(OBJDIR) 

#-------------------------------------------------------------------------------
# Make program 
#-------------------------------------------------------------------------------
PROGRAM = afid 

#Compiling 
all: directories $(PROGRAM) 
$(PROGRAM): $(MOBJS) $(OBJS) 
	$(FC) -o $@ $^ $(LDFLAGS) 

#-------------------------------------------------------------------------------
# Dependencies 
#-------------------------------------------------------------------------------
$(OBJDIR)/param.o: $(SRCDIR)/param.F90
	$(FC) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/AuxiliaryRoutines.o: $(SRCDIR)/AuxiliaryRoutines.F90 
	$(FC) -c -o $@ $< $(LDFLAGS) 
$(OBJDIR)/decomp_2d.o: $(SRCDIR)/decomp_2d.F90
	$(FC) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/decomp_2d_fft.o: $(SRCDIR)/decomp_2d_fft.F90
	$(FC) -c -o $@ $< $(LDFLAGS) 
$(OBJDIR)/ImplicitDecomp.o: $(SRCDIR)/ImplicitDecomp.F90
	$(FC) -c -o $@ $< $(LDFLAGS) 

$(OBJDIR)/%.o: $(SRCDIR)/%.F90 $(MOBJS)
	$(FC) -c -o $@ $< $(LDFLAGS)

#-------------------------------------------------------------------------------
# Clean up 
#-------------------------------------------------------------------------------
clean: 
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*genmod* $(OBJDIR)/*.o obj\
	rm afid

.PHONY: directories
directories: $(OBJDIR) 
$(OBJDIR): 
	mkdir -p ${OBJDIR} 
