
# Specify the file that has the input PIPs/FIs:
PIPfile = obj/BrClH2.f90
dPIPfile = obj/BrClH2-derivatives.f90

# On an example Linux HPC server, these LIBs work:
# LIB = -L/opt/apps/software/numlib/ScaLAPACK/2.1.0-gompi-2021b-fb/lib -lscalapack -lmkl_intel_thread -lmkl_core -liomp5 -lmkl_intel_lp64

# On an example Linux desktop, these LIBs and INLCUDEs work:
MKLROOT = /opt/intel/oneapi/mkl/2023.1.0
LIB = -L${MKLROOT}/lib/intel64 -qmkl -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64
INCLUDE = -I${MKLROOT}/include/intel64/ilp64 -i8  -I"${MKLROOT}/include"

################################################################################

EXE = training.x lookAtBounds.x testEnergyAndOrForces.x

FC = ifort
FFLAGS = -O2 -qopenmp -module mod -traceback
FORTS = obj/pes.f90 $(PIPfile) $(dPIPfile)

OBJ = $(patsubst %.f90, %.o, $(FORTS))

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) $(LIB) $(INCLUDE) \
	$+ $(patsubst %.x, src/%.f90, $@) -o $@

obj/%.o : obj/%.f90
	$(FC) -c $< $(FFLAGS) $(LIBFFTW) $(LIB) -o $@

all: $(EXE)

clean:
	rm -f  mod/*.mod $(EXE) obj/*.o
