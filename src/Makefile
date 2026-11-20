#Makefile

# Shell command for the python interpreter.
PYTHON_CALL ?= python3
# Set to TRUE to compile with debug flags.
DEBUG ?= FALSE

# Fortran + MPI compiler.
ifeq ($(PAWSEY_CLUSTER), magnus)
	FC = ftn
else
	FC = mpifort
endif

# Fortran compiler flags.
ifeq ($(DEBUG),TRUE)
	FCFLAGS = -O0 -g -fcheck=array-temps,bounds,do,mem,pointer,recursion -Wall -Wextra -Wconversion -pedantic
else
	FCFLAGS = -O3 -ffast-math -funroll-loops -mtune=native
endif

ifeq ($(OPENMP),TRUE)
	FCFLAGS := $(FCFLAGS) -fopenmp
endif

EXT_SUFFIX := $(shell $(PYTHON_CALL) -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')

OBJDIR := obj
PYDIR := ../qsw_mpi

$(shell mkdir -p $(OBJDIR))

H5PFC ?= h5pfc
ifeq ($(shell uname -s),Darwin)
# Keep the parallel HDF5 symbols private when h5py loads its own HDF5 library.
HDF5_ENV = HDF5_USE_SHLIB=no LDFLAGS="$(LDFLAGS) -Wl,-exported_symbol,_PyInit_fHDF5"
endif

.PHONY: all parallel_io
all: $(PYDIR)/fMPI$(EXT_SUFFIX)
ifneq ($(shell command -v $(H5PFC)),)
all: parallel_io
endif

parallel_io: $(PYDIR)/fHDF5$(EXT_SUFFIX)

$(PYDIR)/fHDF5$(EXT_SUFFIX): fHDF5.pyf parallel_io.f90
	$(HDF5_ENV) $(PYTHON_CALL) -m numpy.f2py --f90exec="$(H5PFC)" --f77exec="$(H5PFC)" --f90flags="$(FCFLAGS)" -c fHDF5.pyf parallel_io.f90
	mv fHDF5$(EXT_SUFFIX) $(PYDIR)

$(OBJDIR)/%.o: %.f90
	$(FC) $(FCFLAGS) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@ -fPIC

%.o: %.mod

# fMPI f2py shared object library.
$(PYDIR)/fMPI$(EXT_SUFFIX):  $(addprefix $(OBJDIR)/, iso_precisions.o sparse.o operators.o one_norms.o expm.o) $(PYDIR)/foperators$(EXT_SUFFIX) fMPI.f90 fMPI.pyf
	$(PYTHON_CALL) -m numpy.f2py --f90exec="$(FC)" --f90flags="$(FCFLAGS)" -lgomp -c fMPI.pyf iso_precisions.f90 sparse.f90 operators.f90 one_norms.f90 expm.f90 -m fMPI fMPI.f90 
	mv fMPI$(EXT_SUFFIX) $(PYDIR)


# foperators f2py shared object library.
$(PYDIR)/foperators$(EXT_SUFFIX): $(addprefix $(OBJDIR)/, iso_precisions.o sparse.o operators.o) foperators.f90
	$(PYTHON_CALL) -m numpy.f2py --f90exec="$(FC)" --f90flags="-I$(OBJDIR) $(FCFLAGS)" -lgomp -c foperators.pyf iso_precisions.f90 sparse.f90 operators.f90 -m foperators foperators.f90
	mv foperators$(EXT_SUFFIX) $(PYDIR)

# Compile fortran subroutines.
$(OBJDIR)/expm.o: $(addprefix $(OBJDIR)/, iso_precisions.o sparse.o one_norms.o) expm.f90
$(OBJDIR)/one_norms.o: $(addprefix $(OBJDIR)/, iso_precisions.o sparse.o) one_norms.f90
$(OBJDIR)/operators.o: $(addprefix, $(OBJDIR)/, iso_precisions.o sparse.o) operators.f90
$(OBJDIR)/sparse.o: $(addprefix $(OBJDIR)/, iso_precisions.o) sparse.f90
$(OBJDIR)/iso_precisions.o: iso_precisions.f90

####

clean:
	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod

veryclean:
	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod $(PYDIR)/*.so


