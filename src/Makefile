#FC := mpifort -fopenmp -g -fcheck=array-temps,bounds,do,mem,pointer,recursion
#FC := ftn -fopenmp -O2 -march=haswell -ffast-math -funroll-loops
#FC := mpifort -g -fbounds-check
FC := mpifort -fopenmp -O3 -ffast-math -funroll-loops -mtune=native

OBJDIR := obj

$(OBJDIR)/%.o: %.f90
	$(FC) -I$(OBJDIR) -J$(OBJDIR) -c $< -o $@ -fPIC

%.o: %.mod

expm.o: $(addprefix $(OBJDIR)/, iso_precisions.o sparse.o one_norms.o operators.o) expm.f90

one_norms.o: $(addprefix $(OBJDIR)/, iso_precisions.o sparse.o operators.o) one_norms.f90

operators.o: $(addprefix, $(OBJDIR)/, iso_precisions.o sparse.o) operators.f90

sparse.o: $(addprefix $(OBJDIR)/, iso_precisions.o) sparse.f90

iso_precisions.o: iso_precisions.f90

clean:
	rm obj/*.mod obj/*.o

