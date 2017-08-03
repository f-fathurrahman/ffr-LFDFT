#include ../../platform/make.inc.ifort
include ../../platform/make.inc.gfortran
#include ../../platform/make.inc.g95
#include ../../platform/make.inc.pgi
#include ../../platform/make.inc.sun

SRC = \
Build_Kernel.f90 \
fft3d.f90 \
pconv.f90 \
pfft3d.f90 \
PSolver_Kernel.f90 \
scaling_function.f90 \
smooth.f90 \
gequad.f

OBJ = $(SRC:.f90=.o) $(SRC:.f=.o) 

#
# Suffix rule for Fortran 90
#
%.mod :
	@if [! -f $@ ]; then \
		rm $(*F).o; \
		fi
	$(MAKE) $<

%.o : %.f90
	$(F90) $(F90_OPTS) -c -o $(*F).o $<

#
# Fortran 77 sources
# supress warning
.SUFFIXES: .o .f
.f.o:
	$(F77) -c $(F77_OPTS) $<

# Targets
lib: $(OBJ)
	ar rcs libPoissonISF.a *.o

# Targets
main: lib MAIN.f90
	$(F90) -I../../ $(F90_OPTS) MAIN.f90 init_density.f90 -o MAIN.x libPoissonISF.a ../../libmain.a $(LIBS)

# does not delete *.x files
clean:
	rm -rf *.o *.mod libmain.a

# also delete *.x files
cleanx:
	rm -rf *.o *.mod libmain.a *.x



