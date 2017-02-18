#include platform/make.inc.ia32.ifort.qed
#include platform/make.inc.ia32.gfortran
include platform/make.inc.ifort
#include platform/make.inc.gfortran

SRC = \
m_constants.f90 \
m_LF3d.f90 \
init_grid_1d_p.f90 \
init_deriv_matrix_p.f90

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
#
.SUFFIXES: .o .f
.f.o:
	$(F77) $(F77_OPTS) -c $<


# Targets
lib: $(OBJ)
	ar rcs libmain.a *.o

clean:
	rm -rf *.o *.mod libmain.a *.x


