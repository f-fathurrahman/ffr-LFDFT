#include platform/make.inc.ifort
include platform/make.inc.gfortran
#include platform/make.inc.g95

SRC = \
m_constants.f90 \
fft_fftw3.f90 \
m_LF3d.f90 \
init_grid_1d_p.f90 \
init_deriv_matrix_p.f90 \
init_LF3d_p.f90 \
info_LF3d.f90 \
mm_to_nn.f90 \
init_gvec.f90 \
dealloc_LF3d.f90 \
apply_Laplacian.f90 \
solve_poisson_cg.f90 \
solve_poisson_fft.f90

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


