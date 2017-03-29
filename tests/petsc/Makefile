# For the moment, gfortran will used throughout because this is the default
# Fortran compiler used by Open MPI

INC = -I/home/efefer/mysoftwares/petsc-3.7.5/include/ -I../../

MAKE = make
AR = ar
F90 = mpifort
F90_OPTS = -Wall -O3 $(INC) -ffree-form
LIB_LINALG = -lblas -llapack
LIB_FFTW3 = -lfftw3
LIBS = $(LIB_LINALG) $(LIB_FFTW3)

SRC = \
add_nabla2_x.F \
add_nabla2_y.F \
add_nabla2_z.F

OBJ = $(SRC:.F=.o)

#
# Suffix rule for Fortran 90
#
%.mod :
	@if [! -f $@ ]; then \
		rm $(*F).o; \
		fi
	$(MAKE) $<

%.o : %.F
	$(F90) $(F90_OPTS) -c -o $(*F).o $<


# Targets
lib: $(OBJ)
	ar rcs libLFDFT_petsc.a *.o

clean:
	rm -rf *.o *.mod libmain.a *.x


