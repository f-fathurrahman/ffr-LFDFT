include ../../make.inc

POISSON_ISF_SRC = \
Build_Kernel.f90 \
fft3d.f90 \
PSolver_Kernel.f90 \
scaling_function.f90 \
smooth.f90 \
gequad.f


OBJ = $(POISSON_ISF_SRC:.f=.o) $(POISSON_ISF_SRC:.f90=.o)

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

#
# C source
#
.SUFFIXES: .o .c
.c.o:
	$(CC) -c $(CC_OPTS) $<

# Targets
lib: $(OBJ)
	ar rcs libpoisson_ISF.a *.o; mv libpoisson_ISF.a ../../

# does not delete *.x files
clean:
	rm -rf *.o *.mod libpoisson_ISF.a

# special rules to suppress warnings

Build_Kernel.o: Build_Kernel.f90
	$(F90) $(F90_OPTS_NOWARN) -c -o Build_Kernel.o Build_Kernel.f90

fft3d.o: fft3d.f90
	$(F90) $(F90_OPTS_NOWARN) -c -o fft3d.o fft3d.f90

PSolver_Kernel.o: PSolver_Kernel.f90
	$(F90) $(F90_OPTS_NOWARN) -c -o PSolver_Kernel.o PSolver_Kernel.f90

scaling_function.o: scaling_function.f90
	$(F90) $(F90_OPTS_NOWARN) -c -o scaling_function.o scaling_function.f90

smooth.o: smooth.f90
	$(F90) $(F90_OPTS_NOWARN) -c -o smooth.o smooth.f90



