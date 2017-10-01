include ../../make.inc

SPARSKIT_SRC = \
formats.f \
ilut.f \
itaux.f \
iters.f \
unary.f \
blassm.f \
matvec.f

OBJ = $(SPARSKIT_SRC:.f=.o) 

#
# Suffix rule for Fortran 90
#
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
	ar rcs libsparskit.a *.o; mv libsparskit.a ../../

# does not delete *.x files
clean:
	rm -rf *.o *.mod libmain.a




