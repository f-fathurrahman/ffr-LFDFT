#include platform/make.inc.ia32.ifort.qed
#include platform/make.inc.ia32.gfortran
include platform/make.inc.ifort
#include platform/make.inc.gfortran

SRC = m_constants.f90 m_LF1d.f90 m_LF3d.f90 m_globals.f90 \
  init_potential.f90 apply_hamiltonian.f90 davidson_qe.f90 read_arguments.f90 \
  utils.f90 g_psi.f90 solve_diagonalize.f90 solve_minim.f90 get_Etot.f90 \
  write_grad.f90 precond_invKin.f90 get_grad.f90 minim_cg_v1.f90 \
  minim_pcg_v1.f90 minim_pcg_v2.f90 solve_MG.f90 avgPrec.f90 calc_rho.f90 \
  solve_poisson.f90

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

main:
	$(F90) $(F90_OPTS) main.f90 libmain.a -mkl -o sch.x

clean:
	rm -rf *.o *.mod libmain.a

