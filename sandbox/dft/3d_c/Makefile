#include platform/make.inc.ia32.ifort.qed
#include platform/make.inc.ia32.gfortran
include platform/make.inc.ifort
#include platform/make.inc.gfortran

SRC = m_constants.f90 m_LF1d.f90 m_LF3d.f90 \
      m_globals.f90 ortho_gram_schmidt.f90 eig_dsyev.f90 \
      rand_vec.f90 rdiaghg.f90 regterg1.f90 init_pot.f90 \
      apply_Laplacian.f90 apply_Ham.f90 get_Etot.f90 get_grad.f90 \
      minimE_sd.f90 minimE_linmin.f90 minimE_cg.f90 mixlinear.f90 mixadapt.f90 \
      diag_Ham.f90 solve_poisson.f90 get_Rho.f90 print_vars.f90 \
      LDA_VWN.f90 libxc_LDA.f90 \
      rcgdiagg.f90 mixpulay.f90

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
	$(F90) $(F90_OPTS) $(INC_LIBXC) -c -o $(*F).o $<

# not working?
#libxc_LDA.o:
#	$(F90) $(F90_OPTS) $(INC_LIBXC) -c libxc_LDA.f90

#
# Fortran 77 sources
#
.SUFFIXES: .o .f
.f.o:
	$(F77) $(F77_OPTS) -c $<


# Targets
lib: $(OBJ)
	ar rcs libmain.a *.o

main_harm: lib
	$(F90) $(F90_OPTS) main_harm.f90 libmain.a $(LIBS) -o main_harm.x

main_harm_minE: lib
	$(F90) $(F90_OPTS) main_harm_minE.f90 libmain.a $(LIBS) -o main_harm_minE.x

main_Hpseudo: lib
	$(F90) $(F90_OPTS) main_Hpseudo.f90 libmain.a $(LIBS) -o main_Hpseudo.x

clean:
	rm -rf *.o *.mod libmain.a *.x


