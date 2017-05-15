#include platform/make.inc.ifort
#include platform/make.inc.gfortran
#include platform/make.inc.g95
include platform/make.inc.pgi

SRC = \
m_constants.f90 \
m_options.f90 \
m_atoms.f90 \
m_LF3d.f90 \
m_nabla2_sparse.f90 \
m_hamiltonian.f90 \
m_states.f90 \
m_energies.f90 \
m_Ps_HGH.f90 \
m_PsPot.f90 \
init_PsPot.f90 \
dealloc_PsPot.f90 \
fft_fftw3.f90 \
LDA_VWN.f90 \
init_atoms_xyz.f90 \
init_grid_1d_p.f90 \
init_grid_1d_c.f90 \
init_grid_1d_sinc.f90 \
init_deriv_matrix_p.f90 \
init_deriv_matrix_c.f90 \
init_deriv_matrix_sinc.f90 \
init_LF3d_p.f90 \
init_LF3d_c.f90 \
init_LF3d_sinc.f90 \
init_states.f90 \
info_LF3d.f90 \
info_PsPot.f90 \
info_energies.f90 \
info_atoms.f90 \
mm_to_nn.f90 \
init_gvec.f90 \
init_V_ps_loc_G.f90 \
init_V_coul_G.f90 \
dealloc_LF3d.f90 \
op_nabla2.f90 \
init_nabla2_sparse.f90 \
dealloc_nabla2_sparse.f90 \
dealloc_atoms.f90 \
Poisson_solve_cg.f90 \
Poisson_solve_pcg.f90 \
Poisson_solve_fft.f90 \
Poisson_solve_fft_MT.f90 \
op_H.f90 \
op_V_ps_NL.f90 \
init_betaNL.f90 \
calc_betaNL_psi.f90 \
diag_davidson_qe.f90 \
rdiaghg.f90 \
alloc_hamiltonian.f90 \
dealloc_hamiltonian.f90 \
init_V_ps_loc_harmonic.f90 \
ortho_gram_schmidt.f90 \
orthonormalize.f90 \
ortho_check.f90 \
Sch_solve_diag.f90 \
calc_energies.f90 \
calc_Rhoe.f90 \
update_potentials.f90 \
KS_solve_Emin_cg.f90 \
calc_grad.f90 \
calc_dr_periodic.f90 \
calc_dr.f90 \
init_strfact.f90 \
KS_solve_Emin_pcg.f90 \
op_K.f90 \
m_ilu0_prec.f90 \
init_ilu0_prec.f90 \
prec_ilu0.f90 \
dealloc_ilu0_prec.f90 \
diag_davidson.f90 \
diag_lobpcg.f90 \
calc_Ewald.f90 \
mixadapt.f90 \
mixbroyden.f90 \
mixerifc.f90 \
mixlinear.f90 \
shift_atoms.f90 \
bspline.f90 \
init_V_ps_loc_G_interp.f90 \
init_V_coul_G_interp.f90 \
init_strfact_shifted.f90 \
Ylm_real.f90 \
xsf.f90 \
calc_dr_periodic_1pnt.f90 \
eval_LF1d_c.f90 \
eval_LF1d_p.f90 \
eval_LF1d_sinc.f90 \
KS_solve_SCF.f90

SPARSKIT_SRC = \
formats.f \
ilut.f \
itaux.f \
iters.f \
unary.f \
blassm.f \
matvec.f


OBJ = $(SRC:.f90=.o) $(SRC:.f=.o) $(SPARSKIT_SRC:.f=.o)

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
	$(F90) -O3 -c $<


# Targets
lib: $(OBJ)
	ar rcs libmain.a *.o

clean:
	rm -rf *.o *.mod libmain.a *.x


