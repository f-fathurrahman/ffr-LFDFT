#include platform/make.inc.ifort
include platform/make.inc.gfortran
#include platform/make.inc.g95
#include platform/make.inc.pgi
#include platform/make.inc.sun

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
m_input_vars.f90 \
welcome.f90 \
goodbye.f90 \
timestamp.f90 \
compile_info.f90 \
get_month_name.f90 \
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
init_V_ps_loc.f90 \
init_V_ps_loc_G_interp.f90 \
init_V_ps_loc_gaussian_G.f90 \
init_V_ps_loc_gaussian.f90 \
init_V_coul_G_interp.f90 \
init_strfact_shifted.f90 \
Ylm_real.f90 \
xsf.f90 \
calc_dr_periodic_1pnt.f90 \
eval_LF1d_c.f90 \
eval_LF1d_p.f90 \
eval_LF1d_sinc.f90 \
KS_solve_SCF.f90 \
calc_evals.f90 \
read_control.f90 \
read_system.f90 \
read_electrons.f90 \
read_atomic_species.f90 \
read_atomic_positions.f90 \
read_input.f90 \
setup_atoms.f90 \
setup_PsPot.f90 \
setup_LF3d.f90 \
setup_from_input.f90 \
setup_options.f90 \
interp_LF1d_p.f90 \
interp_LF3d_p.f90 \
calc_Ewald_qe.f90 \
rgen.f90 \
hpsort.f90 \
atmlength.f90 \
atom_znucl.f90 \
gen_guess_rho_gaussian.f90 \
gen_random_evecs.f90 \
gen_gaussian_evecs.f90 \
calc_occupations.f90 \
calc_entropy.f90 \
fermi_dirac.f90 \
normalize_rhoe.f90 \
Poisson_solve_ISF.f90 \
m_Faddeeva.f90 \
Poisson_solve_DAGE.f90


SPARSKIT_SRC = \
formats.f \
ilut.f \
itaux.f \
iters.f \
unary.f \
blassm.f \
matvec.f

C_SRC = \
Faddeeva.c

#pconv.f90 
#pfft3d.f90
#POISSON_ISF_SRC = \
Build_Kernel.f90 \
fft3d.f90 \
PSolver_Kernel.f90 \
scaling_function.f90 \
smooth.f90 \
gequad.f


OBJ = $(SRC:.f90=.o) $(SRC:.f=.o) $(SPARSKIT_SRC:.f=.o) \
$(C_SRC:.c=.o) \
$(POISSON_ISF_SRC:.f=.o) $(POISSON_ISF_SRC:.f90=.o)

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
	ar rcs libmain.a *.o

# Targets
main: lib ffr_LFDFT.f90
	$(F90) $(F90_OPTS) ffr_LFDFT.f90 -o $(EXE_MAIN) libmain.a $(LIBS)

# does not delete *.x files
clean:
	rm -rf *.o *.mod libmain.a

# also delete *.x files
cleanx:
	rm -rf *.o *.mod libmain.a *.x



