# Solving Schrodinger equation

This directory contains examples for solving Schrodinger equations using PCG minimization
of band energies.

The related subroutines are: `Sch_solve_Emin_pcg` and `calc_Ebands`.

For systems containing only local potentials, the global variable
`NbetaNL` defined in module `m_PsPot` need to be set to zero.


