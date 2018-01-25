# Solving Schrodinger equation

This directory contains examples for solving Schrodinger equations using various
methods:

- Davidson diagonalization. There are two subroutines: `diag_davidson`
  and `diag_davidson_qe` which is adapted from PWSCF. I found out that
  `diag_davidson_qe` is more unstable than other diagonalization subroutines
  that are available in `ffr-LFDFT`. There may be subtle bugs I'm missing
  when integrating this into `ffr-LFDFT`.
  
- LOBPCG diagonalization: `diag_lobpcg`.

- Minimization of of band energies: `Sch_solve_Emin_pcg`.

Note the different normalization convention used in these subroutines.

For systems containing only local potentials, the global variable
`NbetaNL` defined in module `m_PsPot` need to be set to zero.


