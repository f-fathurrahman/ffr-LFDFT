# Status as of 14 August 2017

- Calculation with Lagrange sinc basis is now functional.

- Two additional Poisson solvers, for isolated systems only, are implemented:
  ISF (Goedecker, *et al*) and DAGE (Sundholm, *et al*).

- Atomic coordinate in input file can be given in both in angstrom and bohr.


# Status as of 22 May 2017

- Calculation with GTH nonlocal pseudopotential is now functional.
  However, the total energy does not exactly coincide with result
  from QE, the difference stems from Ewald energy calculation.
  the electronic energy, however, is quite close.

- There are five supported compilers now: ifort, gfortran, g95,
  pgf90, and sunf5.

- There is still some problem with QE's Davidson subroutine for some
  systems.

- Several systems are quite difficult to converge, although it usually
  will converge eventually.

# Status as of 24 April 2017

- QE's Davidson routine is now working with ILU0 preconditioner.
  On some systems, it is failed for yet known reason.

- SCF is working, however, the convergence needs to be stabilized.

- The process of constructing local potentials is now automatic.

# Status as of 17 April 2017

- ILU0 preconditioner from SPARSKIT can be used now.
  The CG algorithm is now faster due to use of preconditioner.
  Selected files from SPARSKIT is added in the repository
  in order to reduce external dependencies.

- Davidson and LOBPCG methods are working with ILU0 preconditioner.
  However, it is not yet working in QE's Davidson subroutine.

# Status as of 22 February 2017

The program is able to solve Kohn-Sham equations using SCF and direct
minimization, albeit from simple test program and simple local potential
(in this case only harmonic potential).

SCF is a bit shaky but it can converge nevertheless.
