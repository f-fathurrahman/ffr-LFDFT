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


