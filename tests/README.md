This diretory contains various quick and dirty tests for subroutines in
ffr-LFDFT.

Currently, these tests must be compiled and linked manually. Most of
them are quite simple, and only require linking to libmain.a
Several tests however require longer command
to compile and link. The script build.sh can be used to build
these tests.

TODO:

- Explanation about each test
- Instruction to compile the tests
- Reorganization

List of stand-alone programs:

- `Emin_cg_harm`: Kohn-Sham energy minimization of periodic system with
   harmonic potential

- `do_SCF`: Kohn-Sham energy via self-consistent field

- `do_Emin_pcg`: Kohn-Sham energy minimization




