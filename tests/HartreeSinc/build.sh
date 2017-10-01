#!/bin/bash
#g95 -Wall -O -g -fno-second-underscore -ftrace=full \
#ifort -warn -nogen-interfaces -O -g -CB \
gfortran -O3 -Wall -g -fbacktrace \
-I../../ HartreeSinc_main.f90 ../../libmain.a \
HartreeSinc_t_sampling.f90 gauleg.f90 init_density.f90 \
init_density_unnormalized.f90 \
compute_F.f90 construct_F.f90 Cwrap_Faddeeva.o compute_potential.f90 \
-lstdc++ \
-o HartreeSinc_main.x
