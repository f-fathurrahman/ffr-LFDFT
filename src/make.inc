EXE_MAIN = ffr_LFDFT_gfortran.x

MAKE = make
AR = ar

F77 = gfortran -pg -g
F77_OPTS = -O3

F90 = gfortran -pg -g
F90_OPTS = -cpp -Wall -O3 # -fcheck=all # -g -fbacktrace -ffpe-trap=zero,overflow,underflow
F90_OPTS_NOWARN = -cpp -O3

CC = gcc -pg -g
CC_OPTS = -std=c99

LIB_LINALG = -lblas -llapack
LIB_FFTW3 = -lfftw3
LIB_XC = -L/home/efefer/mysoftwares/libxc-3.0.0/lib -lxcf90 -lxc
LIBS = $(LIB_LINALG) $(LIB_FFTW3) $(LIB_XC)

