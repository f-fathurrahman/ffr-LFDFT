EXE_MAIN = ffr_LFDFT_pgi_llvm.x
MAKE = make
AR = ar

F77 = pgfortran-llvm
F77_OPTS = -O3

F90 = pgf90-llvm
F90_OPTS = -Minform=warn -fast -Mpreprocess
F90_OPTS_NOWARN = -fast -Mpreprocess

CC = pgcc
CC_OPTS = -fast -c99 -silent

LIB_LINALG = -lblas -llapack
LIB_FFTW3 = -lfftw3
LIB_XC = -L/home/efefer/mysoftwares/libxc-3.0.0/lib -lxcf90 -lxc
LIBS = $(LIB_LINALG) $(LIB_FFTW3) $(LIB_XC)


