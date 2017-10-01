EXE_MAIN = ffr_LFDFT_pgi.x
MAKE = make
AR = ar

F77 = pgf77
F77_OPTS = -O3

F90 = pgf90
F90_OPTS = -Minform=inform -fast -Mpreprocess
F90_OPTS_NOWARN = -fast -Mpreprocess

LIB_LINALG = -lblas -llapack
LIB_FFTW3 = -lfftw3
LIBS = $(LIB_LINALG) $(LIB_FFTW3)
