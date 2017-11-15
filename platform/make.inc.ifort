EXE_MAIN = ffr_LFDFT_ifort.x
MAKE = make
AR = ar
F77 = ifort
F77_OPTS = -O3

F90 = ifort
F90_OPTS = -fpp -xHost -O2 -fpp -warn -nogen-interfaces
F90_OPTS_NOWARN = -fpp -xHost -O2 -fpp

CC_OPTS = -std=c99

#F90_OPTS = -fpp -O3 # -check uninit -xHost
# -CB for check bound of arrays
#F90_OPTS = -CB -warn -nogen-interfaces
#F90_OPTS = -O3 -funroll-all-loops -vector
#LIB_FFTW3 = /home/efefer/mysoftwares/lib/libfftw3.a

LIB_FFTW3 = -lfftw3
LIB_MKL = -L$(MKLROOT)/lib/intel64 -lmkl_core -lmkl_sequential -lmkl_intel_lp64
LIBS = $(LIB_MKL) $(LIB_FFTW3)

