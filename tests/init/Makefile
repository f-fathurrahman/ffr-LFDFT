TOPDIR = ../../

#include ../../platform/make.inc.ifort
#include ../../platform/make.inc.gfortran
include $(TOPDIR)/platform/make.inc.g95

# Targets
ex_init:
	$(F90) -I$(TOPDIR) $(F90_OPTS) ex_init.f90 -o ex_init.x $(TOPDIR)/libmain.a

clean:
	rm -rv *.x


