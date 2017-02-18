TOPDIR=../..
LIBHDF5="/home/efefer/mysoftwares/hdf5-1.8.16_intel/lib/libhdf5_fortran.a /home/efefer/mysoftwares/hdf5-1.8.16_intel/lib/libhdf5.a -lz"

ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 t_applyInverseKin.f90 -mkl

#ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 t_inverseKin.f90 eig_dsyev.f90 -mkl

#ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 t_inverseKin.f90 eig_dsyev.f90 $TOPDIR/rw_hdf5.o -mkl $LIBHDF5

#ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 sort.f90 t_kin1.f90 -mkl

