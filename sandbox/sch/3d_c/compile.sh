TOPDIR=../../

#ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 t_LF3d_c_v5.f90 -mkl

#ifort -fpp -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 regterg1.f90 rdiaghg.f90 t_LF3d_c_v4.f90 -mkl

#ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 sort.f90 t_LF3d_c_v3.f90 -mkl

ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 t_LF3d_c_v2.f90 -mkl


#ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 $TOPDIR/rw_hdf5.o t_LF3d_c_v1.f90 -mkl ~/mysoftwares/hdf5-1.8.16_intel/lib/libhdf5_fortran.a ~/mysoftwares/hdf5-1.8.16_intel/lib/libhdf5.a -lz

#gfortran -Wall $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 t_LF3d_c_v1.f90 -lblas -llapack
