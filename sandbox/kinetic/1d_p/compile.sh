TOPDIR=../..
ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 t_kin1.f90 -mkl

