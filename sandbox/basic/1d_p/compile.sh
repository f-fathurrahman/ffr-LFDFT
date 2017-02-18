TOPDIR=../..

#ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 c8_inverse.f90 t_transf.f90 -mkl

ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 t_fun1_p.f90

