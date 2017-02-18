TOPDIR=../..

ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 t_fun1_sinc.f90

#gfortran -Wall ../m_constants.f90 ../m_LF1d.f90 t_fun1_sinc.f90

