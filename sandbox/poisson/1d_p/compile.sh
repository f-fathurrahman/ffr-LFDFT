TOPDIR=../..

ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 t_poisson_1d_p.f90 -mkl
