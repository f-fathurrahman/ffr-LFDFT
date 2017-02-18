TOPDIR=../..

ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 t_LF3dp_v3.f90

