TOPDIR=../..


ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 t_sch3d_v2.f90 -mkl

#ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 t_sch3d_v1.f90 -mkl

