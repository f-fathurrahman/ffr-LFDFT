TOPDIR=../..

ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 t_LF3d_c_poisson4.f90 -o p4.x

#ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 t_LF3d_c_poisson3.f90 -o p3.x

#ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 t_LF3d_c_poisson2.f90 -mkl -o p2.x

#ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 t_LF3d_c_poisson1.f90 -o p1.x
