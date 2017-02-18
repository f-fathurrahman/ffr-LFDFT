TOPDIR=../..


ifort -I /home/efefer/intel/mkl/include -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 t_poisson_pardiso.f90 solve_poisson_pardiso.f90 -mkl

#ifort -warn -nogen-interfaces $TOPDIR/m_constants.f90 $TOPDIR/m_LF1d.f90 $TOPDIR/m_LF3d.f90 nabla2_sparse_v1.f90
