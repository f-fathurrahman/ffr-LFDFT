! This module contains several global variables and parameter used in this
! program
MODULE m_globals
  USE m_LF3d
  IMPLICIT NONE
  ! These parameters are similar for x, y, and z directions
  INTEGER :: N
  REAL(8) :: A, B
  REAL(8) :: deltaV
  INTEGER :: Nbasis
  !
  TYPE(LF3d_t) :: LF
  !
  REAL(8), ALLOCATABLE :: Vpot(:)
  REAL(8), ALLOCATABLE :: Vxc(:)
  REAL(8), ALLOCATABLE :: Rho(:), Vhartree(:)
  REAL(8), ALLOCATABLE :: evecs(:,:), evals(:)
  !
  REAL(8) :: Etot
  REAL(8) :: Ekin, Epot, Ehartree, Exc
  !
  INTEGER :: Nstates
  REAL(8), ALLOCATABLE :: Focc(:)
END MODULE

