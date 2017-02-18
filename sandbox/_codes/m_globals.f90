
! This module contains several global variables and parameter used in this
! program
MODULE m_globals
  USE m_LF3d
  IMPLICIT NONE
  ! These parameters are similar for x, y, and z directions
  INTEGER :: N
  REAL(8) :: A, B
  REAL(8) :: h
  !
  TYPE(LF3d_t) :: LF
  !
  REAL(8), ALLOCATABLE :: Vpot(:)
  REAL(8), ALLOCATABLE :: Rho(:)
  REAL(8), ALLOCATABLE :: Vhartree(:)
  REAL(8) :: Ehartree
  !
  REAL(8), ALLOCATABLE :: h_diag(:)
  !
  REAL(8), ALLOCATABLE :: evals(:), evecs(:,:)
  INTEGER :: Nstate
  REAL(8) :: Nelec
  REAL(8), ALLOCATABLE :: Focc(:) ! occupation number
  !
  CHARACTER(16) :: LF_type
  CHARACTER(16) :: Solution_Method
  !
  REAL(8), ALLOCATABLE :: evalsTx(:), evalsTy(:), evalsTz(:), Kprec(:)
  !
  REAL(8), ALLOCATABLE :: eVtau(:), eVtau_old(:)
END MODULE


