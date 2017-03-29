MODULE m_hamiltonian

  IMPLICIT NONE

  REAL(8), ALLOCATABLE :: V_ps_loc(:)
  REAL(8), ALLOCATABLE :: V_Hartree(:)
  REAL(8), ALLOCATABLE :: V_xc(:)

  REAL(8), ALLOCATABLE :: Rhoe(:)

  ! Diagonal elements of Hamiltonian
  REAL(8), ALLOCATABLE :: K_diag(:)

END MODULE 

