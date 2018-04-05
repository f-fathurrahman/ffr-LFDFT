MODULE m_hamiltonian

  IMPLICIT NONE

  REAL(8), ALLOCATABLE :: V_ps_loc(:)
  REAL(8), ALLOCATABLE :: V_ps_loc_long(:)  ! long part of local pseudopotential
  REAL(8), ALLOCATABLE :: V_Hartree(:)
  REAL(8), ALLOCATABLE :: V_xc(:)

  REAL(8), ALLOCATABLE :: Rhoe(:)

  REAL(8), ALLOCATABLE :: betaNL_psi(:,:) ! NbetaNL, Nstates 

END MODULE 

