MODULE m_states

  IMPLICIT NONE 
  
  INTEGER :: Nstates
  INTEGER :: Nstates_occ, Nstates_extra
  INTEGER :: Nspin
  REAL(8) :: Nelectrons

  REAL(8), ALLOCATABLE :: KS_evals(:)
  REAL(8), ALLOCATABLE :: KS_evecs(:,:)

  REAL(8), ALLOCATABLE :: Focc(:)

END MODULE 

