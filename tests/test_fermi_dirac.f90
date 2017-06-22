PROGRAM test_calc_occupations

  USE m_states
  IMPLICIT NONE 
  INTEGER :: ist
  REAL(8) :: Tbeta, efermi
  REAL(8), ALLOCATABLE :: evals(:)

  Nelectrons = 3.d0
  Nstates = 8

  Nstates_occ = Nstates - Nelectrons
  WRITE(*,*) 'Nstates_occ = ', Nstates_occ

  ALLOCATE( Focc(Nstates) )
  Focc(1) = 1.d0
  Focc(2) = 1.d0
  Focc(3) = 1.d0
  Focc(4) = 0.d0
  Focc(5) = 0.d0
  Focc(6) = 0.d0
  Focc(7) = 0.d0
  Focc(8) = 0.d0

  ALLOCATE( evals(Nstates) )
  evals(1) = -5.3d0
  evals(2) = -3.1d0
  evals(3) = -2.02d0
  evals(4) = -1.01d0
  evals(5) = -1.d0
  evals(6) = -0.8d0
  evals(7) = -0.5d0
  evals(8) = -0.1d0

  efermi = 0.5d0*( evals(4) - evals(3) )
  Tbeta = 0.001d0

  WRITE(*,'(1x,A,2F18.10)') 'Fermi, Tbeta = ', efermi, Tbeta
  CALL fermi_dirac( Nstates, evals, efermi, Tbeta, Focc )

  DO ist = 1,Nstates
    WRITE(*,'(1x,I4,2F18.10)') ist, evals(ist), Focc(ist)
  ENDDO 
  WRITE(*,*) 'sum(Focc) = ', sum(Focc)

END PROGRAM 

