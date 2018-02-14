PROGRAM test_calc_occupations

  USE m_states
  IMPLICIT NONE 
  INTEGER :: ist
  REAL(8) :: swidth, efermi
  INTEGER :: i_smear
  REAL(8), ALLOCATABLE :: evals(:)
  REAL(8) :: Nelectrons_calc
  INTEGER :: Nstates_Nocc
  REAL(8) :: Entropy

  Nelectrons = 3.d0
  Nstates = 8

  ! Assume singly-occupied states
  Nstates_Nocc = Nstates - int(Nelectrons)
  Nstates_occ  = int(Nelectrons)

  WRITE(*,*) 'Nstates_occ = ', Nstates_occ

  ALLOCATE( Focc(Nstates) )
  ALLOCATE( evals(Nstates) )

  ! Set eigenvalues manually
  evals(1) = -5.3d0
  evals(2) = -3.1d0
  evals(3) = -2.02d0
  evals(4) = -1.01d0
  evals(5) = -1.d0
  evals(6) = -0.8d0
  evals(7) = -0.5d0
  evals(8) = -0.1d0

  efermi = 0.5d0*( evals(4) + evals(3) )
  WRITE(*,*)
  WRITE(*,*) 'Efermi is set to be ', efermi

  DO i_smear = 1, 10
    swidth = i_smear*0.01d0  ! in Ha
    WRITE(*,*)
    WRITE(*,'(1x,A,F6.3)') 'Smearing width (Ha): ', swidth
    !
    CALL fermi_dirac( .FALSE., Nstates, evals, efermi, swidth, Focc )
    !
    DO ist = 1,Nstates
      WRITE(*,'(1x,I4,2F18.10)') ist, evals(ist), Focc(ist)
    ENDDO 
    WRITE(*,'(1x,A,F18.10)') 'sum(Focc) = ', sum(Focc)
    !
    Nelectrons_calc = 0.d0
    DO ist = 1,Nstates
      IF( evals(ist) <= efermi ) THEN 
        Nelectrons_calc = Nelectrons_calc + Focc(ist)
      ENDIF 
    ENDDO 
    WRITE(*,'(1x,A,F18.10)') 'Nelectrons_calc = ', Nelectrons_calc
    !
    CALL calc_Entropy( .FALSE., Nstates, swidth, Focc, Entropy )
    !
    WRITE(*,'(1x,A,F18.10)') 'Entropy = ', Entropy
  ENDDO 



END PROGRAM 

