SUBROUTINE fermi_dirac( is_spinpol, Nstates, evals, efermi, swidth, Focc )
  IMPLICIT NONE 
  LOGICAL :: is_spinpol
  INTEGER :: Nstates
  REAL(8) :: efermi, swidth
  REAL(8) :: evals(Nstates)
  REAL(8) :: Focc(Nstates)  ! output
  INTEGER :: ist
  !
  IF( is_spinpol ) THEN 
    DO ist = 1,Nstates
      Focc(ist) = 1.d0/( 1.d0 + exp( ( evals(ist) - efermi )/(0.5d0*swidth) ) )
    ENDDO 
  ELSE 
    DO ist = 1,Nstates
      Focc(ist) = 2.d0/( 2.d0 + exp( ( evals(ist) - efermi )/(0.5d0*swidth) ) )
    ENDDO 
  ENDIF 

END SUBROUTINE 

