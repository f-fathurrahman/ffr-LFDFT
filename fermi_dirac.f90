SUBROUTINE fermi_dirac( Nstates, evals, efermi, Tbeta, Focc )
  IMPLICIT NONE 
  INTEGER :: Nstates
  REAL(8) :: efermi, Tbeta
  REAL(8) :: evals(Nstates)
  REAL(8) :: Focc(Nstates)
  INTEGER :: ist
  !
  DO ist = 1,Nstates
    Focc(ist) = 1.d0/( 1.d0 + exp( Tbeta*( evals(ist) - efermi ) ) )
  ENDDO 

END SUBROUTINE 

