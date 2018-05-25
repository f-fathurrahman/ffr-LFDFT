! Calculate electronic entropy
SUBROUTINE calc_Entropy( is_spinpol, Nstates, swidth, Focc, Entropy )
  IMPLICIT NONE 
  LOGICAL :: is_spinpol
  INTEGER :: Nstates
  REAL(8) :: swidth
  REAL(8) :: Focc(Nstates)
  REAL(8) :: Entropy
  INTEGER :: ist
  REAL(8), PARAMETER :: SMALL_OCC = 1.e-10
  REAL(8) :: e

  !
  ! FIXME: No check on value of Focc
  !

  e = 0.d0
  DO ist = 1, Nstates
    IF( Focc(ist) > SMALL_OCC ) THEN 
      e = e + Focc(ist)*log(Focc(ist))
    ELSE
      IF( is_spinpol ) THEN 
        e = e + (1.d0 - Focc(ist))*log(1.d0 - Focc(ist))
      ELSE
        ! spin-degenerate case
        e = e + 2.d0*(1.d0 - 0.5d0*Focc(ist))*log(1.d0 - 0.5d0*Focc(ist))
      ENDIF 
    ENDIF 
  ENDDO
  
  Entropy = -e*2.d0*swidth

END SUBROUTINE 

