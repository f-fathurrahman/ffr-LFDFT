FUNCTION calc_entropy( Focc, Tbeta ) RESULT(e)
  USE m_constants, ONLY : EPS_SMALL
  USE m_states, ONLY : Nstates
  IMPLICIT NONE 
  REAL(8) :: Focc(Nstates)
  REAL(8) :: Tbeta
  INTEGER :: ist
  REAL(8) :: e

  ! nonzero occupations
  !inz = find( occ > EPS_SMALL )
  !occnz = occ(inz)
  !e = sum( occnz*log(occnz) )
  e = 0.d0
  DO ist = 1,Nstates
    ! XXX: probably need to use abs(Focc) here to anticipate negative occupation numbers ?
    IF( Focc(ist) > EPS_SMALL ) THEN 
      ! non-zero occupation
      e = e + Focc(ist)*log(Focc(ist))
    ELSE 
      ! zero occupation
      e = e + (1.d0 - Focc(ist)) * log(1.d0 - Focc(ist))
    ENDIF 
  ENDDO 

  e = e / Tbeta
END FUNCTION 

