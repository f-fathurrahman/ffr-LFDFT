SUBROUTINE normalize_rhoe( Npoints, rhoe )
  USE m_constants, ONLY : EPS_SMALL
  IMPLICIT NONE 
  INTEGER :: Npoints
  REAL(8) :: rhoe(Npoints)
  INTEGER :: ip

  DO ip = 1,Npoints
    IF( rhoe(ip) < EPS_SMALL ) THEN 
      rhoe(ip) = EPS_SMALL
    ENDIF 
  ENDDO 

END SUBROUTINE 

