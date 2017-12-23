SUBROUTINE normalize_rhoe( Npoints, rhoe )
  USE m_constants, ONLY : EPS_SMALL
  USE m_states, ONLY : Nelectrons
  USE m_LF3d, ONLY : dVol => LF3d_dVol
  IMPLICIT NONE 
  INTEGER :: Npoints
  REAL(8) :: Rhoe(Npoints)
  !
  INTEGER :: ip
  REAL(8) :: integRho

  DO ip = 1,Npoints
    IF( rhoe(ip) < EPS_SMALL ) THEN 
      rhoe(ip) = EPS_SMALL
    ENDIF 
  ENDDO 

!!> Need to rescale electron density such that it integrates to number of electrons.
    integRho = sum(Rhoe)*dVol
    IF( abs(integRho - Nelectrons) > 1.0d-6 ) THEN
      WRITE(*,*)
      WRITE(*,'(1x,A,ES18.10)') 'WARNING: diff after mix rho = ', abs(integRho-Nelectrons)
      WRITE(*,*) 'Rescaling Rho'
      Rhoe(:) = Nelectrons/integRho * Rhoe(:)
      integRho = sum(Rhoe)*dVol
      WRITE(*,'(1x,A,F18.10)') 'After rescaling: integRho = ', integRho
    ENDIF 

END SUBROUTINE 

