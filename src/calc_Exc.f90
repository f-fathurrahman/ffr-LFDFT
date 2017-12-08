SUBROUTINE calc_Exc()
  USE m_energies, ONLY : E_xc
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, dVol => LF3d_dVol
  USE m_hamiltonian, ONLY : Rhoe 
  USE m_xc
  IMPLICIT NONE 

  ! these calls should only be done for LDA
  IF( XC_NAME == 'VWN' ) THEN 
    CALL excVWN( Npoints, Rhoe, EPS_XC )
    E_xc = sum( Rhoe(:) * EPS_XC(:) )*dVol
  ENDIF 

END SUBROUTINE 

