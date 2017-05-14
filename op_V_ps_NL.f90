SUBROUTINE op_V_ps_NL( Nstates, Vpsi )
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_PsPot, ONLY : NbetaNL, betaNL, w_NL
  USE m_atoms, ONLY : Natoms
  USE m_states, ONLY : Focc
  USE m_hamiltonian, ONLY : betaNL_psi
  IMPLICIT NONE 
  INTEGER :: Nstates
  REAL(8) :: Vpsi(Npoints,Nstates)
  INTEGER :: ia, ibeta, ist

  IF( NbetaNL <= 0 ) THEN 
    RETURN 
  ENDIF 

  Vpsi(:,:) = 0.d0

  DO ist = 1,Nstates
    DO ia = 1,Natoms
      DO ibeta = 1,NbetaNL
        Vpsi(:,ist) = w_NL(ibeta)*betaNL(:,ibeta)*betaNL_psi(ia,ist,ibeta)
      ENDDO 
    ENDDO 
  ENDDO 

END SUBROUTINE 



SUBROUTINE op_V_ps_NL_1col( ist, Vpsi )
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, dVol => LF3d_dVol
  USE m_PsPot, ONLY : NbetaNL, betaNL, w_NL
  USE m_atoms, ONLY : Natoms
  USE m_states, ONLY : Focc
  USE m_hamiltonian, ONLY : betaNL_psi
  IMPLICIT NONE 
  INTEGER :: ist
  REAL(8) :: Vpsi(Npoints)
  INTEGER :: ia, ibeta

  IF( NbetaNL <= 0 ) THEN 
    RETURN 
  ENDIF 

  Vpsi(:) = 0.d0

  DO ia = 1,Natoms
    DO ibeta = 1,NbetaNL
      Vpsi(:) = Vpsi(:) + w_NL(ibeta)*betaNL(:,ibeta)*betaNL_psi(ia,ist,ibeta)
    ENDDO 
  ENDDO 

END SUBROUTINE 

