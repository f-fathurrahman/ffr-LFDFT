SUBROUTINE calc_betaNL_psi( Nstates, psi )
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_PsPot, ONLY : NbetaNL, betaNL
  USE m_hamiltonian, ONLY : betaNL_psi
  USE m_atoms, ONLY : Natoms
  IMPLICIT NONE 
  INTEGER :: Nstates
  REAL(8) :: psi(Npoints,Nstates)
  INTEGER :: ist, ibeta, ia
  REAL(8) :: ddot

  ! immediate return if no projectors are available
  IF( NbetaNL <= 0 ) THEN
    RETURN 
  ENDIF 

  betaNL_psi(:,:,:) = 0.d0

  DO ia = 1,Natoms
    DO ist = 1,Nstates
      DO ibeta = 1,NbetaNL
        betaNL_psi(ia,ist,ibeta) = ddot( Npoints, betaNL(:,ibeta),1, psi(:,ist),1 ) * dVol
      ENDDO 
    ENDDO 
  ENDDO 

END SUBROUTINE 


