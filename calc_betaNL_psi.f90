SUBROUTINE calc_betaNL_psi( Nstates, psi )
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_PsPot
  IMPLICIT NONE 
  INTEGER :: Nstates
  REAL(8) :: psi(Npoints,Nstates)
  INTEGER :: ist, ibeta, ia
  REAL(8) :: ddot

  betaNL_psi(:,:,:) = 0.d0

  ibeta = 1
  ia = 1
  DO ist = 1,Nstates
    betaNL_psi(ia,ist,ibeta) = ddot( Npoints, betaNL(:,ibeta),1, psi(:,ist),1 ) * dVol
  ENDDO 

END SUBROUTINE 


