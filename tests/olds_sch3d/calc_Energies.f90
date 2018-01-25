! psi are assumed to be orthogonalized
!------------------------------------------------
SUBROUTINE calc_Energies( psi, Ekin, Epot, Etot )
!------------------------------------------------
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_states, ONLY : Nstates
  USE m_hamiltonian, ONLY : V_ps_loc
  IMPLICIT NONE
  !
  REAL(8) :: psi(Npoints, Nstates)
  REAL(8) :: nabla2_psi(Npoints)
  REAL(8) :: Etot, Ekin, Epot
  !
  INTEGER :: ist
  !
  REAL(8) :: ddot

  Etot = 0.d0
  Ekin = 0.d0
  Epot = 0.d0
  ! assume all occupancies are 1.d0
  DO ist = 1, Nstates
    CALL op_nabla2( psi(:,ist), nabla2_psi(:) )
    Ekin = Ekin + (-0.5d0*ddot( Npoints, psi(:,ist),1, nabla2_psi(:),1 ) * dVol)
    Epot = Epot + sum( V_ps_loc(:)*psi(:,ist)**2 ) * dVol
  ENDDO
  Etot = Ekin + Epot
END SUBROUTINE
