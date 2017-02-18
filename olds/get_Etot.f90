
! psi are assumed to be orthogonalized
!-----------------------------------------------
SUBROUTINE get_Etot_v2(Ncol, psi, Ekin, Epot, Etot)
!-----------------------------------------------
  USE m_globals, ONLY : N, Vpsloc, LF, Ehartree
  IMPLICIT NONE
  !
  INTEGER :: Ncol
  REAL(8) :: psi(N**3,Ncol)
  REAL(8) :: nabla2_psi(N**3)
  REAL(8) :: Etot, Ekin(Ncol), Epot
  !
  INTEGER :: ic
  REAL(8) :: deltaV
  !
  REAL(8) :: ddot

  deltaV = LF%LFx%h * LF%LFy%h * LF%LFz%h
  !
  Etot = 0.d0
  Ekin = 0.d0
  Epot = 0.d0
  ! assume all occupancies are 1.d0
  DO ic=1,Ncol
    CALL apply_laplacian( psi(:,ic), nabla2_psi(:) )
    Ekin(ic) = -0.5d0*ddot( N**3, psi(:,ic),1, nabla2_psi(:),1 )
    !
    Epot = Epot + sum( Vpsloc(:)*psi(:,ic)**2 ) ! FIXME: We don't need to multiply to deltaV here
  ENDDO
  Etot = sum(Ekin) + Epot + Ehartree
END SUBROUTINE




! psi are assumed to be orthogonalized
!-----------------------------------------------
SUBROUTINE get_Etot(Ncol, psi, Ekin, Epot, Etot)
!-----------------------------------------------
  USE m_globals, ONLY : N, Vpsloc, LF, Ehartree
  IMPLICIT NONE
  !
  INTEGER :: Ncol
  REAL(8) :: psi(N**3,Ncol)
  REAL(8) :: nabla2_psi(N**3)
  REAL(8) :: Etot, Ekin, Epot
  !
  INTEGER :: ic
  REAL(8) :: deltaV
  !
  REAL(8) :: ddot

  deltaV = LF%LFx%h * LF%LFy%h * LF%LFz%h
  !
  Etot = 0.d0
  Ekin = 0.d0
  Epot = 0.d0
  ! assume all occupancies are 1.d0
  DO ic=1,Ncol
    CALL apply_laplacian( psi(:,ic), nabla2_psi(:) )
    Ekin = Ekin + -0.5d0*ddot( N**3, psi(:,ic),1, nabla2_psi(:),1 )
    !
    Epot = Epot + sum( Vpsloc(:)*psi(:,ic)**2 ) ! FIXME: We don't need to multiply to deltaV here
  ENDDO
  Etot = Ekin + Epot + Ehartree
END SUBROUTINE


