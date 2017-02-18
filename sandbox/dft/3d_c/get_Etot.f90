! psi are assumed to be orthogonalized
!----------------------------------------------------------
SUBROUTINE get_Etot(Nbasis, Nstates, psi, Ekin, Epot, Etot)
!----------------------------------------------------------
  USE m_globals, ONLY : Vpot, LF, Ehartree, Exc, deltaV, Focc
  IMPLICIT NONE
  !
  INTEGER :: Nbasis, Nstates
  REAL(8) :: psi(Nbasis,Nstates)
  REAL(8) :: nabla2_psi(Nbasis)
  REAL(8) :: Etot, Ekin, Epot
  !
  INTEGER :: is
  !
  REAL(8) :: ddot

  Etot = 0.d0
  Ekin = 0.d0
  Epot = 0.d0
  !
  DO is=1,Nstates
    CALL apply_laplacian( psi(:,is), nabla2_psi(:) )
    Ekin = Ekin + -0.5d0*ddot( Nbasis, psi(:,is),1, nabla2_psi(:),1 ) * Focc(is) 
    !
    ! FIXME: We don't need to multiply to deltaV here
    Epot = Epot + sum( Vpot(:) *psi(:,is)**2 )*Focc(is)
  ENDDO
  Etot = Ekin + Epot + Ehartree + Exc
END SUBROUTINE

