! 1-column version
SUBROUTINE apply_Ham(v, Hv)
  USE m_globals, ONLY : N, Vpot, Vhartree, Vxc
  IMPLICIT NONE
  !
  REAL(8) :: v(N**3)
  REAL(8) :: Hv(N**3)

  CALL apply_laplacian(v, Hv)
  !
  Hv(:) = -0.5d0*Hv(:) + ( Vpot(:) + Vhartree(:) + Vxc(:) )*v(:)
END SUBROUTINE


! for use in iterative diagonalization subroutines
! fixme npw is the same as npwx
SUBROUTINE h_psi( npwx, npw, nvec, psi, hpsi )
  IMPLICIT NONE
  INTEGER :: npwx, npw, nvec
  REAL(8) :: psi(npwx,nvec)
  REAL(8) :: hpsi(npwx,nvec)
  !
  INTEGER :: iv

  DO iv=1,nvec
    CALL apply_Ham( psi(:,iv), hpsi(:,iv) )
  ENDDO

END SUBROUTINE

