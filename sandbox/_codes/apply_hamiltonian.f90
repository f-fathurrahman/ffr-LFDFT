
! 1-column version
SUBROUTINE apply_Ham(v, Hv)
  USE m_globals, ONLY : N, Vpot, Vhartree
  IMPLICIT NONE
  !
  REAL(8) :: v(N**3)
  REAL(8) :: Hv(N**3)

  CALL apply_laplacian(v, Hv)
  !
  Hv(:) = -0.5d0*Hv(:) + (Vpot(:) + Vhartree(:))*v(:)
END SUBROUTINE


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


! apply Laplacian to input vector v
SUBROUTINE apply_laplacian(v, nabla2_v)
  USE m_globals, ONLY : N, LF
  IMPLICIT NONE
  REAL(8) :: v(N**3)
  REAL(8) :: nabla2_v(N**3)
  !
  INTEGER :: i, j, k, ip, ii, jj, kk
  !
  ! N**3 should be the same as LF3d%N
  DO ip=1,N**3
    i = LF%lin2xyz(1,ip)
    j = LF%lin2xyz(2,ip)
    k = LF%lin2xyz(3,ip)
    !
    nabla2_v(ip) = 0.d0
    !
    DO ii=1,N
      nabla2_v(ip) = nabla2_v(ip) + LF%LFx%D2jl(ii,i)*v(LF%xyz2lin(ii,j,k))
    ENDDO
    !
    DO jj=1,N
      nabla2_v(ip) = nabla2_v(ip) + LF%LFy%D2jl(jj,j)*v(LF%xyz2lin(i,jj,k))
    ENDDO
    !
    DO kk=1,N
      nabla2_v(ip) = nabla2_v(ip) + LF%LFz%D2jl(kk,k)*v(LF%xyz2lin(i,j,kk))
    ENDDO
  ENDDO

END SUBROUTINE
