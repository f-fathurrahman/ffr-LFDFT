SUBROUTINE init_evalsT()
  USE m_globals, ONLY : LF, evalsTx, evalsTy, evalsTz, N, Kprec
  IMPLICIT NONE
  INTEGER :: Nx, Ny, Nz
  REAL(8), ALLOCATABLE :: Kx(:,:), Ky(:,:), Kz(:,:)

  Nx = LF%LFx%N
  Ny = LF%LFy%N
  Nz = LF%LFz%N

  ALLOCATE( Kx(Nx,Nx), Ky(Ny,Ny), Kz(Nz,Nz) )
  ALLOCATE( evalsTx(Nx), evalsTy(Ny), evalsTz(Nz) )
  !ALLOCATE( Kprec(N**3) )

  Kx = -0.5d0*LF%LFx%D2jl
  Ky = -0.5d0*LF%LFy%D2jl
  Kz = -0.5d0*LF%LFz%D2jl

  CALL eig_dsyev( Kx, evalsTx, Nx )
  CALL eig_dsyev( Ky, evalsTy, Ny )
  CALL eig_dsyev( Kz, evalsTz, Nz )

  DEALLOCATE( Kx, Ky, Kz )
END SUBROUTINE


SUBROUTINE apply_PrecEvals( veci, veco, ic )
  USE m_globals, ONLY : LF, evalsTx, evalsTy, evalsTz, N, eVtau
  IMPLICIT NONE
  REAL(8) :: veci(N**3), veco(N**3)
  !
  INTEGER :: i, j, k, ip, Nx, Ny, Nz, ic

  Nx = LF%LFx%N
  Ny = LF%LFy%N
  Nz = LF%LFz%N
  !
  DO ip = 1, N**3
    i = LF%lin2xyz(1,ip)
    j = LF%lin2xyz(2,ip)
    k = LF%lin2xyz(3,ip)
    !
    !veco(ip) = veci(ip)/( eVtau(ic)/(evalsTx(i) + evalsTy(j) + evalsTz(k)) + 1.d0)
    !veco(ip) = veci(ip)/( eVtau(ic) + 1.d0)
    veco(ip) = veci(ip)/( abs(eVtau(ic) - (evalsTx(i) + evalsTy(j) + evalsTz(k))) + 1.d0)
  ENDDO

END SUBROUTINE




SUBROUTINE apply_PrecKin( N, vec, Kvec )
  USE m_globals, ONLY : LF
  IMPLICIT NONE
  INTEGER :: N
  REAL(8) :: vec(N), Kvec(N)
  !
  INTEGER :: i, j, k, ip, Nx, Ny, Nz
  REAL(8) :: Kinv

  Nx = LF%LFx%N
  Ny = LF%LFy%N
  Nz = LF%LFz%N
  !
  DO ip=1,N
    i = LF%lin2xyz(1,ip)
    j = LF%lin2xyz(2,ip)
    k = LF%lin2xyz(3,ip)
    !
    Kinv = -0.5d0*LF%LFx%D2jl(i,i) + -0.5d0*LF%LFy%D2jl(j,j) + -0.5d0*LF%LFz%D2jl(k,k)
    
    Kvec(ip) = vec(ip)/(Kinv + 1.d0)

  ENDDO
END SUBROUTINE



SUBROUTINE apply_invPrec( N, veci, veco, ic )
  USE m_globals, ONLY : LF, eVtau, Vpot
  IMPLICIT NONE
  INTEGER :: N, ic
  REAL(8) :: veci(N), veco(N)
  !
  INTEGER :: i, j, k, ip, ii, jj, kk, Nx, Ny, Nz

  Nx = LF%LFx%N
  Ny = LF%LFy%N
  Nz = LF%LFz%N
  !
  DO ip=1,N
    i = LF%lin2xyz(1,ip)
    j = LF%lin2xyz(2,ip)
    k = LF%lin2xyz(3,ip)
    !
    veco(ip) = veci(ip)  ! due to unit matrix
    !
    DO ii=1,Nx
      veco(ip) = veco(ip) + -0.5d0/eVtau(ic)*LF%LFx%D2jl(ii,i)*veci(LF%xyz2lin(ii,j,k))
    ENDDO
    !
    DO jj=1,Ny
      veco(ip) = veco(ip) + -0.5d0/eVtau(ic)*LF%LFy%D2jl(jj,j)*veci(LF%xyz2lin(i,jj,k))
    ENDDO
    !
    DO kk=1,Nz
      veco(ip) = veco(ip) + -0.5d0/eVtau(ic)*LF%LFz%D2jl(kk,k)*veci(LF%xyz2lin(i,j,kk))
    ENDDO
  ENDDO

END SUBROUTINE


SUBROUTINE apply_PrecCG( N, rho, phi, NiterMax, ic )
  USE m_globals, ONLY : LF
  IMPLICIT NONE
  INTEGER :: N, NiterMax, ic
  REAL(8) :: rho(N), phi(N)
  !
  REAL(8), ALLOCATABLE :: r(:), p(:) ! residual
  REAL(8), ALLOCATABLE :: nabla2(:)
  INTEGER :: ip, iter
  REAL(8) :: rsold, rsnew, alpha

  ALLOCATE( r(N), p(N), nabla2(N) )

  !
  DO ip=1,N
    phi(ip) = 0.d0
  ENDDO

  CALL apply_invPrec( N, phi, nabla2, ic )
  r(:) = rho(:) - nabla2(:)
  p(:) = r(:)

  !
  rsold = dot_product(r,r)
  !WRITE(*,*) 'rsold = ', rsold

  DO iter=1,NiterMax
    CALL apply_invPrec( N, p, nabla2, ic )
    !
    alpha = rsold/dot_product(p,nabla2)
    !
    phi(:) = phi(:) + alpha*p(:)
    !
    r(:) = r(:) - alpha*nabla2(:)
    !
    rsnew = dot_product(r,r)
    !WRITE(*,*) 'rsnew = ', sqrt(rsnew)
    !
    IF(sqrt(rsnew) < 1.d-8) THEN
    !IF(rsnew < 1.d-10) THEN
      WRITE(*,*) 'Convergence in apply_PrecCG in iter:', iter
      EXIT
    ENDIF
    p(:) = r(:) + (rsnew/rsold)*p(:)
    rsold = rsnew
  ENDDO

  WRITE(*,*) 'Final apply_invPrec: conv = ', sqrt(rsnew)

  DEALLOCATE( r, p, nabla2 )
END SUBROUTINE

