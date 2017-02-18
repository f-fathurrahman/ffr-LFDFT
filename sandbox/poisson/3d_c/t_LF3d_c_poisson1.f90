! efefer 1 January 2016

MODULE gbl_poisson
  USE m_LF3d
  IMPLICIT NONE
  INTEGER :: Nx, Ny, Nz
  REAL(8) :: Lx, Ly, Lz
  TYPE(LF3d_t) :: LF3d
END MODULE


! apply Laplacian to input vector v
SUBROUTINE apply_laplacian(N, v, nabla2_v)
  USE gbl_poisson
  IMPLICIT NONE
  INTEGER :: N
  REAL(8) :: v(N)
  REAL(8) :: nabla2_v(N)
  !
  INTEGER :: i, j, k, ip, ii, jj, kk
  !
  ! N should be the same as LF3d%N
  DO ip=1,N
    i = LF3d%lin2xyz(1,ip)
    j = LF3d%lin2xyz(2,ip)
    k = LF3d%lin2xyz(3,ip)
    !
    nabla2_v(ip) = 0.d0
    !
    DO ii=1,Nx
      nabla2_v(ip) = nabla2_v(ip) + LF3d%LFx%D2jl(ii,i)*v(LF3d%xyz2lin(ii,j,k))
    ENDDO
    !
    DO jj=1,Ny
      nabla2_v(ip) = nabla2_v(ip) + LF3d%LFy%D2jl(jj,j)*v(LF3d%xyz2lin(i,jj,k))
    ENDDO
    !
    DO kk=1,Nz
      nabla2_v(ip) = nabla2_v(ip) + LF3d%LFz%D2jl(kk,k)*v(LF3d%xyz2lin(i,j,kk))
    ENDDO
  ENDDO

END SUBROUTINE


SUBROUTINE solve_poisson_cg(N, rho, phi)
  USE gbl_poisson
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: rho(N), phi(N)
  !
  REAL(8), ALLOCATABLE :: r(:), p(:) ! residual
  REAL(8), ALLOCATABLE :: nabla2_phi(:)
  INTEGER :: ip, iter
  REAL(8) :: rsold, rsnew, alpha

  ALLOCATE( r(N), p(N), nabla2_phi(N) )

  !
  DO ip=1,N
    CALL random_number( phi(ip) )
  ENDDO

  CALL apply_laplacian( N, phi, nabla2_phi )
  r(:) = rho(:) - nabla2_phi(:)
  p(:) = r(:)

  !
  rsold = dot_product(r,r)
  WRITE(*,*) 'rsold = ', rsold

  DO iter=1,N
    CALL apply_laplacian( N, p, nabla2_phi )
    !
    alpha = rsold/dot_product(p,nabla2_phi)
    !
    phi(:) = phi(:) + alpha*p(:)
    !
    r(:) = r(:) - alpha*nabla2_phi(:)
    !
    rsnew = dot_product(r,r)
    WRITE(*,*) 'rsnew = ', sqrt(rsnew)
    !
    IF(sqrt(rsnew) < 1.d-10) THEN
    !IF(rsnew < 1.d-10) THEN
      WRITE(*,*) 'Convergence in solve_poisson_cg in iter:', iter
      EXIT
    ENDIF
    p(:) = r(:) + (rsnew/rsold)*p(:)
    rsold = rsnew
  ENDDO

  DEALLOCATE( r, p, nabla2_phi )
END SUBROUTINE



PROGRAM t_poisson
  USE m_constants
  USE m_LF3d
  USE gbl_poisson
  IMPLICIT NONE
  !

  REAL(8), ALLOCATABLE :: rho(:) ! density
  REAL(8), ALLOCATABLE :: phi(:) ! potential
  !
  REAL(8) :: sigma1, sigma2, r, x0, y0, z0, deltaV
  INTEGER :: ip
  REAL(8) :: Uana, Unum

  Nx = 37
  Ny = 37
  Nz = 37
  !
  Lx = 16.d0
  Ly = 16.d0
  Lz = 16.d0
  !
  CALL init_LF3d_c( LF3d, (/Nx,Ny,Nz/), (/0.d0,0.d0,0.d0/), (/Lx,Ly,Lz/) )
  
  ALLOCATE( rho(Nx*Ny*Nz) )
  ALLOCATE( phi(Nx*Ny*Nz) )

  ! center of the box
  x0 = Lx/2.d0
  y0 = Ly/2.d0
  z0 = Lz/2.d0
  ! Initialize
  sigma1 = 0.75d0
  sigma2 = 0.50d0
  DO ip = 1, LF3d%N
    r = norm2( LF3d%lingrid(:,ip) - (/x0,y0,z0/) )
    rho(ip) = exp(-r**2/(2*sigma2**2))/(2*pi*sigma2**2)**1.5d0 - &
              exp(-r**2/(2*sigma1**2))/(2*pi*sigma1**2)**1.5d0
    !WRITE(*,'(1x,I5,2F18.10)') ip, r, rho(ip)
  ENDDO

  ! For cluster LF
  deltaV = LF3d%LFx%h * LF3d%LFy%h * LF3d%LFz%h
  WRITE(*,*) 'deltaV = ', deltaV
  WRITE(*,*) sum( rho(:) )*deltaV

  ! Solve Poisson equation
  CALL solve_poisson_cg( LF3d%N, -4.d0*PI*rho, phi )

  !
  Unum = 0.5d0*sum( rho(:)*phi(:) )*deltaV
  Uana = ( (1.d0/sigma1 + 1.d0/sigma2)/2.d0 - sqrt(2.d0)/sqrt(sigma1**2 + sigma2**2) )/sqrt(PI)
  WRITE(*,*) 'Unum = ', Unum
  WRITE(*,*) 'Uana = ', Uana

  DEALLOCATE( rho, phi )
END PROGRAM

