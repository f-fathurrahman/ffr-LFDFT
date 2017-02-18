! efefer 9 May 2016

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
    !CALL random_number( phi(ip) )
    phi(ip) = 0.d0
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


FUNCTION eval_rho(r) RESULT(rho)
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  REAL(8) :: r, rho
  !
  !rho = 8.d0*PI*( 3.d0 - 2.d0*r**2 ) * exp(-r**2)
  rho = ( 4.d0*r**2 - 6.d0 ) * exp(-r**2)
END FUNCTION


FUNCTION eval_phi(r) RESULT(phi)
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  REAL(8) :: r, phi
  phi = -4.d0*pi*exp(-r**2)
END


PROGRAM t_poisson
  USE m_constants
  USE m_LF3d
  USE gbl_poisson
  IMPLICIT NONE
  !

  REAL(8), ALLOCATABLE :: rho(:) ! density
  REAL(8), ALLOCATABLE :: phi(:) ! potential
  REAL(8), ALLOCATABLE :: phi_ana(:)
  !
  REAL(8) :: sigma, r, x0, y0, z0, deltaV
  INTEGER :: ip, ix, iy, iz
  REAL(8) :: Unum
  REAL(8) :: eval_rho, eval_phi

  Nx = 47
  Ny = 47
  Nz = 47
  !
  Lx = 14.d0
  Ly = 14.d0
  Lz = 14.d0
  !
  CALL init_LF3d_c( LF3d, (/Nx,Ny,Nz/), (/0.d0,0.d0,0.d0/), (/Lx,Ly,Lz/) )
  
  ALLOCATE( rho(Nx*Ny*Nz) )
  ALLOCATE( phi(Nx*Ny*Nz) )
  ALLOCATE( phi_ana(Nx*Ny*Nz) )

  ! center of the box
  x0 = Lx/2.d0
  y0 = Ly/2.d0
  z0 = Lz/2.d0

  ! Initialize charge density and analytic solution
  sigma = 1.d0
  DO ip = 1, LF3d%N
    r = norm2( LF3d%lingrid(:,ip) - (/x0,y0,z0/) )
    rho(ip) = eval_rho(r)
    phi_ana(ip) = eval_phi(r)
  ENDDO

  ! For cluster LF
  deltaV = LF3d%LFx%h * LF3d%LFy%h * LF3d%LFz%h
  WRITE(*,*) 'deltaV = ', deltaV
  WRITE(*,*) sum( rho(:) )*deltaV

  ! Solve Poisson equation
  CALL solve_poisson_cg( LF3d%N, -4.d0*PI*rho, phi )

  iy = Ny/2
  iz = Nz/2
  DO ix = 1, Nx
    ip = LF3d%xyz2lin(ix,iy,iz)
    WRITE(444,'(1x,4F18.10)') LF3d%lingrid(1,ip), rho(ip), phi(ip), phi_ana(ip)
  ENDDO
  WRITE(*,*) 'sum(phi-phi_ana)', sum(phi-phi_ana)

  !
  Unum = 0.5d0*sum( rho(:)*phi(:) )*deltaV
  WRITE(*,*) 'Unum = ', Unum
  WRITE(*,'(1x,A,2I5)') 'iy, iz = ', iy, iz

  DEALLOCATE( rho, phi, phi_ana )
END PROGRAM

