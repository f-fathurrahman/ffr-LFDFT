! efefer 7 January 2016

! Testing solution of Poisson equation using interpolating scaling method
! proposed by L. Genovese

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


SUBROUTINE solve_poisson_ISF( rho_LF )
  USE m_constants
  USE gbl_poisson
  IMPLICIT NONE
  !
  REAL(8) :: rho_LF(3,Nx*Ny*Nz)
  !
  REAL(8), ALLOCATABLE :: tmpgrid(:,:,:)
  REAL(8), ALLOCATABLE :: rhophi(:,:,:)
  INTEGER :: ip
  INTEGER :: i, j, k
  REAL(8) :: Lx, Ly, Lz


  ALLOCATE( tmpgrid(Nx,Ny,Nz) )
  ALLOCATE( rhophi(Nx,Ny,Nz) )

  DO k=1,Nz
    DO j=1,Ny
      DO i=1,Nx
        tmpgrid(i,j,k) = 0.5d0*Lx*(2*i-1)/Nx + 0.5d0*Ly*(2*j-1)/Ny
                         0.5d0*Lz*(2*k-1)/Nz
      ENDDO
    ENDDO
  ENDDO

  DO ip = 1,LF3d%N
    i = lin2xyz(1,ip)
    j = lin2xyz(2,ip)
    k = lin2xyz(3,ip)
    rhophi(i,j,k) = rho_LF(ip)
  ENDDO

  DEALLOCATE( rhophi, tmpgrid )

END SUBROUTINE



PROGRAM t_poisson
  USE m_constants
  USE m_LF3d
  USE gbl_poisson
  IMPLICIT NONE
  !

  REAL(8), ALLOCATABLE :: rho_LF(:), rho(:) ! density
  REAL(8), ALLOCATABLE :: phi(:) ! potential
  !
  REAL(8) :: sigma1, sigma2, r, x0, y0, z0, deltaV
  INTEGER :: ip
  REAL(8) :: Uana, Unum

  ! Setup system
  Nx = 39
  Ny = 39
  Nz = 39
  !
  Lx = 16.d0
  Ly = 16.d0
  Lz = 16.d0
  !
  CALL init_LF3d_p( LF3d, Nx,Ny,Nz, Lx,Ly,Lz )
  
  ALLOCATE( rho_LF(Nx*Ny*Nz) )
  ALLOCATE( phi(Nx*Ny*Nz) )

  ! center of the box
  x0 = Lx/2.d0
  y0 = Ly/2.d0
  z0 = Lz/2.d0
  ! Initialize
  sigma1 = 0.75d0
  sigma2 = 0.50d0
  DO ip = 1, LF3d%N
    !FIXME add some small positive value to avoid zero distance
    r = norm2( LF3d%lingrid(:,ip) - (/x0,y0,z0/) )
    rho_LF(ip) = exp(-r**2/(2*sigma2**2))/(2*pi*sigma2**2)**1.5d0 - &
              exp(-r**2/(2*sigma1**2))/(2*pi*sigma1**2)**1.5d0
  ENDDO

  deltaV = LF3d%LFx%h * LF3d%LFy%h * LF3d%LFz%h
  WRITE(*,*) 'deltaV = ', deltaV
  WRITE(*,*) sum( rho(:) )*deltaV

  ! Solve Poisson equation
  !CALL solve_poisson_pcg( LF3d%N, -4.d0*PI*rho, phi )

  !
  Unum = 0.5d0*sum( rho_LF(:)*phi(:) )*deltaV
  Uana = ( (1.d0/sigma1 + 1.d0/sigma2)/2.d0 - sqrt(2.d0)/sqrt(sigma1**2 + sigma2**2) )/sqrt(PI)
  WRITE(*,*) 'Unum = ', Unum
  WRITE(*,*) 'Uana = ', Uana

  DEALLOCATE( rho, phi )
END PROGRAM

