SUBROUTINE calc_Ewald( )

  USE m_constants, ONLY : PI
  USE m_atoms, ONLY : Natoms, Nspecies, atm2species, &
                      Zv => AtomicValences, &
                      strf => StructureFactor
  USE m_LF3d, ONLY : G2 => LF3d_G2, &
                     Npoints => LF3d_Npoints, &
                     NN => LF3d_NN, &
                     LL => LF3d_LL, &
                     lingrid => LF3d_lingrid, &
                     dVol => LF3d_dVol
  USE m_energies, ONLY : E_nn
  IMPLICIT NONE 
  !
  REAL(8), ALLOCATABLE :: sigma(:)
  REAL(8), ALLOCATABLE :: dr(:)
  REAL(8) :: cx, cy, cz, dx, dy, dz
  REAL(8) :: c1, cc1
  INTEGER :: Nx, Ny, Nz
  INTEGER :: ip, isp, ia
  REAL(8) :: gchg, intrho
  REAL(8), ALLOCATABLE :: rho_is(:,:)
  REAL(8), ALLOCATABLE :: Rho(:), phi(:)
  COMPLEX(8), ALLOCATABLE :: ctmp(:)
  REAL(8) :: E_H, E_self
  
  WRITE(*,*) 'Calculating Ewald energy'

  ALLOCATE( sigma(Nspecies) )
  sigma(:) = 0.25d0   !!! DEFAULT !!!

  cx = 0.5d0*LL(1)
  cy = 0.5d0*LL(2)
  cz = 0.5d0*LL(3)

  ALLOCATE( dr(Npoints) )
  DO ip = 1,Npoints
    dx = lingrid(1,ip) - cx
    dy = lingrid(2,ip) - cy
    dz = lingrid(3,ip) - cz
    dr(ip) = sqrt( dx**2 + dy**2 + dz**2 )
  ENDDO

  ALLOCATE( rho_is(Npoints,Nspecies) )
  ALLOCATE( Rho(Npoints) )
  ALLOCATE( phi(Npoints) )
  ALLOCATE( ctmp(Npoints) )

  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  Rho(:) = 0.d0

  DO isp = 1, Nspecies
    !
    c1 = 2.d0*sigma(isp)**2
    cc1 = sqrt(2.d0*PI*sigma(isp)**2)**3
    !
    DO ip = 1,Npoints
      gchg = Zv(isp) * exp( -dr(ip)**2/c1 ) / cc1
      ctmp(ip) = cmplx( gchg, 0.d0, kind=8 )
    ENDDO
    intrho = real( sum(ctmp)*dVol )
    WRITE(*,*) 'ctmp initial: ', isp, intrho
    !
    CALL fft_fftw3( ctmp, Nx, Ny, Nz, .false. )  ! to G-space
    !
    DO ip = 1,Npoints
      ctmp(ip) = ctmp(ip)*strf(ip,isp)
    ENDDO 
    !
    CALL fft_fftw3( ctmp, Nx, Ny, Nz, .true. )
    DO ip = 1, Npoints
      rho_is(ip,isp) = real( ctmp(ip), kind=8 )
    ENDDO 
    !
    intrho = sum( rho_is(:,isp) ) * dVol
    WRITE(*,*) 'isp, intrho = ', isp, intrho
    !
    Rho(:) = Rho(:) + rho_is(:,isp)
  ENDDO 

  intrho = sum(Rho)*dVol
  WRITE(*,*) 'Total intrho = ', intrho

  ! Solve Poisson equation
  DO ip = 1,Npoints
    ctmp(ip) = cmplx( Rho(ip), 0.d0, kind=8 )
  ENDDO 
  CALL fft_fftw3( ctmp, Nx, Ny, Nz, .false. )
  !
  ctmp(1) = cmplx(0.d0,0.d0,kind=8)
  DO ip = 2,Npoints
    ctmp(ip) = 4.d0*PI*ctmp(ip)/G2(ip)
  ENDDO 
  !
  CALL fft_fftw3( ctmp, Nx, Ny, Nz, .true. )
  DO ip = 1,Npoints
    phi(ip) = real( ctmp(ip), kind=8 )
  ENDDO 
  E_H = 0.5*sum( phi(:)*Rho(:) ) * dVol

  E_self = 0.d0
  DO ia = 1,Natoms
    isp = atm2species(ia)
    E_self = E_self + Zv(isp)**2 / (2.d0*sqrt(PI)) * (1.d0/sigma(isp))
  ENDDO 

  E_nn = E_H - E_self

  WRITE(*,*) 'E_nn = ', E_nn

  DEALLOCATE( phi )
  DEALLOCATE( ctmp )
  DEALLOCATE( Rho )
  DEALLOCATE( rho_is )
  DEALLOCATE( dr )
  DEALLOCATE( sigma )

END SUBROUTINE 

