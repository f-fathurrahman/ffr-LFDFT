! The same as init_V_ps_loc_H_hgh, but initialize first in G-space
! the transformed to real space via FFT
!
SUBROUTINE init_V_ps_loc_H_hgh_G( Npoints, r, V )
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : G2 => LF3d_G2, &
                     Gv => LF3d_Gv, &
                     LL => LF3d_LL, &
                     NN => LF3d_NN
  IMPLICIT NONE
  INTEGER :: Npoints
  REAL(8) :: r(Npoints)
  REAL(8) :: V(Npoints)
  !
  INTEGER :: ig
  REAL(8) :: Vol
  INTEGER :: Ng, Nx, Ny, Nz
  REAL(8) :: z_val
  REAL(8) :: rlocal
  REAL(8) :: c(2)
  REAL(8) :: pre1, pre2, Gr, expGr2
  ! temporary array for FFT
  COMPLEX(8), ALLOCATABLE :: ctmp(:)
  ! structure factor
  COMPLEX(8), ALLOCATABLE :: strf(:,:)
  ! some hardcoded parameters, for structure factor calculation
  INTEGER, PARAMETER :: Nspecies = 1
  INTEGER, PARAMETER :: Natoms = 1
  REAL(8) :: Xpos(3,Natoms)
  INTEGER :: atm2species(Natoms)

  Xpos(:,1) = 0.5d0*LL
  WRITE(*,*) 'Xpos = ', Xpos

  atm2species(1) = 1

  ALLOCATE( ctmp(Npoints) )
  ALLOCATE( strf(Npoints,Nspecies) )

  Ng = Npoints
  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  CALL calc_strfact( Natoms, Xpos, Nspecies, atm2species, Ng, Gv, strf )
  WRITE(*,*) 'sum(strf) = ', sum(strf)

  z_val = 1.d0
  rlocal = 0.2d0
  c(1) = -4.180237d0
  c(2) = 0.725075d0

  Vol = LL(1)*LL(2)*LL(3)
  pre1 = -4.d0*PI*z_val/Vol
  pre2 = sqrt(8.d0*PI**3)* (rlocal**3) / Vol

  ! In the current implementatation Ng = Npoints
  ctmp(1) = cmplx(0.0,0.0,kind=8)
  DO ig = 2, Npoints
    Gr = sqrt( G2(ig) )*rlocal
    expGr2 = exp(-0.5*Gr**2)
    ctmp(ig) = pre1/G2(ig)*expGr2 + pre2*expGr2 * ( c(1) + c(2)*(3d0 - Gr**2) )
  ENDDO
  
  DO ig = 1, Npoints
    ctmp(ig) = ctmp(ig)*strf(ig,1)
  ENDDO
  
  ! inverse FFT
  CALL fft_fftw3( ctmp, Nx, Ny, Nz, .true. )

  ! XXX Need Npoints ?
  V(:) = real(ctmp(:), kind=8)

  DEALLOCATE( ctmp )
  DEALLOCATE( strf )

END SUBROUTINE

