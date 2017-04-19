!! PURPOSE:
!!
!!   This subroutine has the same purpose as `init_V_ps_loc_H_hgh`,
!!   but initialize first in G-space instead of R-space.
!!   The usual real space representation is obtained by the transforming
!!   it to real space via FFT.
!! 
!! AUTHOR:
!!
!!   Fadjar Fathurrahman
!!
!! NOTE:
!!
!!   This should be used only in the periodic case
!!
SUBROUTINE init_V_ps_loc_H_hgh_G( Npoints, V )

  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : G2 => LF3d_G2, &
                     LL => LF3d_LL, &
                     NN => LF3d_NN
  USE m_atoms, ONLY : strf => StructureFactor
  IMPLICIT NONE
  !! Number of points
  INTEGER :: Npoints

  !! The potential (in real space)
  REAL(8) :: V(Npoints)

  INTEGER :: ig
  REAL(8) :: Vol  ! unit cell volume
  INTEGER :: Nx, Ny, Nz
  REAL(8) :: z_val
  REAL(8) :: rlocal
  REAL(8) :: c(2)
  REAL(8) :: pre1, pre2, Gr, expGr2
  ! temporary array for FFT
  COMPLEX(8), ALLOCATABLE :: ctmp(:)

  ALLOCATE( ctmp(Npoints) )

  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  z_val = 1.d0
  rlocal = 0.2d0
  c(1) = -4.180237d0
  c(2) = 0.725075d0

  Vol = LL(1)*LL(2)*LL(3)
  pre1 = -4.d0*PI*z_val/Vol
  pre2 = sqrt(8.d0*PI**3)* (rlocal**3) / Vol

  ! In the current implementatation Ng = Npoints
  ctmp(1) = cmplx(0.0,0.0,kind=8)
  DO ig = 1, Npoints
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

