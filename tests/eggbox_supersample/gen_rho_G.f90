SUBROUTINE gen_rho_G( center, length, rho_gauss )

  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     G2 => LF3d_G2, &
                     NN => LF3d_NN, &
                     LL => LF3d_LL, &
                     Gv => LF3d_Gv, &
                     grid_x => LF3d_grid_x, &
                     grid_y => LF3d_grid_y, &
                     grid_z => LF3d_grid_z
  IMPLICIT NONE 
  !
  REAL(8) :: center(3)
  REAL(8) :: length
  REAL(8) :: rho_gauss(Npoints)
  !
  INTEGER :: ip, Nx, Ny, Nz, ig
  REAL(8) :: Omega, ff
  REAL(8) :: GX, shiftx, shifty, shiftz
  INTEGER :: Ngvec
  COMPLEX(8), ALLOCATABLE :: ctmp(:)
  COMPLEX(8), ALLOCATABLE :: strf(:)

  Ngvec = Npoints
  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  ! Cell volume
  Omega = LL(1) * LL(2) * LL(3)

  ALLOCATE( ctmp(Ngvec) )
  ALLOCATE( strf(Ngvec) )

  shiftx = 0.5d0*( grid_x(2) - grid_x(1) )
  shifty = 0.5d0*( grid_y(2) - grid_y(1) )
  shiftz = 0.5d0*( grid_z(2) - grid_z(1) )

  strf(:) = cmplx(0.d0,0.d0,kind=8)
  DO ig = 1,Ngvec
    GX = (center(1)-shiftx)*Gv(1,ig) + (center(2)-shifty)*Gv(2,ig) + (center(3)-shiftz)*Gv(3,ig)
    strf(ig) = strf(ig) + cmplx( cos(GX), -sin(GX), kind=8 )
  ENDDO 

  rho_gauss(:) = 0.d0

  ctmp(:) = cmplx(0.d0,0.d0,kind=8)
  DO ig = 2,Ngvec
    ctmp(ig) = exp(-length**2*G2(ig)) * strf(ig) / Omega
  ENDDO

  ! inverse FFT: G -> R
  CALL fft_fftw3( ctmp, Nx, Ny, Nz, .true. )
  ! XXX: Move this outside isp loop ?
  DO ip = 1,Npoints
    ff = real(ctmp(ip),kind=8)
    IF( ff > 0.d0 ) THEN 
      rho_gauss(ip) = ff
    ELSE 
      rho_gauss(ip) = 0.d0
    ENDIF 
  ENDDO 

  DEALLOCATE( ctmp )
  DEALLOCATE( strf )

END SUBROUTINE 

