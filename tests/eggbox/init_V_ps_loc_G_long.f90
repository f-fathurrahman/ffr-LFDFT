SUBROUTINE init_V_ps_loc_G_long( V_long )
  
  USE m_PsPot, ONLY : Ps => Ps_HGH_Params

  USE m_Ps_HGH, ONLY : hgh_eval_Vloc_G_long

  USE m_atoms, ONLY : strf => StructureFactor, Nspecies

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     G2 => LF3d_G2, &
                     NN => LF3d_NN, &
                     LL => LF3d_LL

  IMPLICIT NONE 
  !
  REAL(8) :: V_long(Npoints)
  !
  INTEGER :: ip, isp, Nx, Ny, Nz
  REAL(8) :: Gm, Omega
  COMPLEX(8), ALLOCATABLE :: ctmp(:)

  WRITE(*,*)
  WRITE(*,*) 'V_ps_loc_long is initialized via G-space'

  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  ! Cell volume
  Omega = LL(1) * LL(2) * LL(3)

  ALLOCATE( ctmp(Npoints) )

  V_long(:) = 0.d0

  DO isp = 1,Nspecies

    ctmp(:) = cmplx(0.d0,0.d0,kind=8)
    DO ip = 1,Npoints
      Gm = sqrt(G2(ip))
      ctmp(ip) = hgh_eval_Vloc_G_long( Ps(isp), Gm ) * strf(ip,isp) / Omega
    ENDDO

    ! inverse FFT: G -> R
    CALL fft_fftw3( ctmp, Nx, Ny, Nz, .true. )

    ! XXX: Move this outside isp loop ?
    DO ip = 1,Npoints
      V_long(ip) = V_long(ip) + real( ctmp(ip), kind=8 )
    ENDDO 

  ENDDO 

  DEALLOCATE( ctmp )

END SUBROUTINE 

