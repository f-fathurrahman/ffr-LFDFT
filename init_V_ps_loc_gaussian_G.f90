! Gaussian potential:
! V(r) = sum_{i} A_{i} \exp( -\alpha_{i} r^{2} )
!
SUBROUTINE init_V_ps_loc_gaussian_G( Nparams, A, alpha )
  USE m_constants, ONLY : PI
  USE m_hamiltonian, ONLY : V_ps_loc
  USE m_atoms, ONLY : strf => StructureFactor, Nspecies
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     G2 => LF3d_G2, &
                     NN => LF3d_NN, &
                     LL => LF3d_LL

  IMPLICIT NONE 
  !
  INTEGER :: Nparams
  REAL(8) :: A(Nparams)
  REAL(8) :: alpha(Nparams)
  !
  INTEGER :: ip, isp, Nx, Ny, Nz
  REAL(8) :: Gm, Omega
  COMPLEX(8), ALLOCATABLE :: ctmp(:)

  IF( Nspecies /= Nparams ) THEN 
    WRITE(*,*) 'ERROR in constructing Gaussian potential in G-space'
    WRITE(*,'(1x,A,2I4)') 'Nspecies /= Nparams : ', Nspecies, Nparams
    STOP 
  ENDIF 

  WRITE(*,*)
  WRITE(*,*) 'V_ps_loc is Gaussian potentials constructed in G-space'

  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  ! Cell volume
  Omega = LL(1) * LL(2) * LL(3)

  ALLOCATE( ctmp(Npoints) )

  V_ps_loc(:) = 0.d0

  DO isp = 1,Nspecies

    ctmp(:) = cmplx(0.d0,0.d0,kind=8)
    DO ip = 1,Npoints
      Gm = sqrt(G2(ip))
      ctmp(ip) = -A(isp)*(PI/alpha(isp))**1.5d0 * exp(-0.25d0*G2(ip)/alpha(isp)) &
                 * strf(ip,isp) / Omega
    ENDDO

    ! inverse FFT: G -> R
    CALL fft_fftw3( ctmp, Nx, Ny, Nz, .true. )

    ! XXX: Move this outside isp loop ?
    DO ip = 1,Npoints
      V_ps_loc(ip) = V_ps_loc(ip) + real( ctmp(ip), kind=8 )
    ENDDO 

  ENDDO 

  WRITE(*,*) 'sum(V_ps_loc) = ', sum(V_ps_loc)

  flush(6)

  DEALLOCATE( ctmp )

END SUBROUTINE 

