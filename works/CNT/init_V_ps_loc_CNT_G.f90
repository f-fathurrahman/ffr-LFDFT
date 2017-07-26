! Local pseudopotential proposed by Mayer (2004)
!
SUBROUTINE init_V_ps_loc_CNT_G()
  USE m_constants, ONLY : PI, Ry2eV
  USE m_hamiltonian, ONLY : V_ps_loc
  USE m_atoms, ONLY : strf => StructureFactor, Nspecies
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     G2 => LF3d_G2, &
                     NN => LF3d_NN, &
                     LL => LF3d_LL

  IMPLICIT NONE 
  !
  INTEGER, PARAMETER :: Nparams = 3
  REAL(8) :: A(Nparams)
  REAL(8) :: alpha(Nparams)
  !
  INTEGER :: ip, isp, Nx, Ny, Nz, i
  REAL(8) :: Gm, Omega
  COMPLEX(8), ALLOCATABLE :: ctmp(:)

  !
  A(1) = 10.607d0/( Ry2eV*2.d0 )
  A(2) = 29.711d0/( Ry2eV*2.d0 )
  A(3) = -98.911d0/( Ry2eV*2.d0 )
  !
  alpha(1) = 0.12126d0
  alpha(2) = 1.9148d0
  alpha(3) = 0.60078d0

  WRITE(*,*)
  WRITE(*,*) 'V_ps_loc is Gaussian potentials constructed in G-space'

  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  ! Cell volume
  Omega = LL(1) * LL(2) * LL(3)

  ALLOCATE( ctmp(Npoints) )

  V_ps_loc(:) = 0.d0

  isp = 1 ! only one species
  DO ip = 1,Npoints
    ctmp(ip) = cmplx(0.d0,0.d0,kind=8)
    Gm = sqrt(G2(ip))
    DO i = 1,Nparams
      ctmp(ip) = -A(i)*(PI/alpha(i))**1.5d0 * exp(-0.25d0*G2(ip)/alpha(i)) + ctmp(ip)
    ENDDO 
    ctmp(ip) = ctmp(ip) * strf(ip,isp) / Omega
  ENDDO

  ! inverse FFT: G -> R
  CALL fft_fftw3( ctmp, Nx, Ny, Nz, .true. )

  ! XXX: Move this outside isp loop ?
  DO ip = 1,Npoints
    V_ps_loc(ip) = V_ps_loc(ip) + real( ctmp(ip), kind=8 )
  ENDDO 

  WRITE(*,*) 'sum(V_ps_loc) = ', sum(V_ps_loc)

  flush(6)

  DEALLOCATE( ctmp )

END SUBROUTINE 

