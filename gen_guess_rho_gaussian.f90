SUBROUTINE gen_guess_rho_gaussian()
  
  USE m_hamiltonian, ONLY : Rhoe
  USE m_atoms, ONLY : strf => StructureFactor, &
                      Nspecies, AtomicValences, SpeciesSymbols
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     G2 => LF3d_G2, &
                     NN => LF3d_NN, &
                     LL => LF3d_LL, &
                     dVol => LF3d_dVol
  USE m_states, ONLY : Nelectrons
  USE m_constants, ONLY : EPS_SMALL
  IMPLICIT NONE 
  INTEGER :: ip, isp, Nx, Ny, Nz, Ngvec, ig
  REAL(8) :: Omega, ff
  REAL(8) :: length, znucl, zion, integRho
  COMPLEX(8), ALLOCATABLE :: ctmp(:)
  REAL(8) :: atom_znucl

  Ngvec = Npoints

  WRITE(*,*)
  WRITE(*,*) 'Generating guess density'

  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  ! Cell volume
  Omega = LL(1) * LL(2) * LL(3)

  ALLOCATE( ctmp(Npoints) )

  Rhoe(:) = EPS_SMALL

  DO isp = 1,Nspecies
    
    zion = AtomicValences(isp)
    znucl = atom_znucl(SpeciesSymbols(isp))
    CALL atmlength( 0.d0, length, zion, znucl )

    WRITE(*,'(1x,A,A,3F7.3)') 'zion, znucl, length = ', &
               trim(SpeciesSymbols(isp)), zion, znucl, length

    ctmp(:) = cmplx(0.d0,0.d0,kind=8)
    DO ig = 1,Ngvec
      ctmp(ig) = exp(-length**2*G2(ig)) * strf(ig,isp) / Omega
    ENDDO

    ! inverse FFT: G -> R
    CALL fft_fftw3( ctmp, Nx, Ny, Nz, .true. )

    ! XXX: Move this outside isp loop ?
    DO ip = 1,Npoints
      ff = real( ctmp(ip), kind=8 )
      IF( ff > 0.d0 ) THEN 
        Rhoe(ip) = Rhoe(ip) + ff
      ELSE
        Rhoe(ip) = Rhoe(ip) + abs(ff)
      ENDIF 
    ENDDO 

  ENDDO 

  integRho = sum(Rhoe)*dVol
  WRITE(*,*) 'Initial: integRho = ', integRho
  ! scale
  Rhoe(:) = Rhoe(:)/integRho*Nelectrons
  integRho = sum(Rhoe)*dVol
  WRITE(*,*) 'After scaling: integRho = ', integRho

  DEALLOCATE( ctmp )

  flush(6)

END SUBROUTINE 

