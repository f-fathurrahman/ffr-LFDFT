! efefer, 7 May 2016
!
! Solution of Kohn-Sham equation


!------------------------------------------------------------------------------
PROGRAM t_LF3d_c
!------------------------------------------------------------------------------
  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : init_LF3d_sinc
  USE m_globals, ONLY : LF, N, A, B, deltaV
  USE m_globals, ONLY : evecs, evals
  USE m_globals, ONLY : Etot, Ekin, Epot, Ehartree, Exc
  USE m_globals, ONLY : Rho, Vhartree, Vpot, Vxc
  IMPLICIT NONE
  ! Parameter for harmonic potential
  REAL(8), PARAMETER :: omega=2.d0
  INTEGER :: Nstates, Nbasis
  REAL(8) :: sigma1, integRho, dr
  INTEGER :: ip
  REAL(8), ALLOCATABLE :: Rho_old(:)
  REAL(8) :: dEtot, Etot_old
  REAL(8) :: betaMix
  INTEGER :: MaxSCFIters, iterSCF
  REAL(8), ALLOCATABLE :: epsxc(:), depsxc(:)
  REAL(8) :: h

  N = 24
  h = 0.3d0

  CALL init_LF3d_sinc(LF, (/N,N,N/), (/h,h,h/) )
  ! FIXME Need these?
  A = LF%LFx%A
  B = LF%LFx%B

  deltaV = LF%LFx%h * LF%LFy%h * LF%LFz%h
  WRITE(*,*) 'deltaV = ', deltaV

  Nbasis  = N**3
  Nstates = 4

  ! Set up potential
  ALLOCATE( Vpot(Nbasis) )
  CALL init_pot_harmonic( omega, Vpot )

  ! Initial rho, Hartree and XC potential
  ALLOCATE( Rho(Nbasis) )
  ALLOCATE( Vhartree(Nbasis) )
  ALLOCATE( Vxc(Nbasis) )

  ALLOCATE( epsxc(Nbasis) )
  ALLOCATE( depsxc(Nbasis) )

  ALLOCATE( Rho_old(Nbasis) )

  ! Eigenvectors and eigenvalues
  ALLOCATE( evecs(Nbasis,Nstates) )
  ALLOCATE( evals(Nstates) )

  !
  sigma1 = 0.4d0 ! tune this such that delta Integrated rho is small
  DO ip = 1, Nbasis
    dr = norm2( LF.lingrid(:,ip) )
    Rho(ip) = 2.d0*Nstates*exp( -dr**2/(2.d0*sigma1**2) )/(2.d0*PI*sigma1**2)**1.5d0
    !WRITE(*,*) ip, dr, Rho(ip)
  ENDDO
  integRho = sum(Rho(:))*deltaV
  WRITE(*,*) 'delta Integrated Rho: ', abs(dble(2.d0*Nstates)-integRho)
  !STOP


  Etot_old = 0.d0 ! an unlikely number for total energy (?)
  MaxSCFIters = 50
  betaMix = 0.7d0
  Rho_old(:) = Rho(:)

  WRITE(*,'(1x,A,F18.10)') 'betaMix = ', betaMix

  DO iterSCF = 1, MaxSCFIters
    !
    CALL solve_poisson_cg( Nbasis, -4.d0*PI*Rho(:), Vhartree )
    Ehartree = 0.5d0*sum( Rho(:) * Vhartree(:) ) * deltaV
    !
    CALL excVWN( Nbasis, Rho, epsxc )
    Exc      = sum( Rho(:) * epsxc(:) ) * deltaV
    !
    CALL excpVWN( Nbasis, Rho, depsxc )
    Vxc      = epsxc(:) + Rho(:) * depsxc(:)

    CALL diag_Ham( Nbasis, Nstates, evecs, evals )
    !CALL print_evals( Nstates, evals )

    CALL get_Etot( Nstates, evecs, Ekin, Epot, Etot )
    CALL print_energies()

    CALL get_Rho( Nbasis, Nstates, evecs, Rho )

    dEtot = abs(Etot - Etot_old)
    IF( dEtot < 1.d-6 ) THEN
      WRITE(*,*) 'Convergence achieved'
      EXIT
    ENDIF

    WRITE(*,*)
    WRITE(*,'(1x,A,I5,F18.10,F18.10)') 'SCF: ', iterSCF, Etot, dEtot

    ! new rho
    Rho(:) = betaMix*Rho(:) + (1.d0-betaMix)*Rho_old(:)
    integRho = sum(Rho(:))*deltaV
    WRITE(*,*) 'MIXING: delta Integrated Rho: ', abs(dble(2.d0*Nstates)-integRho)

    Rho_old(:) = Rho(:)
    Etot_old = Etot

  ENDDO

  CALL print_energies()
  CALL print_evals(Nstates,evals)

  DEALLOCATE( Rho )
  DEALLOCATE( Vhartree )
  DEALLOCATE( Rho_old )
  DEALLOCATE( evecs )
  DEALLOCATE( evals )
  DEALLOCATE( Vpot )
  DEALLOCATE( epsxc )
  DEALLOCATE( depsxc )

END PROGRAM



