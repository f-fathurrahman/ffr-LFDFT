! efefer, 7 May 2016
!
! Solution of Kohn-Sham equation


!------------------------------------------------------------------------------
PROGRAM t_LF3d
!------------------------------------------------------------------------------
  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : init_LF3d_c
  USE m_globals, ONLY : LF, N, A, B, deltaV
  USE m_globals, ONLY : evecs, evals
  USE m_globals, ONLY : Etot, Ekin, Epot, Ehartree, Exc
  USE m_globals, ONLY : Rho, Vhartree, Vpot, Vxc
  USE m_globals, ONLY : Nstates, Focc
  IMPLICIT NONE
  ! Parameter for harmonic potential
  REAL(8), PARAMETER :: omega=2.d0
  INTEGER :: Nbasis
  REAL(8) :: sigma1, integRho, dr
  INTEGER :: ip
  REAL(8), ALLOCATABLE :: Rho_old(:)
  REAL(8) :: dEtot, Etot_old
  REAL(8) :: beta0, dmix, betamax
  INTEGER :: MaxSCFIters, iterSCF
  REAL(8), ALLOCATABLE :: epsxc(:), depsxc(:)
  REAL(8), ALLOCATABLE :: fmix(:), betamix(:)
  INTEGER :: DIM_PULAY
  REAL(8), ALLOCATABLE :: fmix_pulay(:,:), mu_pulay(:,:)
  REAL(8) :: r0(3)
  REAL(8) :: Nelectrons

  N = 16
  A = 0.d0
  B = 6.d0

  r0 = A + (B-A)/2d0

  CALL init_LF3d_c(LF, (/N,N,N/), (/A,A,A/), (/B,B,B/) )

  deltaV = LF%LFx%h * LF%LFy%h * LF%LFz%h
  WRITE(*,*) 'deltaV = ', deltaV

  ! For convenience
  Nbasis  = N**3
  !
  Nstates = 1
  ALLOCATE( Focc(Nstates) )
  Focc(:) = 2.d0
  Nelectrons = sum(Focc)

  ! Set up potential
  ALLOCATE( Vpot(Nbasis) )
  CALL init_pot_harmonic( omega, Vpot )
  !CALL init_pspot_H( Vpot )

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

  ! Build initial charge density
  ! any N-representable charge density should work
  sigma1 = 0.3d0 ! tune this such that delta Integrated rho is small
  DO ip = 1, Nbasis
    dr = norm2( LF%lingrid(:,ip) - r0(:) )
    Rho(ip) = Nelectrons*exp( -dr**2/(2.d0*sigma1**2) )/(2.d0*PI*sigma1**2)**1.5d0
  ENDDO
  integRho = sum(Rho(:))*deltaV
  Rho = Rho/integRho*Nelectrons  ! renormalize
  integRho = sum(Rho(:))*deltaV
  WRITE(*,*) 'delta Integrated Rho: ', abs(Nelectrons-integRho)
  !STOP


  Etot_old = 0.d0 ! an unlikely number for total energy (?)
  MaxSCFIters = 150
  beta0 = 0.001d0
  betamax = 0.5d0
  Rho_old(:) = Rho(:)

  ! for adaptive mixing
  ALLOCATE( fmix(Nbasis) )
  ALLOCATE( betamix(Nbasis) )

  ! for Pulay mixing
  !DIM_PULAY = 4
  !ALLOCATE( fmix_pulay(Nbasis,DIM_PULAY) ); fmix_pulay = 0.d0
  !ALLOCATE( mu_pulay(Nbasis,DIM_PULAY) ); mu_pulay = 0.d0

  WRITE(*,'(1x,A,F18.10)') 'beta0 = ', beta0

  DO iterSCF = 0, MaxSCFIters
    !
    CALL solve_poisson_cg( Nbasis, -4.d0*PI*Rho(:), Vhartree )
    Ehartree = 0.5d0*sum( Rho(:) * Vhartree(:) ) * deltaV
    !Vhartree(:) = 0.d0
    !Ehartree    = 0.d0
    !
    CALL excVWN( Nbasis, Rho, epsxc )
    Exc      = sum( Rho(:) * epsxc(:) ) * deltaV
    !
    CALL excpVWN( Nbasis, Rho, depsxc )
    Vxc      = epsxc(:) + Rho(:) * depsxc(:)
    !Exc    = 0.d0
    !Vxc(:) = 0.d0

    CALL diag_Ham( Nbasis, Nstates, evecs, evals )
    !CALL print_evals( Nstates, evals )

    CALL get_Etot( Nbasis, Nstates, evecs, Ekin, Epot, Etot )
    CALL print_energies()

    CALL get_Rho( Nbasis, Nstates, Focc, evecs, Rho )

    dEtot = abs(Etot - Etot_old)
    IF( dEtot < 1.d-6 ) THEN
      WRITE(*,*) 'Convergence achieved: Etot, dEtot', Etot, dEtot
      EXIT
    ENDIF

    WRITE(*,*)
    WRITE(*,'(1x,A,I5,3F18.10)') 'SCF: ', iterSCF, Etot, dEtot, norm2(Rho-Rho_old)

    ! new rho
    !Rho = beta0*Rho + (1.d0-beta0)*Rho_old
    !integRho = sum(Rho(:))*deltaV
    !Rho = Rho/integRho*2.d0*Nstates
    !integRho = sum(Rho(:))*deltaV
    ! 
    ! Choose of the the following charge density mixing
    CALL mixlinear( iterSCF, beta0, Nbasis, Rho, Rho_old, dmix )
    !CALL mixadapt( iterSCF, beta0, betamax, Nbasis, Rho, Rho_old, betamix, fmix, dmix )
    !CALL mixpulay( iterSCF, Nbasis, DIM_PULAY, Rho, mu_pulay, fmix_pulay, dmix )
    ! Renormalize charge density
    integRho = sum(Rho(:))*deltaV
    Rho = Rho/integRho*Nelectrons
    integRho = sum(Rho(:))*deltaV
    WRITE(*,*) 'MIXING: delta Integrated Rho: ', abs(Nelectrons-integRho)
    WRITE(*,*) 'MIXING: dmix = ', dmix
    !STOP

    Rho_old(:) = Rho(:)

    Etot_old = Etot

  ENDDO

  CALL print_energies()
  CALL print_evals(Nstates,evals)

  !DEALLOCATE( mu_pulay )
  !DEALLOCATE( fmix_pulay )
  DEALLOCATE( Focc )
  DEALLOCATE( Rho )
  DEALLOCATE( Vhartree )
  DEALLOCATE( Rho_old )
  DEALLOCATE( evecs )
  DEALLOCATE( evals )
  DEALLOCATE( Vpot )
  DEALLOCATE( epsxc )
  DEALLOCATE( depsxc )

END PROGRAM



