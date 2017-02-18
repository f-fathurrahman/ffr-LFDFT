! efefer, 7 May 2016
!
! Solution of Kohn-Sham equation


!------------------------------------------------------------------------------
PROGRAM t_LF3d
!------------------------------------------------------------------------------
  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : init_LF3d_c
  USE m_globals, ONLY : LF, N, A, B, deltaV, Nbasis
  USE m_globals, ONLY : evecs, evals
  USE m_globals, ONLY : Etot, Ekin, Epot, Ehartree, Exc
  USE m_globals, ONLY : Rho, Vhartree, Vpot, Vxc
  USE m_globals, ONLY : Nstates, Focc
  IMPLICIT NONE
  ! Parameter for harmonic potential
  REAL(8), PARAMETER :: omega=2.d0
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
  Nstates = 4
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

  ! Eigenvectors and eigenvalues
  ALLOCATE( evecs(Nbasis,Nstates) )
  ALLOCATE( evals(Nstates) )

  CALL minimE_sd( 1.d-4, 5000, .TRUE. )

  CALL print_energies()
  !CALL print_evals(Nstates,evals)

  DEALLOCATE( Focc )
  DEALLOCATE( Rho )
  DEALLOCATE( Vhartree )
  DEALLOCATE( evecs )
  DEALLOCATE( evals )
  DEALLOCATE( Vpot )

END PROGRAM



