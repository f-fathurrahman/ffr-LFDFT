! efefer, 15 January 2016
!
! Solution of Schrodinger equation

! Using iterative (partial) diagonalization


!------------------------------------------------------------------------------
PROGRAM t_LF3d_sch
!------------------------------------------------------------------------------
  USE m_constants, ONLY: PI
  USE m_LF3d
  USE m_globals
  IMPLICIT NONE
  !
  INTEGER :: iscf
  REAL(8) :: Ekin, Epot, Etot, dEtot, Etot_old
  REAL(8), ALLOCATABLE :: Rho_old(:)

  CALL read_arguments()

  IF( LF_type == 'sinc' ) THEN
    CALL init_LF3d_sinc(LF, (/N,N,N/), (/h,h,h/) )
    A = LF%LFx%A
    B = LF%LFx%B
  ELSEIF( LF_type == 'box' ) THEN
    CALL init_LF3d_c( LF, (/N,N,N/), (/A,A,A/), (/B,B,B/) )
  ELSEIF( LF_type == 'per' ) THEN
    CALL init_LF3d_p( LF, (/N,N,N/), (/A,A,A/), (/B,B,B/) )
  ENDIF

  ! Set up potential
  ALLOCATE( Vpot(N**3) )
  ALLOCATE( Rho(N**3) )
  ALLOCATE( Rho_old(N**3) )
  ALLOCATE( Vhartree(N**3) )

  Nstate = 1
  Nelec = 1.d0
  ALLOCATE( Focc(Nstate) )
  Focc(1) = 1.d0
  ALLOCATE( evecs(N**3,Nstate), evals(Nstate) )


  !CALL init_pot_harmonic( 2.d0, Vpot )
  !CALL init_pspot_H( Vpot )
  !CALL init_pspot_H2( Vpot )
  !CALL init_pot_softCoulomb( 1.d0, Vpot )
  CALL init_pot_coulomb(1.d0, Vpot)

  CALL init_rho_gaussian()
  STOP

  CALL calc_hartree()

!  IF( Solution_Method == 'diag') THEN
!    CALL solve_diagonalize()
!  ELSEIF( Solution_Method == 'cg' .OR. Solution_Method == 'pcg' ) THEN
!    CALL solve_minim()
!  ELSEIF( Solution_Method == 'MG' ) THEN
!    CALL solve_MG() ! not implemented yet
!  ELSE
!    WRITE(*,*) 'Unrecognized method: ', trim(Solution_Method)
!    STOP
!  ENDIF

  Etot_old = 0.d0
  Rho_old = Rho

  DO iscf = 1, 20

    WRITE(*,*)
    WRITE(*,'(1x,A,I8)') 'SCF cycle #', iscf
    CALL solve_diagonalize()
    CALL get_Etot( Nstate, evecs, Ekin, Epot, Etot )

    dEtot = abs(Etot - Etot_old)
    WRITE(*,'(1x,A,I5,2F18.10)') 'iscf, Etot, dEtot = ', iscf, Etot, dEtot
    Etot_old = Etot

    IF( dEtot < 1.d-6 ) THEN
      WRITE(*,*) 'Convergence achieved'
      EXIT
    ENDIF

    CALL calc_rho()
    Rho = 0.9*Rho + 0.1*Rho_old  ! linear mix
    WRITE(*,*) 'Check Rho after mix:', sum(Rho)

    CALL calc_hartree()
    Rho_old = Rho

  ENDDO


  DEALLOCATE( Focc )
  DEALLOCATE( Rho )
  DEALLOCATE( Vhartree )
  DEALLOCATE( Vpot )
  DEALLOCATE( evecs, evals )

END PROGRAM



