PROGRAM ffr_LFDFT

  USE m_constants, ONLY : Ry2eV
  USE m_input_vars, ONLY : startingwfc, assume_isolated
  USE m_options, ONLY : FREE_NABLA2, I_KS_Solve, I_POISSON_SOLVE
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs

  IMPLICIT NONE 
  INTEGER :: Narg
  CHARACTER(64) :: filein
  INTEGER :: ist
  INTEGER :: iargc  ! pgf90 
  INTEGER :: tstart, counts_per_second, tstop

  Narg = iargc()
  IF( Narg /= 1 ) THEN 
    WRITE(*,*) 'ERROR: exactly one arguments must be given: input file path'
    STOP 
  ENDIF 

  CALL getarg( 1, filein )

  CALL system_clock( tstart, counts_per_second )

  CALL welcome()

  CALL read_input( filein )
  CALL setup_from_input()
  CALL setup_options()

  CALL info_atoms()
  CALL info_PsPot()
  CALL info_LF3d()
 
  ! needs LF3d to be initialized first
  CALL init_betaNL()

  ! Initialize occupation numbers
  CALL init_states()

  ! FIXME: needed anyway even for sinc ??
  ! Structure factor, shifted to FFT grid
  CALL init_strfact_shifted()

  ! Ewald energy
  IF( assume_isolated /= 'sinc' ) THEN 
    !
    CALL calc_Ewald_qe()
  ENDIF 
  ! FIXME: Need subroutine to calculation ion-ion energy

  ! Memory for potentials
  CALL alloc_hamiltonian()

  ! Local pseudopotential
  IF( assume_isolated == 'sinc' ) THEN 
    CALL init_V_ps_loc()
  ELSE 
    CALL init_V_ps_loc_G()
  ENDIF 

  ! Laplacian matrix
  CALL init_nabla2_sparse()
  ! ILU0 preconditioner based on kinetic matrix
  CALL init_ilu0_prec()

  IF( FREE_NABLA2 ) THEN 
    CALL dealloc_nabla2_sparse()
  ENDIF 

  IF( I_POISSON_SOLVE == 1 ) THEN 
    CALL init_Poisson_solve_ISF()
  ELSEIF( I_POISSON_SOLVE == 2 ) THEN 
    CALL init_Poisson_solve_DAGE()
  ENDIF 

  ! Guess density
  IF( startingwfc /= 'random' ) THEN 
    CALL gen_guess_rho_gaussian()
  ENDIF 

  ! Manually allocate KS eigenvectors and eigenvalues
  ALLOCATE( evecs(Npoints,Nstates), evals(Nstates) )

  CALL gen_random_evecs()  ! also needed for initial diagonalization routine

  IF( startingwfc /= 'random' ) THEN 
    ! This will call diagonalization routine
    CALL gen_gaussian_evecs()
  ELSE
    WRITE(*,*)
    WRITE(*,*) 'Using random starting wavefunction'
  ENDIF 

  IF( I_KS_SOLVE == 1 ) THEN 

    CALL KS_solve_Emin_pcg( 3.d-5, .FALSE. )
    CALL info_energies()
    CALL calc_evals( Nstates, Focc, evecs, evals )

  ELSEIF( I_KS_SOLVE == 2 ) THEN 

    IF( startingwfc == 'random' ) THEN 
      ! Initial Rhoe and potentials
      CALL calc_rhoe( Focc, evecs )
      CALL update_potentials()
    ENDIF 

    CALL KS_solve_SCF()
    CALL info_energies()

  ENDIF 
  
  WRITE(*,*)
  WRITE(*,*) 'Final eigenvalues (Ha and eV)'
  WRITE(*,*)
  DO ist = 1,Nstates
    WRITE(*,'(1x,I8,2F18.10)') ist, evals(ist), evals(ist)*2.d0*Ry2eV
  ENDDO

  !
  DEALLOCATE( evecs, evals )
  DEALLOCATE( Focc )

  CALL dealloc_nabla2_sparse()
  CALL dealloc_ilu0_prec()
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()
  CALL dealloc_PsPot()
  CALL dealloc_atoms()

  CALL system_clock( tstop )

  WRITE(*,*)
  WRITE(*,'(1x,A,ES18.10,A)') 'Total elapsed time: ', &
           dble(tstop - tstart)/counts_per_second, ' second.'

  CALL goodbye()

END PROGRAM

