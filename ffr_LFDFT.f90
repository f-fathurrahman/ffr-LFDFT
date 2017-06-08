PROGRAM ffr_LFDFT

  USE m_constants, ONLY : Ry2eV
  USE m_options, ONLY : FREE_NABLA2, I_KS_Solve
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

  CALL system_clock( tstart, counts_per_second )

  Narg = iargc()
  IF( Narg /= 1 ) THEN 
    WRITE(*,*) 'ERROR: exactly one arguments must be given: input file path'
    STOP 
  ENDIF 

  CALL getarg( 1, filein )

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

  ! Structure factor, shifted to FFT grid
  CALL init_strfact_shifted()

  ! Ewald energy
  CALL calc_Ewald_qe()

  ! Memory for potentials
  CALL alloc_hamiltonian()

  ! Local pseudopotential
  CALL init_V_ps_loc_G()

  ! Laplacian matrix
  CALL init_nabla2_sparse()
  ! ILU0 preconditioner based on kinetic matrix
  CALL init_ilu0_prec()

  IF( FREE_NABLA2 ) THEN 
    CALL dealloc_nabla2_sparse()
  ENDIF 

  ! Guess density
  CALL gen_guess_rho_gaussian()

  ! Manually allocate KS eigenvectors and eigenvalues
  ALLOCATE( evecs(Npoints,Nstates), evals(Nstates) )

  CALL gen_random_evecs()
  CALL gen_gaussian_evecs()

  IF( I_KS_SOLVE == 1 ) THEN 

    CALL KS_solve_Emin_pcg( 3.d-5, .FALSE. )
    !CALL KS_solve_Emin_pcg( 3.d-5, 1000, .TRUE. )
    CALL info_energies()
    CALL calc_evals( Nstates, Focc, evecs, evals )
    WRITE(*,*)
    WRITE(*,*) 'Final eigenvalues (Ha and eV)'
    WRITE(*,*)
    DO ist = 1,Nstates
      WRITE(*,'(1x,I8,2F18.10)') ist, evals(ist), evals(ist)*2.d0*Ry2eV
    ENDDO

  ELSEIF( I_KS_SOLVE == 2 ) THEN 
    ! Initial Rhoe and potentials
    CALL calc_rhoe( Focc, evecs )
    CALL update_potentials()
    CALL KS_solve_SCF()
    !
    CALL info_energies()

  ENDIF 

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
  WRITE(*,*) 'Total elapsed time: ', dble(tstop - tstart)/counts_per_second, ' seconds.'
  WRITE(*,*)

END PROGRAM

