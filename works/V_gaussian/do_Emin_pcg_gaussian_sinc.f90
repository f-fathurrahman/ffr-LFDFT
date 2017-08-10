PROGRAM do_Emin_pcg_gaussian

  USE m_constants, ONLY : Ry2eV
  USE m_options, ONLY : FREE_NABLA2, I_POISSON_SOLVE
  USE m_input_vars, ONLY : assume_isolated
  USE m_atoms, ONLY : Natoms, Nspecies, AtomicCoords, atm2species, &
                      SpeciesSymbols, Zv => AtomicValences
  USE m_states, ONLY : Nstates, Focc, Nelectrons, Nstates_occ
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, dVol => LF3d_dVol
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  USE m_PsPot, ONLY : NbetaNL
  IMPLICIT NONE 
  INTEGER :: Narg
  INTEGER :: NN(3)
  REAL(8) :: hh(3)
  CHARACTER(64) :: filexyz, arg_tmp
  INTEGER :: ip, ist, N_in
  !
  INTEGER :: Nparams
  REAL(8), ALLOCATABLE :: alpha(:), A(:)
  REAL(8) :: alpha_in, A_in
  !
  INTEGER :: iargc  ! pgf90 
  INTEGER :: tstart, counts_per_second, tstop

  Narg = iargc()
  IF( Narg /= 3 ) THEN 
    WRITE(*,*) 'ERROR: exactly three arguments must be given:'
    WRITE(*,*) '       N A alpha'
    WRITE(*,*)
    WRITE(*,*) ' A and alpha is Gaussian parameter:'
    WRITE(*,*) '   f(r) = A*exp( - alpha*r^2 )'
    STOP 
  ENDIF 

  CALL getarg( 1, arg_tmp )
  READ(arg_tmp, *) N_in
  
  CALL getarg( 2, arg_tmp )
  READ(arg_tmp, *) A_in
  
  CALL getarg( 3, arg_tmp )
  READ(arg_tmp, *) alpha_in

  ! Force assume_isolated to 'sinc'
  assume_isolated = 'sinc'

  ! Start timing
  CALL system_clock( tstart, counts_per_second )

  ! Initialize states and occupation numbers MANUALLY
  Nstates = 1
  Nstates_occ = Nstates
  Nelectrons = 2.d0*Nstates
  ALLOCATE( Zv(1) )
  Zv(1) = Nelectrons
  ALLOCATE( Focc(Nstates) )
  Focc(:) = 2.d0

  ! 'Atomic' positions
  Nspecies = 1
  Natoms = 1
  ALLOCATE( AtomicCoords(3,Natoms) )
  AtomicCoords(:,1) = (/ 0.d0, 0.d0, 0.d0 /)
  ALLOCATE( atm2species(Natoms) )
  atm2species(1) = 1
  ALLOCATE( SpeciesSymbols(Nspecies) )
  SpeciesSymbols(1) = 'X'

  !
  NN(:) = N_in
  hh(:) = 16.d0/(N_in-1)
  CALL init_LF3d_sinc( NN, hh )
  CALL info_LF3d()

  ! Memory for potentials
  CALL alloc_hamiltonian()

  ! Local pseudopotential
  Nparams = Nspecies
  ALLOCATE( A(Nparams) )
  ALLOCATE( alpha(Nparams) )
  A(1) = A_in
  alpha(1) = alpha_in
  !
  CALL init_V_ps_loc_gaussian( Nparams, A, alpha ) 

  ! Set this explicitly in order to skip any term involving nonlocal pseudopotentials
  NbetaNL = 0

  CALL calc_E_NN()

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

  ! Manually allocate KS eigenvectors and eigenvalues
  ALLOCATE( evecs(Npoints,Nstates), evals(Nstates) )

  ! Initialize to random wavefunction
  DO ist = 1, Nstates
    DO ip = 1, Npoints
      CALL random_number( evecs(ip,ist) )
    ENDDO
  ENDDO
  CALL orthonormalize( Nstates, evecs )

  CALL ortho_check( Npoints, Nstates, dVol, evecs )

  CALL KS_solve_Emin_pcg( 3.d-5, .FALSE. )

  CALL info_energies()

  CALL calc_evals( Nstates, Focc, evecs, evals )
  WRITE(*,*)
  WRITE(*,*) 'Final eigenvalues (Ha and eV)'
  WRITE(*,*)
  DO ist = 1,Nstates
    WRITE(*,'(1x,I8,2F18.10)') ist, evals(ist), evals(ist)*2.d0*Ry2eV
  ENDDO 

  !
  DEALLOCATE( evecs, evals )
  DEALLOCATE( Focc )

  IF( I_POISSON_SOLVE == 1 ) THEN 
    CALL dealloc_Poisson_solve_ISF()
  ELSEIF( I_POISSON_SOLVE == 2 ) THEN 
    CALL dealloc_Poisson_solve_DAGE()
  ENDIF 

  CALL dealloc_nabla2_sparse()
  CALL dealloc_ilu0_prec()
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()
  CALL dealloc_atoms()

  CALL system_clock( tstop )

  WRITE(*,*)
  WRITE(*,'(1x,A,ES18.10,A)') 'Total elapsed time: ', &
           dble(tstop - tstart)/counts_per_second, ' second.'
  WRITE(*,*)

END PROGRAM

