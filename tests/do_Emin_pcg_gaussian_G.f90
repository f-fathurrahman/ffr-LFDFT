PROGRAM do_Emin_pcg_gaussian_G

  USE m_constants, ONLY : Ry2eV
  USE m_options, ONLY : FREE_NABLA2
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
  REAL(8) :: AA(3), BB(3)
  CHARACTER(64) :: filexyz, arg_N
  INTEGER :: ip, ist, N_in
  !
  INTEGER :: Nparams
  REAL(8), ALLOCATABLE :: alpha(:), A(:)
  !
  INTEGER :: iargc  ! pgf90 
  INTEGER :: tstart, counts_per_second, tstop

  CALL system_clock( tstart, counts_per_second )

  Narg = iargc()
  IF( Narg /= 1 ) THEN 
    WRITE(*,*) 'ERROR: exactly one arguments must be given:'
    WRITE(*,*) '       N'
    STOP 
  ENDIF 

  CALL getarg( 1, arg_N )
  READ(arg_N, *) N_in

  ! Initialize states and occupation numbers MANUALLY
  Nstates = 1
  Nstates_occ = 1
  Nelectrons = 2.d0
  ALLOCATE( Zv(1) )
  Zv(1) = Nelectrons
  ALLOCATE( Focc(Nstates) )
  Focc(1) = 2.d0

  ! 'Atomic' positions
  Nspecies = 1
  Natoms = 1
  ALLOCATE( AtomicCoords(3,Natoms) )
  AtomicCoords(:,1) = 0.d0
  ALLOCATE( atm2species(Natoms) )
  atm2species(1) = 1
  ALLOCATE( SpeciesSymbols(Nspecies) )
  SpeciesSymbols(1) = 'X'

  !
  NN = (/ N_in, N_in, N_in /)
  AA = (/ -8.d0, -8.d0, -8.d0 /)
  BB = (/ 8.d0, 8.d0, 8.d0 /)
  CALL init_LF3d_p( NN, AA, BB )
  CALL info_LF3d()

  ! Structure factor, shifted to FFT grid
  CALL init_strfact_shifted()

  ! Memory for potentials
  CALL alloc_hamiltonian()

  ! Local pseudopotential
  Nparams = Nspecies
  ALLOCATE( A(Nparams) )
  ALLOCATE( alpha(Nparams) )
  A(1) = 1.d0
  alpha(1) = 3.d0
  !
  CALL init_V_ps_loc_gaussian_G( Nparams, A, alpha )

  ! 
  NbetaNL = 0
  CALL calc_Ewald_qe()

  ! Laplacian matrix
  CALL init_nabla2_sparse()

  ! ILU0 preconditioner based on kinetic matrix
  CALL init_ilu0_prec()

  IF( FREE_NABLA2 ) THEN 
    CALL dealloc_nabla2_sparse()
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

  CALL dealloc_nabla2_sparse()
  CALL dealloc_ilu0_prec()
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()
!  CALL dealloc_atoms()

  CALL system_clock( tstop )

  WRITE(*,*)
  WRITE(*,*) 'Total elapsed time: ', dble(tstop - tstart)/counts_per_second, ' second.'
  WRITE(*,*)

END PROGRAM

