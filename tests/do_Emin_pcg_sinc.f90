PROGRAM do_Emin_pcg

  USE m_constants, ONLY : Ry2eV
  USE m_options, ONLY : FREE_NABLA2, I_POISSON_SOLVE
  USE m_PsPot, ONLY : PsPot_Dir, NbetaNL
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs

  IMPLICIT NONE 
  INTEGER :: Narg
  INTEGER :: NN(3)
  REAL(8) :: hh(3)
  CHARACTER(64) :: filexyz, arg_N
  INTEGER :: ip, ist, N_in
  INTEGER :: iargc  ! pgf90 
  INTEGER :: tstart, counts_per_second, tstop

  CALL system_clock( tstart, counts_per_second )

  Narg = iargc()
  IF( Narg /= 2 ) THEN 
    WRITE(*,*) 'ERROR: exactly two arguments must be given:'
    WRITE(*,*) '       N and path to structure file'
    STOP 
  ENDIF 

  CALL getarg( 1, arg_N )
  READ(arg_N, *) N_in

  CALL getarg( 2, filexyz )

  CALL init_atoms_xyz(filexyz)

  ! Override PsPot_Dir
  PsPot_Dir = '../pseudopotentials/pade_gth/'
  CALL init_PsPot()

  !
  NN = (/ N_in, N_in, N_in /)
  hh(:) = (/1.d0, 1.d0, 1.d0/)*(16.d0/(NN(1)-1))
  CALL init_LF3d_sinc( NN, hh )

  CALL info_atoms()
  CALL info_PsPot()
  CALL info_LF3d()

!  CALL init_betaNL()
  NbetaNL = 0

  ! Initialize occupation numbers
  CALL init_states()

  ! Memory for potentials
  CALL alloc_hamiltonian()

  ! Local pseudopotential
  CALL init_V_ps_loc()

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

  I_POISSON_SOLVE = 2

  IF( I_POISSON_SOLVE == 1 ) THEN
    CALL init_Poisson_solve_ISF()
  ELSEIF( I_POISSON_SOLVE == 2 ) THEN
    CALL init_Poisson_solve_DAGE()
  ENDIF

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
  CALL dealloc_PsPot()
  CALL dealloc_atoms()

  CALL system_clock( tstop )

  WRITE(*,*)
  WRITE(*,*) 'Total elapsed time: ', dble(tstop - tstart)/counts_per_second, ' second.'
  WRITE(*,*)

END PROGRAM

