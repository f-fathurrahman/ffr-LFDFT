PROGRAM test_V_NL

  USE m_PsPot, ONLY : PsPot_Dir
  USE m_states, ONLY : Nstates, &
                       evecs => KS_evecs, &
                       evals => KS_evals
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  IMPLICIT NONE 
  INTEGER :: Narg
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  CHARACTER(64) :: filexyz, arg_N
  INTEGER :: N_in, ip, ist
  INTEGER :: iargc  ! pgf90 

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
  PsPot_Dir = '../HGH/'
  CALL init_PsPot()

  WRITE(*,*) 'Pass here 41'

  !
  NN = (/ N_in, N_in, N_in /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)
  CALL init_LF3d_p( NN, AA, BB )

  CALL info_atoms()
  CALL info_PsPot()
  CALL info_LF3d()

  CALL init_betaNL()

  ! Initialize occupation numbers
  CALL init_states()

  ! Structure factor, shifted to FFT grid
  CALL init_strfact_shifted()

  ! Ewald energy
  CALL calc_Ewald()

  ! Memory for potentials
  CALL alloc_hamiltonian()

  ! Local pseudopotential
  CALL init_V_ps_loc_G()

  ! Laplacian matrix
  CALL init_nabla2_sparse()
  ! ILU0 preconditioner based on kinetic matrix
  CALL init_ilu0_prec()

  ! Manually allocate KS eigenvectors and eigenvalues
  ALLOCATE( evecs(Npoints,Nstates), evals(Nstates) )

  ! Initialize to random wavefunction
  DO ist = 1, Nstates
    DO ip = 1, Npoints
      CALL random_number( evecs(ip,ist) )
    ENDDO
  ENDDO
  CALL orthonormalize( Nstates, evecs )

  CALL KS_solve_Emin_pcg( 3.d-5, 1000, .FALSE. )

  CALL dealloc_atoms()
  CALL dealloc_PsPot()
  CALL dealloc_LF3d()
END PROGRAM 

