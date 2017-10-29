!!> Some initial steps
SUBROUTINE setup_ffr_LFDFT()

  USE m_options, ONLY : FREE_NABLA2, I_POISSON_SOLVE
  USE m_input_vars, ONLY : assume_isolated
  ! 
  INTEGER :: Narg
  INTEGER :: iargc  ! needed for several compilers
  CHARACTER(64) :: filein

  Narg = iargc()
  IF( Narg /= 1 ) THEN
    WRITE(*,*)
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

  ! FIXME: needed anyway even for sinc ??
  ! Structure factor, shifted to FFT grid
  CALL init_strfact_shifted()

  ! Ewald energy
  IF( assume_isolated == 'sinc' ) THEN
    CALL calc_E_NN()
  ELSE
    CALL calc_Ewald_qe()
  ENDIF

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
END SUBROUTINE 
