! efefer, 15 March 2016
!
! Solution of Schrodinger equation
!
! Using iterative (partial) diagonalization


!------------------------------------------------------------------------------
PROGRAM LFdft
!------------------------------------------------------------------------------
  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : init_LF3d_sinc, init_LF3d_c, init_LF3d_p
  USE m_globals, ONLY : LF, Npoints, Vpsloc, Rho, A,B,h,N, LF_type, Solution_Method, &
     evecs, evals, Nstate, Focc, Vhartree
  IMPLICIT NONE
  !
  !REAL(8) :: Ekin, Epot, Etot

  CALL read_arguments()

  CALL init_system()

  WRITE(*,'(/,1x,A)') 'Initializing grids and basis functions:'
  WRITE(*,*)          '---------------------------------------'
  IF( LF_type == 'sinc' ) THEN
    CALL init_LF3d_sinc(LF, (/N,N,N/), (/h,h,h/) )
    A = LF%LFx%A
    B = LF%LFx%B
  ELSEIF( LF_type == 'box' ) THEN
    CALL init_LF3d_c( LF, (/N,N,N/), (/A,A,A/), (/B,B,B/) )
  ELSEIF( LF_type == 'per' ) THEN
    ! TODO: This is not yet supported
    CALL init_LF3d_p( LF, (/N,N,N/), (/A,A,A/), (/B,B,B/) )
  ENDIF
  
  ! FIXME
  Npoints = N**3


  ! Set up potential
  ALLOCATE( Vpsloc(Npoints) )
  ALLOCATE( Rho(Npoints) )
  ALLOCATE( evecs(Npoints,Nstate), evals(Nstate) )
  ALLOCATE( Vhartree(Npoints) )

  IF( LF_type == 'per' ) THEN
    CALL init_rho_gaussian_per()
  ELSE
    CALL init_rho_gaussian()
  ENDIF

  CALL test_write_xsf()
  STOP 

  CALL calc_hartree()


  CALL init_Vpsloc()


  IF( Solution_Method == 'diag') THEN
    CALL solve_diagonalize()
  ELSEIF( Solution_Method == 'cg' .OR. Solution_Method == 'pcg' ) THEN
    CALL solve_minim()
  ELSEIF( Solution_Method == 'MG' ) THEN
    CALL solve_MG() ! FIXME not implemented yet
  ELSE
    WRITE(*,*) 'ERROR: Unrecognized method: ', trim(Solution_Method)
    STOP
  ENDIF

  CALL calc_rho()

  DEALLOCATE( Focc )
  DEALLOCATE( Rho )
  DEALLOCATE( Vpsloc )
  DEALLOCATE( evecs, evals )

END PROGRAM



