PROGRAM test_sch

  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, dVol => LF3d_dVol
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  USE m_atoms, ONLY : Nspecies, atpos => AtomicCoords, Natoms, atm2species
  USE m_PsPot, ONLY : NbetaNL
  IMPLICIT NONE
  !
  INTEGER :: ist, ip
  INTEGER :: NN(3)
  REAL(8) :: hh(3)
  REAL(8) :: ddot
  !
  REAL(8), ALLOCATABLE :: A(:), alpha(:)
  REAL(8) :: A_in, alpha_in
  INTEGER :: Nparams
  INTEGER :: N_in

  N_in = 35
  NN(:) = N_in
  hh(:) = 16.d0/(N_in-1)
  CALL init_LF3d_sinc( NN, hh )

  CALL info_LF3d()

  ! Set up potential
  CALL alloc_hamiltonian()

  CALL init_nabla2_sparse()
  CALL init_ilu0_prec()

  ! Local pseudopotential
  alpha_in = 2.5d0
  A_in = 1.d0
  ! A = 10, alpha=0.45 --> small band gap
  ! band gap is proportional to sqrt(A_in)/alpha_in
  Nspecies = 1
  Nparams = Nspecies
  ALLOCATE( A(Nparams) )
  ALLOCATE( alpha(Nparams) )
  !
  A(1) = A_in/(2.d0*PI*alpha_in**2)**1.5d0
  alpha(1) = 0.5d0/alpha_in**2
  !
  Natoms = 1
  ALLOCATE( atpos(3,Natoms) )
  ALLOCATE( atm2species(Natoms) )
  atpos(:,1) = (/ 0.d0, 0.d0, 0.d0 /)
  atm2species(1) = 1
  !
  CALL init_V_ps_loc_gaussian( Nparams, A, alpha )

  ! States initialization
  Nstates = 4
  ALLOCATE( Focc(Nstates) )
  Focc(:) = 1.d0

  NbetaNL = 0

  ! Set up initial guess
  ALLOCATE( evecs(Npoints,Nstates), evals(Nstates) )
  DO ist = 1, Nstates
    DO ip = 1, Npoints
      CALL random_number( evecs(ip,ist) )
    ENDDO
  ENDDO
  CALL orthonormalize( Nstates, evecs )

  CALL Sch_solve_Emin_pcg( 2, 3.d-5, .FALSE., 1d-4, 100, .TRUE. )

  WRITE(*,*)
  WRITE(*,*) 'Final eigenvalues:'
  DO ist = 1, Nstates
    WRITE(*,'(1x,I4,F18.10)') ist, evals(ist)
  ENDDO

  WRITE(*,*)
  WRITE(*,*) 'Check normalization:'
  DO ist = 1, Nstates
    WRITE(*,'(1x,I4,F18.10)') ist, ddot( Npoints, evecs(:,ist), 1, evecs(:,ist), 1 )*dVol
  ENDDO 

  DEALLOCATE( evecs, evals )
  DEALLOCATE( Focc )
  CALL dealloc_atoms()
  CALL dealloc_nabla2_sparse()
  CALL dealloc_ilu0_prec()
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()

END PROGRAM


