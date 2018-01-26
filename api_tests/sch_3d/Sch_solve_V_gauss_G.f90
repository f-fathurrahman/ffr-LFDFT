PROGRAM test_Sch_solve

  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, dVol => LF3d_dVol
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  USE m_atoms, ONLY : Nspecies, atpos => AtomicCoords, Natoms, atm2species
  USE m_PsPot, ONLY : NbetaNL
  USE m_options, ONLY : I_CG_BETA
  IMPLICIT NONE
  !
  INTEGER :: ist, ip
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  !
  INTEGER, ALLOCATABLE :: btype(:)
  INTEGER :: notcnv, dav_iter
  !
  INTEGER, EXTERNAL :: iargc
  INTEGER :: Nargs
  !
  REAL(8), ALLOCATABLE :: A(:), alpha(:)
  REAL(8) :: A_in, alpha_in
  INTEGER :: Nparams
  INTEGER :: N_in
  CHARACTER(64) :: diag_method
  !
  REAL(8), EXTERNAL :: ddot

  I_CG_BETA = 1 ! Fletcher-Reeves

  CALL setup_args()
  
  NN(:) = N_in
  AA(:) = (/ 0.d0, 0.d0, 0.d0 /)
  BB(:) = (/ 16.d0, 16.d0, 16.d0 /)
  CALL init_LF3d_p( NN, AA, BB )

  ! Set up potential
  CALL alloc_hamiltonian()

  CALL init_nabla2_sparse()
  CALL init_ilu0_prec()

  ! Local pseudopotential
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
  CALL init_strfact_shifted()
  CALL init_V_ps_loc_gaussian_G( Nparams, A, alpha )

  Nstates = 4
  ALLOCATE( Focc(Nstates) )
  Focc(:) = 1.d0

  NbetaNL = 0

  ALLOCATE( evecs(Npoints,Nstates), evals(Nstates) )
  DO ist = 1, Nstates
    DO ip = 1, Npoints
      CALL random_number( evecs(ip,ist) )
    ENDDO
  ENDDO

  ! Setup proper orthonormalization
  !
  SELECT CASE(trim(diag_method))
  CASE( 'davidson', 'Emin1', 'Emin2' )
    CALL orthonormalize( Nstates, evecs )
  CASE( 'lobpcg', 'davidson_qe2', 'davidson_qe3', 'davidson_qe4' )
    CALL ortho_gram_schmidt( evecs, Npoints, Npoints, Nstates )
  CASE DEFAULT
    ! Should not reach here
    WRITE(*,*) 'ERROR 82 Unknown diag_method = ', trim(diag_method)
    STOP 
  END SELECT 

  !
  ! Call the diagonalization subroutine
  !
  IF( trim(diag_method) == 'davidson_qe2' ) THEN 
    !
    ALLOCATE( btype(Nstates) )
    btype(:) = 1
    CALL diag_davidson_qe( Npoints, Nstates, 2*Nstates, evecs, 1.0d-4, evals, &
                           btype, notcnv, dav_iter, .TRUE. )
    !
  ELSEIF( trim(diag_method) == 'davidson_qe3' ) THEN 
    !
    ALLOCATE( btype(Nstates) )
    btype(:) = 1
    CALL diag_davidson_qe( Npoints, Nstates, 3*Nstates, evecs, 1.0d-4, evals, &
                           btype, notcnv, dav_iter, .TRUE. )
    !
  ELSEIF( trim(diag_method) == 'davidson_qe4' ) THEN 
    !
    ALLOCATE( btype(Nstates) )
    btype(:) = 1
    CALL diag_davidson_qe( Npoints, Nstates, 4*Nstates, evecs, 1.0d-4, evals, &
                           btype, notcnv, dav_iter, .TRUE. )
    !
  ELSEIF( trim(diag_method) == 'davidson' ) THEN 
    !
    CALL diag_davidson( evals, evecs, 1.0d-4, .TRUE. )
    !
  ELSEIF( trim(diag_method) == 'lobpcg' ) THEN 
    !
    CALL diag_lobpcg( evals, evecs, 1.0d-4, .TRUE. )
    !
  ELSEIF( trim(diag_method) == 'Emin1' ) THEN 
    !
    CALL Sch_solve_Emin_pcg( 1, 3.d-5, .FALSE., 1d-4, 100, .TRUE. )
    !
  ELSEIF( trim(diag_method) == 'Emin2' ) THEN 
    !
    CALL Sch_solve_Emin_pcg( 2, 3.d-5, .FALSE., 1d-4, 100, .TRUE. )
    !
  ENDIF 

  WRITE(*,*)
  WRITE(*,*) 'Final eigenvalues:'
  DO ist = 1, Nstates
    WRITE(*,'(1x,I4,F18.10)') ist, evals(ist)
  ENDDO

  ! Remember to renormalize eigenvectors for diag_davidson_qe and diag_lobpcg
  IF( trim(diag_method) == 'davidson_qe2' .OR. &
      trim(diag_method) == 'davidson_qe3' .OR. &
      trim(diag_method) == 'davidson_qe4' .OR. &
      trim(diag_method) == 'lobpcg' ) THEN 
    evecs(:,:) = evecs(:,:)/sqrt(dVol)
  ENDIF 

  WRITE(*,*)
  WRITE(*,*) 'Check normalization:'
  DO ist = 1, Nstates
    WRITE(*,'(1x,I4,F18.10)') ist, ddot( Npoints, evecs(:,ist), 1, evecs(:,ist), 1 )*dVol
  ENDDO 

  IF( allocated(btype) ) DEALLOCATE( btype )
  DEALLOCATE( evecs, evals )
  DEALLOCATE( Focc )
  CALL dealloc_atoms()
  CALL dealloc_nabla2_sparse()
  CALL dealloc_ilu0_prec()
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()


CONTAINS 

!----------------------
SUBROUTINE setup_args()
!----------------------
  INTEGER :: iargc
  CHARACTER(64) :: arg_tmp

  Nargs = iargc()
  IF( Nargs /= 4 ) THEN 
    WRITE(*,*) 'ERROR: exactly the following arguments must be given:'
    WRITE(*,*) '       N A_in alpha_in diag_method'
    WRITE(*,*)
    WRITE(*,*) 'A and alpha is Gaussian parameter:'
    WRITE(*,*) '   f(r) = A*exp( - alpha*r^2 )'
    WRITE(*,*)
    WRITE(*,*) 'with'
    WRITE(*,*) '   A     = A_in/(2.d0*PI*alpha_in**2)**1.5d0'
    WRITE(*,*) '   alpha = 0.5d0/alpha_in**2'
    WRITE(*,*)
    WRITE(*,*) 'Example call:'
    WRITE(*,*) '   ./progname.x 35 1.0 0.15 davidson_qe'
    WRITE(*,*)
    STOP 
  ENDIF 

  CALL getarg( 1, arg_tmp )
  READ(arg_tmp, *) N_in
  
  CALL getarg( 2, arg_tmp )
  READ(arg_tmp, *) A_in
  
  CALL getarg( 3, arg_tmp )
  READ(arg_tmp, *) alpha_in

  CALL getarg( 4, diag_method )
  IF( trim(diag_method) == 'davidson_qe2' ) THEN 
    WRITE(*,*)
    WRITE(*,*) 'Using davidson_qe2'
    !
  ELSEIF( trim(diag_method) == 'davidson_qe3' ) THEN 
    WRITE(*,*)
    WRITE(*,*) 'Using davidson_qe3'
    !
  ELSEIF( trim(diag_method) == 'davidson_qe4' ) THEN 
    WRITE(*,*)
    WRITE(*,*) 'Using davidson_qe4'
    !
  ELSEIF( trim(diag_method) == 'davidson' ) THEN 
    WRITE(*,*)
    WRITE(*,*) 'Using davidson'
    !
  ELSEIF( trim(diag_method) == 'lobpcg' ) THEN 
    WRITE(*,*)
    WRITE(*,*) 'Using LOBPCG'
    !
  ELSEIF( trim(diag_method) == 'Emin1' ) THEN 
    WRITE(*,*)
    WRITE(*,*) 'Using Emin1'
    !
  ELSEIF( trim(diag_method) == 'Emin2' ) THEN 
    WRITE(*,*)
    WRITE(*,*) 'Using Emin2'
    !
  ELSE 
    WRITE(*,*)
    WRITE(*,*) 'Unknown diag_method = ', trim(diag_method)
    STOP 
  ENDIF 


END SUBROUTINE 



END PROGRAM


