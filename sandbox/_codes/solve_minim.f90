
SUBROUTINE solve_minim()
  USE m_globals, ONLY : N, evals, evecs, Nstate, Solution_Method, eVtau, eVtau_old
  IMPLICIT NONE
  INTEGER :: ist
  REAL(8), ALLOCATABLE :: Hv(:,:), Hred(:,:)

  IF( Solution_Method == 'cg' ) THEN
    CALL minimE_cg( 3.d-5, 1000, .FALSE. )
  ELSEIF( Solution_Method == 'pcg' ) THEN  
    ! for preconditioning
    ALLOCATE( eVtau(Nstate) )
    ALLOCATE( eVtau_old(Nstate) )
    eVtau = 10d0*0.03674930946797074 ! 10 eV
    eVtau_old = 0.d0
    !eVtau(:) = 10.d0
    !
    CALL minimE_pcg( 3.d-5, 1000, .FALSE. )
    !CALL minimE_pcg_Gan( 100, .FALSE. )
  ENDIF

  ALLOCATE( Hv(N**3,Nstate) )
  ALLOCATE( Hred(Nstate,Nstate) )

  DO ist = 1, Nstate
    CALL apply_Ham( evecs(:,ist), Hv(:,ist) )
  ENDDO

  CALL dgemm('T','N', Nstate,Nstate, N**3, 1.d0, evecs, N**3, Hv, N**3, 0.d0, Hred,Nstate)

  !WRITE(*,*) Hred
  CALL eig_dsyev( Hred, evals, Nstate )
  DO ist = 1, Nstate
    WRITE(*,'(1x,I5,F18.10)') ist, evals(ist)
  ENDDO

  DEALLOCATE( Hv )

END SUBROUTINE




