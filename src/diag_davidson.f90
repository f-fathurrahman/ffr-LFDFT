! Diagonalization of Hamiltonian via Davidson method.
! This subroutine is based on Davidson subroutine used in SOCORRO.

SUBROUTINE diag_davidson( evals, v, TOLERANCE, verbose )
  
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints , dVol => LF3d_dVol
  USE m_states, ONLY : Nstates
  IMPLICIT NONE
  ! Arguments
  REAL(8) :: v(Npoints,Nstates)
  REAL(8) :: evals(Nstates)
  REAL(8) :: TOLERANCE
  LOGICAL :: verbose
  ! Local variable
  REAL(8) :: RNORM
  INTEGER :: ist, iter, MAX_DIR
  LOGICAL :: IS_CONVERGED
  REAL(8) :: MACHINE_ZERO
  REAL(8), ALLOCATABLE :: RES_TOL(:), RES_NORM(:), evals_red(:)
  REAL(8), ALLOCATABLE :: cmat(:,:), H_MAT(:,:), O_MAT(:,:), evecs(:,:)
  REAL(8), ALLOCATABLE :: HV(:,:), R(:,:), HR(:,:), xtemp(:,:)
  REAL(8), ALLOCATABLE :: evals_old(:)
  ! BLAS function
  REAL(8) :: ddot
  !
  REAL(8) :: Ebands, Ebands_old, diff_Ebands

  ALLOCATE( RES_TOL(Nstates) )
  ALLOCATE( RES_NORM(Nstates) )
  ALLOCATE( cmat(Nstates,Nstates) )

  ALLOCATE( H_MAT(2*Nstates,2*Nstates) )
  ALLOCATE( O_MAT(2*Nstates,2*Nstates) )
  ALLOCATE( evecs(2*Nstates,2*Nstates) )
  ALLOCATE( evals_red(2*Nstates) )
  ALLOCATE( evals_old(Nstates) )

  ALLOCATE( HV(Npoints,Nstates) )
  ALLOCATE( R(Npoints,Nstates) )
  ALLOCATE( HR(Npoints,Nstates) )
  ALLOCATE( xtemp(Npoints,Nstates) )

  ! Apply Hamiltonian
  CALL calc_betaNL_psi( Nstates, V )
  CALL op_H( Nstates, V, HV ) 

  ! Calculate Rayleigh quotient
  DO ist=1,Nstates
    evals(ist) = ddot(Npoints, V(:,ist),1, HV(:,ist),1)*dVol
  ENDDO

  ! Calculate matrix of residual vector
  DO ist=1,Nstates
    R(:,ist) = evals(ist)*V(:,ist) - HV(:,ist)
    RES_TOL(ist) = SQRT( ddot(Npoints, R(:,ist),1, R(:,ist),1)*dVol )
  ENDDO

  iter = 1
  IS_CONVERGED = .FALSE.
  MAX_DIR = 100
  MACHINE_ZERO = 2.220446049250313D-16
  
  RNORM = 1.D0

  evals_old(:) = 0.d0
  Ebands = sum( evals(1:Nstates) )
  Ebands_old = Ebands

  DO WHILE ( (iter <= MAX_DIR) .AND. (.NOT.IS_CONVERGED) )
    
    RES_NORM = 1.D0

    DO ist = 1,Nstates
      IF(MACHINE_ZERO < RES_TOL(ist)) RES_NORM(ist) = 1.D0/RES_TOL(ist)
    END DO

    ! Scale the residual vectors
    DO ist=1,Nstates
      R(:,ist) = RES_NORM(ist)*R(:,ist)
    ENDDO

    ! Apply preconditioner
    DO ist=1,Nstates
      call prec_ilu0_inplace( R(:,ist) )
    ENDDO

! Construct the reduced hamiltonian. The reduced hamiltonian has dimensions
!  2nb x 2nb and is constructed by filling in four nb x nb blocks one at a time:
! __ 
!|  |
!| <v|H|v>   <v|H|r>  |
!    h_mat = |  |
!| *******   <r|H|r>  | 
!|__|

    CALL calc_betaNL_psi( Nstates, R )
    CALL op_H( Nstates, R, HR )

    IF(iter == 1) THEN
      CALL DGEMM('T','N',Nstates,Nstates,Npoints, 1.d0, V, Npoints, HV,Npoints, 0.d0,cmat,Nstates)
      H_MAT(1:Nstates,1:Nstates) = cmat*dVol
    ELSE
      H_MAT(1:Nstates,1:Nstates) = 0.d0
      DO ist = 1,Nstates
        H_MAT(ist,ist) = evals(ist)  ! times dVol ??
      ENDDO
    ENDIF

    ! <v|H|r> --> cmat
    CALL DGEMM('T','N',Nstates,Nstates,Npoints, 1.d0,V,Npoints, HR,Npoints, 0.d0,cmat,Nstates)
    H_MAT(1:Nstates,Nstates+1:2*Nstates) = cmat*dVol
    CALL DGEMM('T','N',Nstates,Nstates,Npoints, 1.d0,HR,Npoints, V,Npoints, 0.d0,cmat,Nstates)
    H_MAT(Nstates+1:2*Nstates,1:Nstates) = cmat*dVol
    ! <r|H|r> --> cmat
    CALL DGEMM('T','N',Nstates,Nstates,Npoints, 1.d0,R,Npoints, HR,Npoints, 0.d0,cmat,Nstates)
    H_MAT(Nstates+1:2*Nstates,Nstates+1:2*Nstates) = cmat*dVol


! Construct the reduced overlap matrix which has dimenstions 2nb x 2nb
!   and is constructed by filling in four nb x nb blocks one at a time:
! _   _ 
!|     |
!|  <v|v>   <v|r>  |
!    o_mat = |     |
!|  *****   <r|r>  | 
!|_   _|

    O_MAT(1:Nstates,1:Nstates) = 0.d0
    DO ist = 1,Nstates
      O_MAT(ist,ist) = 1.d0
    END DO
    ! <v|r> --> cmat
    CALL DGEMM('T','N',Nstates,Nstates,Npoints, 1.d0,V,Npoints, R,Npoints, 0.d0,cmat,Nstates)
    O_MAT(1:Nstates,Nstates+1:2*Nstates) = cmat*dVol
    CALL DGEMM('T','N',Nstates,Nstates,Npoints, 1.d0,R,Npoints, V,Npoints, 0.d0,cmat,Nstates)
    O_MAT(Nstates+1:2*Nstates,1:Nstates) = cmat*dVol
    ! <r|r> --> cmat
    CALL DGEMM('T','N',Nstates,Nstates,Npoints, 1.d0,R,Npoints, R,Npoints, 0.d0,cmat,Nstates)
    O_MAT(Nstates+1:2*Nstates,Nstates+1:2*Nstates) = cmat*dVol

    CALL rdiaghg( 2*Nstates, 2*Nstates, H_MAT, O_MAT, 2*Nstates, evals_red, evecs )

    evals = evals_red(1:Nstates)
    cmat = evecs(1:Nstates,1:Nstates)
    ! V*cmat --> V
    CALL DGEMM('N','N',Npoints,Nstates,Nstates, 1.d0,V,Npoints, cmat,Nstates, 0.d0,xtemp,Npoints)
    V = xtemp
    ! HV = HV*cmat
    CALL DGEMM('N','N',Npoints,Nstates,Nstates, 1.d0,HV,Npoints, cmat,Nstates, 0.d0,xtemp,Npoints)
    HV = xtemp
    !
    cmat = evecs(Nstates+1:2*Nstates,1:Nstates)
    ! V = V + R*cmat
    CALL DGEMM('N','N',Npoints,Nstates,Nstates, 1.d0,R,Npoints, cmat,Nstates, 1.d0,V,Npoints)
    ! HV = HV + HR*cmat
    CALL DGEMM('N','N',Npoints,Nstates,Nstates, 1.d0,HR,Npoints, cmat,Nstates, 1.d0,HV,Npoints)

    ! Calculate matrix of residual vector
    DO ist=1,Nstates
      R(:,ist) = evals(ist)*V(:,ist) - HV(:,ist)
      RES_TOL(ist) = SQRT( ddot(Npoints, R(:,ist),1, R(:,ist),1)*dVol )
    ENDDO

    RNORM = SUM( abs(evals - evals_old) )/REAL(Nstates, kind=8)
    !
    Ebands = sum( evals(1:Nstates) )
    diff_Ebands = abs( Ebands - Ebands_old )
    !
    IF( verbose ) THEN 
      WRITE(*,*)
      WRITE(*,'(1x,A,I8,ES18.10)') 'Davidson: iter, RNORM ', iter, RNORM
      WRITE(*,'(1x,A,I8,F18.10,ES18.10)') 'Davidson Ebands ', iter, Ebands, diff_Ebands
      WRITE(*,*) 'Eigenvalues convergence:'
      DO ist = 1,Nstates
        WRITE(*,'(1X,I5,F18.10,ES18.10)') ist, evals(ist), abs( evals(ist)-evals_old(ist) )
      ENDDO 
    ENDIF 

    IS_CONVERGED = RNORM <= TOLERANCE

    evals_old(:) = evals(:)
    Ebands_old = Ebands

    iter = iter + 1
    !
    IF( verbose ) THEN 
      flush(6)
    ENDIF 
  END DO

  DEALLOCATE(RES_TOL)
  DEALLOCATE(RES_NORM)
  DEALLOCATE(cmat)
  DEALLOCATE(H_MAT)
  DEALLOCATE(O_MAT)
  DEALLOCATE(evecs)
  DEALLOCATE(evals_red)
  DEALLOCATE(HV)
  DEALLOCATE(R)
  DEALLOCATE(HR)
  DEALLOCATE(xtemp)
END SUBROUTINE


