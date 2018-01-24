! eFeFeR

!------------------------------------------------
SUBROUTINE diag_lobpcg( LAMBDA, X, tolerance, verbose )
!------------------------------------------------
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_states, ONLY : Nstates
  IMPLICIT NONE 
  ! arguments
  REAL(8), INTENT(inout) :: lambda(Nstates)
  REAL(8), INTENT(inout) :: X(Npoints,Nstates)
  REAL(8), INTENT(in) :: tolerance
  LOGICAL :: verbose
  ! Local variables
  INTEGER, PARAMETER :: maxIter=100
  REAL(8), PARAMETER :: tfudge=1.d10
  ! Allocatable arrays
  REAL(8), ALLOCATABLE :: Q(:,:), HQ(:,:)
  REAL(8), ALLOCATABLE :: temp1(:,:), T(:,:), G(:,:), tempX(:,:), U(:,:)
  REAL(8), ALLOCATABLE :: resnrm(:)
  REAL(8), ALLOCATABLE :: evals_T(:)
  REAL(8), ALLOCATABLE :: lambda_old(:)
  !
  LOGICAL :: IS_CONVERGED
  INTEGER :: iter
  INTEGER :: Nstates2,Nstates3
  INTEGER :: nconv, ilock,nlock
  INTEGER :: info
  REAL(8) :: mem
  REAL(8) :: Ebands, Ebands_old, diff_Ebands
  REAL(8) :: RNORM
  ! Iterator
  INTEGER :: i
  !
  REAL(8), ALLOCATABLE :: IMat(:,:)


  IS_CONVERGED = .FALSE.

  ALLOCATE( IMat(Nstates,Nstates) )
  IMat(:,:) = 0.d0
  DO i = 1,Nstates
    IMat(i,i) = 1.d0
  ENDDO 

  Nstates2 = Nstates*2
  Nstates3 = Nstates*3

  ! Allocate memory
  ALLOCATE(Q(Npoints,Nstates3)); Q(:,:) = 0.d0
  ALLOCATE(HQ(Npoints,Nstates3)); HQ(:,:) = 0.d0
  ALLOCATE(temp1(Nstates,Nstates)); temp1(:,:) = 0.d0
  ALLOCATE(T(Nstates3,Nstates3)); T(:,:) = 0.d0
  ALLOCATE(G(Nstates3,Nstates3)); G(:,:) = 0.d0
  ALLOCATE(tempX(Npoints,Nstates)); tempX(:,:) = 0.d0
  ALLOCATE(U(Nstates3,Nstates3)); U(:,:) = 0.d0
  ALLOCATE(resnrm(Nstates)); resnrm(:) = 0.d0
  ALLOCATE(lambda_old(Nstates)); lambda_old(:) = 0.d0

  ALLOCATE( evals_T(Nstates3) ); evals_T(:) = 0.d0

  mem = (7.d0*Npoints*Nstates3 + Nstates*Nstates + 3.d0*Nstates3*Nstates3)*16.0
  mem = mem + Nstates*8.0
  IF( verbose ) THEN 
    WRITE(*,'(1x,A,F10.3)') 'Allocated dynamic memory in LOBPCG = ', mem/1024.d0/1024.d0
  ENDIF 

  ! Initial wavefunction
  Q(1:Npoints,1:Nstates) = X(:,:)

!
! Apply Hamiltonian
!
  CALL calc_betaNL_psi( Nstates, Q(:,1:Nstates) )
  CALL op_H( Nstates, Q(:,1:Nstates), HQ(:,1:Nstates) )

!-----------------------------------------
! First iteration, pulled out of the loop
!-----------------------------------------
  iter = 1
  ! XHX <-- X* HX
  CALL dgemm('T','N',Nstates,Nstates,Npoints,1.d0,Q,Npoints,HQ,Npoints,0.d0,temp1,Nstates)

  ! Calculate residual vectors: W <-- HX - XHX
  CALL dcopy(Npoints*Nstates, HQ,1, Q(1,Nstates+1), 1) ! W <-- HX
  CALL dgemm('N','N',Npoints,Nstates,Nstates,-1.d0,Q,Npoints,temp1,Nstates, &
            1.d0,Q(1,Nstates+1),Npoints)
  ! Diagonalize
  lambda_old(:) = lambda(:)
  CALL rdiaghg( Nstates, Nstates, temp1, IMat, Nstates, lambda, temp1 )
  
!
! Check convergence
!
  nconv = 0 ! Reset nconv
  nlock = 0

  DO i=1,Nstates
    resnrm(i) = abs(lambda(i)-lambda_old(i))
    IF(resnrm(i) < tolerance) nconv = nconv + 1
    IF(resnrm(i) < tolerance/TFUDGE) nlock = nlock + 1
  ENDDO

  RNORM = SUM( abs(lambda - lambda_old) )/REAL(Nstates, kind=8)
  !
  Ebands = sum( lambda(1:Nstates) )
  diff_Ebands = 1.d0
  Ebands_old = Ebands
  !
  IF( verbose ) THEN 
    WRITE(*,*)
    WRITE(*,'(1x,A,I8,ES18.10)') 'LOBPCG: iter, rnorm', 1, RNORM
    WRITE(*,'(1x,A,I8,F18.10,ES18.10)') 'LOBPCG Ebands ', iter, Ebands, diff_Ebands
    WRITE(*,*) 'Eigenvalues convergence:'
    DO i = 1,Nstates
      WRITE(*,'(I8,F18.10,ES18.10)') i, lambda(i), resnrm(i)
    ENDDO 
  ENDIF 

  IF(nlock > 0) THEN
    WRITE(*,*) 'WARNING: nlock=',nlock
  ENDIF

  IS_CONVERGED = RNORM <= TOLERANCE
  !
  IF(nconv >= Nstates) THEN
    WRITE(*,*)
    WRITE(*,*) 'LOBPCG: Convergence achieved based on nconv'
    GOTO 10
  ELSEIF( IS_CONVERGED ) THEN 
    WRITE(*,*)
    WRITE(*,*) 'LOBPCG: Convergence achieved based on RNORM'
    GOTO 10
  ENDIF 


  ! Apply preconditioner
  DO i=1,Nstates
    CALL prec_ilu0_inplace( Q(:,Nstates+i) )
  ENDDO
  
!
! Apply Hamiltonian
!
  CALL calc_betaNL_psi( Nstates, Q(:,Nstates+1:Nstates2) )
  CALL op_H( Nstates, Q(:,Nstates+1:Nstates2), HQ(:,Nstates+1:Nstates2) )

  ! C <-- W* W
  CALL dgemm('T','N',Nstates,Nstates,Npoints, 1.d0,Q(1,Nstates+1),Npoints, &
       Q(1,Nstates+1),Npoints, 0.d0,temp1,Nstates)
  temp1 = ( temp1 + transpose(temp1) )*0.5d0
  
  ! Cholesky decomposition
  CALL dpotrf('U',Nstates,temp1,Nstates,info)
  IF(info /= 0) THEN 
    WRITE(*,'(1x,A,I4)') 'ERROR calculating Cholesky decomposition : info ', info
    STOP 
  ENDIF

  ! Solve linear equations
  ! W <-- W/C
  CALL dtrsm('R', 'U', 'N', 'N', Npoints,Nstates, 1.d0,temp1,Nstates, Q(1,Nstates+1),Npoints)
  ! HW <-- HW/C
  CALL dtrsm('R', 'U', 'N', 'N', Npoints,Nstates, 1.d0,temp1,Nstates, HQ(1,Nstates+1),Npoints)


  ! T <-- Q* HQ
  CALL dgemm('T','N',Nstates2,Nstates2,Npoints, 1.d0,Q,Npoints, HQ,Npoints, &
     0.d0,T,Nstates3)
  T = ( T + transpose(T) )*0.5d0
  ! G <-- Q* Q
  CALL dgemm('T','N',Nstates2,Nstates2,Npoints, 1.d0,Q,Npoints, Q,Npoints, &
     0.d0,G,Nstates3)
  G = (G + transpose(G) )*0.5d0
  
  !call eig_zhegv_f90(T,Nstates3, G,Nstates3, U,Nstates3, Nstates2)
  CALL rdiaghg( Nstates2, Nstates2, T, G, Nstates3, evals_T, U )


  ! X <-- Q U
  CALL dgemm('N','N',Npoints,Nstates,Nstates2, 1.d0,Q,Npoints, U,Nstates3, &
     0.d0,tempX,Npoints)
  CALL dcopy(Npoints*Nstates, tempX,1, Q,1)

  ! HX <-- HQ U
  CALL dgemm('N','N',Npoints,Nstates,Nstates2, 1.d0,HQ,Npoints, U,Nstates3, &
     0.d0,tempX,Npoints)
  CALL dcopy(Npoints*Nstates, tempX,1, HQ,1)

  ! P <-- W
  CALL dcopy(Npoints*Nstates, Q(1,Nstates+1),1, Q(1,Nstates2+1),1)
  ! HP <-- HW
  CALL dcopy(Npoints*Nstates, HQ(1,Nstates+1),1, HQ(1,Nstates2+1),1)


!-----------------------------------------
! Begin of LOBPCG main iteration
!-----------------------------------------

  DO iter=2,maxIter
    ! XHX <-- X* HX
    CALL dgemm('T','N',Nstates,Nstates,Npoints, 1.d0,Q,Npoints, HQ,Npoints, 0.d0,temp1,Nstates)
    ! Calculate residual vectors
    CALL dcopy(Npoints*Nstates, HQ,1, Q(1,Nstates+1),1) ! W <-- HX
    CALL dgemm('N','N',Npoints,Nstates,Nstates, -1.d0,Q,Npoints, temp1,Nstates, &
      1.d0,Q(1,Nstates+1),Npoints)

    !call eig_zheevd_f90(temp1,Nstates,lambda,Nstates)
    lambda_old(:) = lambda(:)
    CALL rdiaghg( Nstates, Nstates, temp1, IMat, Nstates, lambda, temp1 )
       
    ! Check convergence
    nconv = 0 ! reset nconv
    nlock = 0
    DO i = 1,Nstates
      ! NOTE: This is different from the original LOBPCG in KSSOLV/Knyazev
      resnrm(i) = abs(lambda(i)-lambda_old(i))
      IF(resnrm(i) < tolerance) nconv = nconv + 1
      IF(resnrm(i) < tolerance/TFUDGE) ilock = ilock + 1
    ENDDO
    !
    RNORM = SUM( abs(lambda - lambda_old) )/REAL(Nstates, kind=8)
    Ebands = sum( lambda(1:Nstates) )
    diff_Ebands = abs( Ebands - Ebands_old )
    Ebands_old = Ebands
    !
    IF( verbose ) THEN 
      WRITE(*,*)
      WRITE(*,'(1x,A,I8,ES18.10)') 'LOBPCG: iter, RNORM ', iter, RNORM
      WRITE(*,'(1x,A,I8,F18.10,ES18.10)') 'LOBPCG Ebands ', iter, Ebands, diff_Ebands
      WRITE(*,*) 'Eigenvalues convergence:'
      DO i = 1,Nstates
        WRITE(*,'(I8,F18.10,ES18.10)') i, lambda(i), resnrm(i)
      ENDDO 
    ENDIF 

    IS_CONVERGED = RNORM <= TOLERANCE
    !
    IF(nconv >= Nstates) THEN
      WRITE(*,*)
      WRITE(*,*) 'LOBPCG: Convergence achieved based on nconv'
      GOTO 10
    ELSEIF( IS_CONVERGED ) THEN 
      WRITE(*,*)
      WRITE(*,*) 'LOBPCG: Convergence achieved based on RNORM'
      GOTO 10
    ENDIF 

    ! Apply preconditioner
    DO i = 1,Nstates
      CALL prec_ilu0_inplace( Q(:,Nstates+i) )
    ENDDO 
  
    IF( nlock > 0 ) THEN
      WRITE(*,*)
      WRITE(*,'(1x,A,I8)') 'LOBPCG: Warning: nlock = ',nlock
    ENDIF 
  
!
! Apply Hamiltonian
!
    CALL calc_betaNL_psi( Nstates, Q(:,Nstates+1:Nstates2) )
    CALL op_H( Nstates, Q(:,Nstates+1:Nstates2), HQ(:,Nstates+1:Nstates2) )

    ! C <-- W* W
    CALL dgemm('T','N',Nstates,Nstates,Npoints, 1.d0,Q(1,Nstates+1),Npoints, &
         Q(1,Nstates+1),Npoints, 0.d0,temp1,Nstates)
    temp1  = (temp1 + transpose(temp1))*0.5d0
    
    ! Cholesky decomposition
    CALL dpotrf('U',Nstates,temp1,Nstates,info)
    IF(info /= 0) THEN 
      WRITE(*,*)
      WRITE(*,'(1x,A,I4)') 'LOBPCG: ERROR calculating Cholesky decomposition : info ', info
      STOP 
    ENDIF

    ! Solve linear equations
    ! W <-- W/C
    CALL dtrsm('R', 'U', 'N', 'N', Npoints,Nstates, 1.d0,temp1,Nstates, Q(1,Nstates+1),Npoints)
    ! HW <-- HW/C
    CALL dtrsm('R', 'U', 'N', 'N', Npoints,Nstates, 1.d0,temp1,Nstates, HQ(1,Nstates+1),Npoints)

    ! T <-- Q* HQ
    CALL dgemm('T','N',Nstates3,Nstates3,Npoints, 1.d0,Q,Npoints, HQ,Npoints, &
      0.d0,T,Nstates3)
    T = ( T + transpose(T) )*0.5d0
    ! G <-- Q* Q
    CALL dgemm('T','N',Nstates3,Nstates3,Npoints, 1.d0,Q,Npoints, Q,Npoints, &
      0.d0,G,Nstates3)
    G = ( G + transpose(G) )*0.5d0

    CALL rdiaghg( Nstates3, Nstates3, T,G, Nstates3, evals_T, U )

    ! X <-- Q U
    CALL dgemm('N','N',Npoints,Nstates,Nstates3, 1.d0,Q,Npoints, U,Nstates3, &
       0.d0,tempX,Npoints)
    CALL dcopy(Npoints*Nstates, tempX,1, Q,1)
    ! HX <-- HQ U
    CALL dgemm('N','N',Npoints,Nstates,Nstates3, 1.d0,HQ,Npoints, U,Nstates3, &
      0.d0,tempX,Npoints)
    CALL dcopy(Npoints*Nstates, tempX,1, HQ,1)

    ! P
    CALL dgemm('N','N',Npoints,Nstates,Nstates, 1.d0,Q(1,Nstates2+1),Npoints, &
      U(Nstates2+1,1),Nstates3, 0.d0,tempX,Npoints)
    CALL dgemm('N','N',Npoints,Nstates,Nstates, 1.d0,Q(1,Nstates+1),Npoints, &
      U(Nstates+1,1),Nstates3, 0.d0,Q(1,Nstates2+1),Npoints)
    CALL daxpy(Npoints*Nstates, 1.d0,tempX,1, Q(1,Nstates2+1),1)

    ! HP
    CALL dgemm('N','N',Npoints,Nstates,Nstates, 1.d0,HQ(1,Nstates2+1),Npoints, &
      U(Nstates2+1,1),Nstates3, 0.d0,tempX,Npoints)
    CALL dgemm('N','N',Npoints,Nstates,Nstates, 1.d0,HQ(1,Nstates+1),Npoints, &
      U(Nstates+1,1),Nstates3, 0.d0,HQ(1,Nstates2+1),Npoints)
    CALL daxpy(Npoints*Nstates, 1.d0,tempX,1, HQ(1,Nstates2+1),1)

    ! C = P* P
    CALL dgemm('T','N',Nstates,Nstates,Npoints, &
      1.d0,Q(1,Nstates2+1),Npoints, Q(1,Nstates2+1),Npoints, 0.d0,temp1,Nstates)
    temp1 = ( temp1 + transpose(temp1) )*0.5d0
    !
    ! Cholesky decomposition
    CALL dpotrf('U',Nstates,temp1,Nstates,info)
    IF(info /= 0) THEN 
      WRITE(*,*)
      WRITE(*,'(1x,A,I4)') 'LOBPCG: ERROR calculating Cholesky decomposition : info = ', info
      STOP
    ENDIF 
    ! P = P/C
    CALL dtrsm('R','U','N','N', Npoints,Nstates, 1.d0,temp1,Nstates, &
      Q(1,Nstates2+1),Npoints)
    ! HP = HP/C
    CALL dtrsm('R','U','N','N', Npoints,Nstates, 1.d0,temp1,Nstates, &
      HQ(1,Nstates2+1),Npoints)

    IF( verbose ) THEN 
      flush(6)
    ENDIF 

  ENDDO ! end of iteration loop

10 CONTINUE

  X = Q(1:Npoints,1:Nstates)

  DEALLOCATE(lambda_old)
  DEALLOCATE(Q)
  DEALLOCATE(HQ)
  DEALLOCATE(temp1)
  DEALLOCATE(T)
  DEALLOCATE(G)
  DEALLOCATE(tempX)
  DEALLOCATE(U)
  DEALLOCATE(resnrm)
END SUBROUTINE 



