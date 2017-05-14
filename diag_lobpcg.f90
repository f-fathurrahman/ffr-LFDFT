! eFeFeR

!------------------------------------------------
SUBROUTINE diag_lobpcg( LAMBDA, X, tolerance )
!------------------------------------------------
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, dVol => LF3d_dVol
  USE m_states, ONLY : Nstates
  IMPLICIT NONE 
  ! arguments
  REAL(8) :: lambda(Nstates)
  REAL(8) :: X(Npoints,Nstates)
  REAL(8) :: tolerance
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
  INTEGER :: iter
  INTEGER :: Nstates2,Nstates3
  INTEGER :: nconv, ilock,nlock
  INTEGER :: info
  real(8) :: mem
  ! Iterator
  INTEGER :: i
  ! Functions
  REAL(8) :: ddot
  !
  REAL(8), ALLOCATABLE :: IMat(:,:)

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
  WRITE(*,*) 'Allocated dynamic memory in LOBPCG = ', mem/1024.d0/1024.d0

  ! Initial wavefunction
  Q(1:Npoints,1:Nstates) = X(:,:)*sqrt(dVol)

!
! Apply Hamiltonian
!
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
  WRITE(*,*)
  WRITE(*,*) 'Eigenvalue convergence:'
  DO i=1,Nstates
    ! TODO: use BLAS
    !resnrm(i) = sqrt( ddot(Npoints, Q(1,Nstates+i),1, Q(1,Nstates+i),1) )
    resnrm(i) = abs(lambda(i)-lambda_old(i))
    WRITE(*,*) i, lambda(i), resnrm(i)
    IF(resnrm(i) < tolerance) nconv = nconv + 1
    IF(resnrm(i) < tolerance/TFUDGE) nlock = nlock + 1
  ENDDO

  IF(nlock > 0) THEN
    WRITE(*,*) 'WARNING: nlock=',nlock
  ENDIF

  IF(nconv >= Nstates) GOTO 10

  ! Apply preconditioner
  DO i=1,Nstates
    CALL prec_ilu0_inplace( Q(:,Nstates+i) )
  ENDDO
  
!
! Apply Hamiltonian
!
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
    WRITE(*,*)
    WRITE(*,*) 'Eigenvalues convergence:'
    DO i=1,Nstates
      !resnrm(i) = sqrt( ddot(Npoints, Q(1,Nstates+i),1, Q(1,Nstates+i),1) )
      resnrm(i) = abs(lambda(i)-lambda_old(i))
      WRITE(*,*) i, lambda(i), resnrm(i)
      IF(resnrm(i) < tolerance) nconv = nconv + 1
      IF(resnrm(i) < tolerance/TFUDGE) ilock = ilock + 1
    ENDDO
    WRITE(*,*) 'LOBPCG iter = ', iter, 'nconv = ', nconv
  
    IF(nconv >= Nstates) GOTO 10

    ! Apply preconditioner
    DO i=1,Nstates
      CALL prec_ilu0_inplace( Q(:,Nstates+i) )
    ENDDO 
  
    IF(nlock > 0) THEN
      WRITE(*,*) 'Warning: nlock=',nlock
    ENDIF 
  
!
! Apply Hamiltonian
!
    CALL op_H( Nstates, Q(:,Nstates+1:Nstates2), HQ(:,Nstates+1:Nstates2) )

    ! C <-- W* W
    CALL dgemm('T','N',Nstates,Nstates,Npoints, 1.d0,Q(1,Nstates+1),Npoints, &
         Q(1,Nstates+1),Npoints, 0.d0,temp1,Nstates)
    temp1  = (temp1 + transpose(temp1))*0.5d0
    
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
      WRITE(*,'(1x,A,I4)') 'ERROR calculating Cholesky decomposition : info = ', info
      STOP
    ENDIF 
    ! P = P/C
    CALL dtrsm('R','U','N','N', Npoints,Nstates, 1.d0,temp1,Nstates, &
      Q(1,Nstates2+1),Npoints)
    ! HP = HP/C
    CALL dtrsm('R','U','N','N', Npoints,Nstates, 1.d0,temp1,Nstates, &
      HQ(1,Nstates2+1),Npoints)
  ENDDO ! end of iteration loop

10 CONTINUE

  ! XHX = X* HX
  !call dgemm('T','N',Nstates,Nstates,Npoints,ONE,X,Npoints,HX,Npoints,ZERO,XHX_temp,Nstates)
  ! XHX = (XHX + XHX*)/2
  !call mkl_zomatadd('Col','N','T',Nstates,Nstates,HALF,XHX_temp,Nstates,&
  !    HALF,XHX_temp,Nstates, XHX,Nstates)
  ! Calculate the eigenvalues and eigenvectors
  !call eig_zheev(XHX,lambda,Nstates)
  !W = X ! save X to W, W must not be used again ...
  !call dgemm('N','N',Npoints,Nstates,Nstates, ONE,W,Npoints, XHX,Npoints, ZERO,X,Npoints)

!  WRITE(*,*) 'Number of converged eigenvalues:', nconv
  
!  DO i=1,Nstates
!    WRITE(*,'(1x,I6,F18.10,ES18.10)') i, lambda(i), resnrm(i)
!  ENDDO

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



