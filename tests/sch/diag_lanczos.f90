SUBROUTINE diag_lanczos( Ncol, evals, evecs )
  USE m_globals, ONLY : N
  IMPLICIT NONE
  !
  INTEGER :: Ncol
  REAL(8) :: evals(Ncol)
  REAL(8) :: evecs(N**3,Ncol)
  !
  REAL(8), ALLOCATABLE :: Vm(:,:), W(:)  ! Lanczos vectors
  INTEGER :: j
  REAL(8) :: r1
  REAL(8) :: alpha(Ncol), beta(Ncol)  ! alpha and beta of the triangular matrix
  REAL(8) :: evecsT(Ncol,Ncol)  ! the eigenvector of triangular matrix
  REAL(8), ALLOCATABLE :: work(:)
  INTEGER :: info
  !
  REAL(8) :: ddot

  ALLOCATE( Vm(N**3,Ncol), W(N**3) )
  ALLOCATE( work(2*Ncol-1) )

  CALL r8_rand_vec( N**3,Vm(:,1) )
  r1 = sqrt( ddot( N**3, Vm(:,1),1, Vm(:,1),1 ) )
  Vm(:,1) = Vm(:,1)/r1

  r1 = sqrt( ddot(N**3,Vm(:,1),1,Vm(:,1),1) )

  WRITE(*,*) 'r1 = ', r1
  
  DO j = 1, Ncol-1
    CALL apply_Ham( Vm(:,j), W(:) )
    alpha(j) = ddot( N**3, W(:),1, Vm(:,j),1 )
    !WRITE(*,*) 'j, alpha = ', j, alpha(j)
    IF(j>1) THEN
      W(:) = W(:) - alpha(j)*Vm(:,j) - beta(j-1)*Vm(:,j-1)
    ELSE
      W(:) = W(:) - alpha(j)*Vm(:,j) 
    ENDIF
    beta(j) = sqrt( ddot( N**3, W(:),1, W(:),1 ) )
    Vm(:,j+1) = W(:)/beta(j)
    r1 = sqrt( ddot(N**3, Vm(:,j+1),1, Vm(:,j+1),1 ) )
    WRITE(*,*) 'j, r1 = ', r1
  ENDDO

  CALL apply_Ham( Vm(:,Ncol), W(:) )
  alpha(Ncol) = ddot( N**3, W(:),1, Vm(:,ncol),1 )
  !WRITE(*,*) 'Ncol, alpha = ', Ncol,alpha(Ncol)

  IF( Ncol > 1 ) THEN
    CALL dstev('V',Ncol, alpha, beta, evecsT,Ncol, work, info )
    WRITE(*,*) 'info = ', info
    evals(:) = alpha(:)
    WRITE(*,*) 'Pass here'
  ELSE
    evals(1) = alpha(1)
  ENDIF

  WRITE(*,*) evals(:)

  DEALLOCATE( Vm, W )
  DEALLOCATE( work )


END SUBROUTINE
