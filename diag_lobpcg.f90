! eFeFeR

!------------------------------------------------
SUBROUTINE diag_lobpcg( Nstates, LAMBDA, X )
!------------------------------------------------
  use m_LF3d, only : Npoints => LF3d_Npoints
  implicit none
  ! arguments
  integer :: Nstates
  integer :: ik
  real(8) :: lambda(Nstates)
  real(8) :: X(Npoints,Nstates)
  ! Local variables
  integer, parameter :: maxIter=100
  real(8), parameter :: tolerance=5.d-5
  real(8), parameter :: tfudge=1.d10
  ! Allocatable arrays
  real(8), allocatable :: Q(:,:), HQ(:,:)
  real(8), allocatable :: temp1(:,:), T(:,:), G(:,:), tempX(:,:), U(:,:)
  real(8), allocatable :: resnrm(:)
  REAL(8), ALLOCATABLE :: evals_T(:)
  !
  integer :: iter
  integer :: Nstates2,Nstates3
  integer :: nconv, ilock,nlock
  integer :: info
  real(8) :: mem
  ! Iterator
  integer :: i,j
  ! For timing
  real(8) :: time1,time2
  ! Functions
  real(8) :: ddot
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
  allocate(Q(Npoints,Nstates3)); Q(:,:) = 0.d0
  allocate(HQ(Npoints,Nstates3)); HQ(:,:) = 0.d0
  allocate(temp1(Nstates,Nstates)); temp1(:,:) = 0.d0
  allocate(T(Nstates3,Nstates3)); T(:,:) = 0.d0
  allocate(G(Nstates3,Nstates3)); G(:,:) = 0.d0
  allocate(tempX(Npoints,Nstates)); tempX(:,:) = 0.d0
  allocate(U(Nstates3,Nstates3)); U(:,:) = 0.d0
  allocate(resnrm(Nstates)); resnrm(:) = 0.d0

  ALLOCATE( evals_T(Nstates3) ); evals_T(:) = 0.d0

  mem = (7.d0*Npoints*Nstates3 + Nstates*Nstates + 3.d0*Nstates3*Nstates3)*16.0
  mem = mem + Nstates*8.0
  write(*,*) 'Allocated dynamic memory in LOBPCG = ', mem/1024.d0/1024.d0

  ! Initial wavefunction
  Q(1:Npoints,1:Nstates) = X(:,:)

!
! Apply Hamiltonian
!
  call op_H( Nstates, Q(:,1:Nstates), HQ(:,1:Nstates) )

!-----------------------------------------
! First iteration, pulled out of the loop
!-----------------------------------------
  iter = 1
  ! XHX <-- X* HX
  call dgemm('T','N',Nstates,Nstates,Npoints,1.d0,Q,Npoints,HQ,Npoints,0.d0,temp1,Nstates)

  ! Calculate residual vectors: W <-- HX - XHX
  call dcopy(Npoints*Nstates, HQ,1, Q(1,Nstates+1), 1) ! W <-- HX
  call dgemm('N','N',Npoints,Nstates,Nstates,-1.d0,Q,Npoints,temp1,Nstates, &
            1.d0,Q(1,Nstates+1),Npoints)
  ! Diagonalize
  CALL rdiaghg( Nstates, Nstates, temp1, IMat, Nstates, lambda, temp1 )
  
!
! Check convergence
!
  nconv = 0 ! Reset nconv
  nlock = 0
  do i=1,Nstates
    ! TODO: use BLAS
    resnrm(i) = sqrt( ddot(Npoints, Q(1,Nstates+i),1, Q(1,Nstates+i),1) )
    write(*,*) i, lambda(i), resnrm(i)
    if(resnrm(i) < tolerance) nconv = nconv + 1
    if(resnrm(i) < tolerance/TFUDGE) nlock = nlock + 1
  enddo

  if(nconv >= Nstates) GOTO 10

  ! Apply preconditioner
  do i=1,Nstates
!    Q(:,Nstates+i) = prec*Q(:,Nstates+i) ! element-wise multiplication
    CALL prec_ilu0_inplace( Q(:,Nstates+i) )
  ENDDO
  
  if(nlock > 0) then
    write(*,*) 'WARNING: nlock=',nlock
  ENDIF
  
!
! Apply Hamiltonian
!
  call op_H( Nstates, Q(:,Nstates+1:Nstates2), HQ(:,Nstates+1:Nstates2) )

  ! C <-- W* W
  call dgemm('T','N',Nstates,Nstates,Npoints, 1.d0,Q(1,Nstates+1),Npoints, &
       Q(1,Nstates+1),Npoints, 0.d0,temp1,Nstates)
  temp1 = ( temp1 + transpose(temp1) )*0.5d0
  
  ! Cholesky decomposition
  call dpotrf('U',Nstates,temp1,Nstates,info)
  if(info /= 0) then
    write(*,'(1x,A,I4)') 'ERROR calculating Cholesky decomposition : info ', info
    stop
  ENDIF

  ! Solve linear equations
  ! W <-- W/C
  call dtrsm('R', 'U', 'N', 'N', Npoints,Nstates, 1.d0,temp1,Nstates, Q(1,Nstates+1),Npoints)
  ! HW <-- HW/C
  call dtrsm('R', 'U', 'N', 'N', Npoints,Nstates, 1.d0,temp1,Nstates, HQ(1,Nstates+1),Npoints)


  ! T <-- Q* HQ
  call dgemm('T','N',Nstates2,Nstates2,Npoints, 1.d0,Q,Npoints, HQ,Npoints, &
     0.d0,T,Nstates3)
  T = ( T + transpose(T) )*0.5d0
  ! G <-- Q* Q
  call dgemm('T','N',Nstates2,Nstates2,Npoints, 1.d0,Q,Npoints, Q,Npoints, &
     0.d0,G,Nstates3)
  G = (G + transpose(G) )*0.5d0
  
  !call eig_zhegv_f90(T,Nstates3, G,Nstates3, U,Nstates3, Nstates2)
  CALL rdiaghg( Nstates2, Nstates2, T, G, Nstates3, evals_T, U )


  ! X <-- Q U
  call dgemm('N','N',Npoints,Nstates,Nstates2, 1.d0,Q,Npoints, U,Nstates3, &
     0.d0,tempX,Npoints)
  call dcopy(Npoints*Nstates, tempX,1, Q,1)

  ! HX <-- HQ U
  call dgemm('N','N',Npoints,Nstates,Nstates2, 1.d0,HQ,Npoints, U,Nstates3, &
     0.d0,tempX,Npoints)
  call dcopy(Npoints*Nstates, tempX,1, HQ,1)

  ! P <-- W
  call dcopy(Npoints*Nstates, Q(1,Nstates+1),1, Q(1,Nstates2+1),1)
  ! HP <-- HW
  call dcopy(Npoints*Nstates, HQ(1,Nstates+1),1, HQ(1,Nstates2+1),1)


!-----------------------------------------
! Begin of LOBPCG main iteration
!-----------------------------------------

  do iter=2,maxIter
    ! XHX <-- X* HX
    call dgemm('T','N',Nstates,Nstates,Npoints, 1.d0,Q,Npoints, HQ,Npoints, 0.d0,temp1,Nstates)
    ! Calculate residual vectors
    call dcopy(Npoints*Nstates, HQ,1, Q(1,Nstates+1),1) ! W <-- HX
    call dgemm('N','N',Npoints,Nstates,Nstates, -1.d0,Q,Npoints, temp1,Nstates, &
      1.d0,Q(1,Nstates+1),Npoints)

    !call eig_zheevd_f90(temp1,Nstates,lambda,Nstates)
    CALL rdiaghg( Nstates, Nstates, temp1, IMat, Nstates, lambda, temp1 )
       
    ! Check convergence
    nconv = 0 ! reset nconv
    nlock = 0
    do i=1,Nstates
      resnrm(i) = sqrt( ddot(Npoints, Q(1,Nstates+i),1, Q(1,Nstates+i),1) )
      if(resnrm(i) < tolerance) nconv = nconv + 1
      if(resnrm(i) < tolerance/TFUDGE) ilock = ilock + 1
    enddo
    write(*,*) 'iter = ', iter, 'nconv = ', nconv
  
    if(nconv >= Nstates-3) GOTO 10

    ! Apply preconditioner
    do i=1,Nstates
      !W(:,i) = prec*rV(:,i)
      CALL prec_ilu0_inplace( Q(:,Nstates+i) )
    enddo
  
    if(nlock > 0) then
      write(*,*) 'Warning: nlock=',nlock
    endif
  
!
! Apply Hamiltonian
!
    call op_H( Nstates, Q(:,Nstates+1:Nstates2), HQ(:,Nstates+1:Nstates2) )

    ! C <-- W* W
    call dgemm('T','N',Nstates,Nstates,Npoints, 1.d0,Q(1,Nstates+1),Npoints, &
         Q(1,Nstates+1),Npoints, 0.d0,temp1,Nstates)
    temp1  = (temp1 + transpose(temp1))*0.5d0
    
    ! Cholesky decomposition
    call dpotrf('U',Nstates,temp1,Nstates,info)
    if(info /= 0) then
      write(*,'(1x,A,I4)') 'ERROR calculating Cholesky decomposition : info ', info
      stop
    endif

    ! Solve linear equations
    ! W <-- W/C
    call dtrsm('R', 'U', 'N', 'N', Npoints,Nstates, 1.d0,temp1,Nstates, Q(1,Nstates+1),Npoints)
    ! HW <-- HW/C
    call dtrsm('R', 'U', 'N', 'N', Npoints,Nstates, 1.d0,temp1,Nstates, HQ(1,Nstates+1),Npoints)

    ! T <-- Q* HQ
    call dgemm('T','N',Nstates3,Nstates3,Npoints, 1.d0,Q,Npoints, HQ,Npoints, &
      0.d0,T,Nstates3)
    T = ( T + transpose(T) )*0.5d0
    ! G <-- Q* Q
    call dgemm('T','N',Nstates3,Nstates3,Npoints, 1.d0,Q,Npoints, Q,Npoints, &
      0.d0,G,Nstates3)
    G = ( G + transpose(G) )*0.5d0

    !call eig_zhegv_f90(T,Nstates3, G,Nstates3, U,Nstates3, Nstates3)
    CALL rdiaghg( Nstates3, Nstates3, T,G, Nstates3, evals_T, U )

    ! X <-- Q U
    call dgemm('N','N',Npoints,Nstates,Nstates3, 1.d0,Q,Npoints, U,Nstates3, &
       0.d0,tempX,Npoints)
    call dcopy(Npoints*Nstates, tempX,1, Q,1)
    ! HX <-- HQ U
    call dgemm('N','N',Npoints,Nstates,Nstates3, 1.d0,HQ,Npoints, U,Nstates3, &
      0.d0,tempX,Npoints)
    call dcopy(Npoints*Nstates, tempX,1, HQ,1)

    ! P
    call dgemm('N','N',Npoints,Nstates,Nstates, 1.d0,Q(1,Nstates2+1),Npoints, &
      U(Nstates2+1,1),Nstates3, 0.d0,tempX,Npoints)
    call dgemm('N','N',Npoints,Nstates,Nstates, 1.d0,Q(1,Nstates+1),Npoints, &
      U(Nstates+1,1),Nstates3, 0.d0,Q(1,Nstates2+1),Npoints)
    call daxpy(Npoints*Nstates, 1.d0,tempX,1, Q(1,Nstates2+1),1)

    ! HP
    call dgemm('N','N',Npoints,Nstates,Nstates, 1.d0,HQ(1,Nstates2+1),Npoints, &
      U(Nstates2+1,1),Nstates3, 0.d0,tempX,Npoints)
    call dgemm('N','N',Npoints,Nstates,Nstates, 1.d0,HQ(1,Nstates+1),Npoints, &
      U(Nstates+1,1),Nstates3, 0.d0,HQ(1,Nstates2+1),Npoints)
    call daxpy(Npoints*Nstates, 1.d0,tempX,1, HQ(1,Nstates2+1),1)

    ! C = P* P
    call dgemm('T','N',Nstates,Nstates,Npoints, &
      1.d0,Q(1,Nstates2+1),Npoints, Q(1,Nstates2+1),Npoints, 0.d0,temp1,Nstates)
    temp1 = ( temp1 + transpose(temp1) )*0.5d0
    !
    ! Cholesky decomposition
    call dpotrf('U',Nstates,temp1,Nstates,info)
    if(info /= 0) then
      write(*,'(1x,A,I4)') 'ERROR calculating Cholesky decomposition : info = ', info
      stop
    endif
    ! P = P/C
    call dtrsm('R','U','N','N', Npoints,Nstates, 1.d0,temp1,Nstates, &
      Q(1,Nstates2+1),Npoints)
    ! HP = HP/C
    call dtrsm('R','U','N','N', Npoints,Nstates, 1.d0,temp1,Nstates, &
      HQ(1,Nstates2+1),Npoints)
  enddo ! end of iteration loop

10 continue

  ! XHX = X* HX
  !call dgemm('T','N',Nstates,Nstates,Npoints,ONE,X,Npoints,HX,Npoints,ZERO,XHX_temp,Nstates)
  ! XHX = (XHX + XHX*)/2
  !call mkl_zomatadd('Col','N','T',Nstates,Nstates,HALF,XHX_temp,Nstates,&
  !    HALF,XHX_temp,Nstates, XHX,Nstates)
  ! Calculate the eigenvalues and eigenvectors
  !call eig_zheev(XHX,lambda,Nstates)
  !W = X ! save X to W, W must not be used again ...
  !call dgemm('N','N',Npoints,Nstates,Nstates, ONE,W,Npoints, XHX,Npoints, ZERO,X,Npoints)

  write(*,*) 'Number of converged eigenvalues:', nconv
  
  do i=1,Nstates
    write(*,'(1x,I6,F18.10,ES18.10)') i, lambda(i), resnrm(i)
  enddo

  deallocate(Q)
  deallocate(HQ)
  deallocate(temp1)
  deallocate(T)
  deallocate(G)
  deallocate(tempX)
  deallocate(U)
  deallocate(resnrm)
end subroutine


