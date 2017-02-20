! eFeFeR

!------------------------------------------------
SUBROUTINE LOBPCG_CPU(ik,LAMBDA,X,NBASIS,NSTATES,NGWX,rhoe,nnr)
!------------------------------------------------
  implicit none
  ! arguments
  integer :: nbasis
  integer :: nstates
  integer :: ik
  real(8) :: lambda(nstates)
  complex(8) :: X(nbasis,nstates)
  logical :: t_prec
  integer :: NGWX
  integer :: NNR
  real(8) :: rhoe(nnr)
  ! Local variables
  integer, parameter :: maxIter=100
  real(8), parameter :: tolerance=5.d-5
  real(8), parameter :: tfudge=1.d10
  complex(8), parameter :: Z_ZERO=cmplx(0.d0,0.d0)
  complex(8), parameter :: Z_ONE=cmplx(1.d0,0.d0)
  complex(8), parameter :: MZ_ONE=cmplx(-1.d0,0.d0)
  complex(8), parameter :: Z_HALF=cmplx(0.5d0,0.d0)
  ! Allocatable arrays
  complex(8), allocatable :: Q(:,:), HQ(:,:)
  complex(8), allocatable :: temp1(:,:), T(:,:), G(:,:), tempX(:,:), U(:,:)
  real(8), allocatable :: resnrm(:)
  !
  real(8) :: prec(nbasis) ! TODO: data types? complex or real?
  integer :: iter
  integer :: nstates2,nstates3
  integer :: nconv, ilock,nlock
  integer :: info
  real(8) :: mem
  complex(8) :: tz
  ! Iterator
  integer :: i,j
  ! For timing
  real(8) :: time1,time2
  ! Functions
  complex(8) :: zdotc

  nstates2 = nstates*2
  nstates3 = nstates*3

! Allocate memory
  allocate(Q(nbasis,nstates3)); Q(:,:) = Z_ZERO
  allocate(HQ(nbasis,nstates3)); HQ(:,:) = Z_ZERO
  allocate(temp1(nstates,nstates)); temp1(:,:) = Z_ZERO
  allocate(T(nstates3,nstates3)); T(:,:) = Z_ZERO
  allocate(G(nstates3,nstates3)); G(:,:) = Z_ZERO
  allocate(tempX(nbasis,nstates)); tempX(:,:) = Z_ZERO
  allocate(U(nstates3,nstates3)); U(:,:) = Z_ZERO
  allocate(resnrm(nstates)); resnrm(:) = 0.d0

  mem = (7.d0*nbasis*nstates3 + nstates*nstates + 3.d0*nstates3*nstates3)*16.0
  mem = mem + nstates*8.0
  write(*,*) 'Allocated dynamic memory in LOBPCG = ', mem/1024.d0/1024.d0

  ! Initial wavefunction
  Q(1:nbasis,1:nstates) = cmplx(0.d0,0.d0)
  do i=1,nstates
    Q(i,i) = cmplx(10.d0,0.d0)
  enddo

  call mytimer(time1)

  ! Orthogonalize the initial wavefunctions
  call ortho_qr_f90(Q(1:nbasis,1:nstates),nbasis,nstates)
  write(*,*) 'sum(Q(1:nbasis,1:nstates)) = ', sum(Q(1:nbasis,1:nstates))

!
! Apply Hamiltonian
!
  do i=1,nstates
   call apply_Ham(Q(:,i), HQ(:,i), nbasis, ik, rhoe, nnr) 
  enddo
  write(*,*) 'After apply_Ham:'
  write(*,*) 'sum(X) = ', sum(Q(1:nbasis,1:nstates))
  write(*,*) 'sum(HX) = ', sum(HQ(1:nbasis,1:nstates))

  !stop 'breakpoint 80'

  if(t_prec) then
    !write(*,*) 'Generating preconditioner ...'
    !call genprec(prec,nbasis)
    write(*,*) 'No preconditioning ...'
  else
    prec = 1.d0 ! Not a useful preconditioner ;-)
  endif

!-----------------------------------------
! First iteration, pulled out of the loop
!-----------------------------------------
  iter = 1
  ! XHX <-- X* HX
  call zgemm('C','N',nstates,nstates,nbasis,Z_ONE,Q,nbasis,HQ,nbasis,Z_ZERO,temp1,nstates)
  write(*,*) 'Pass 1'
  write(*,*) 'sum(X) = ', sum(Q(1:nbasis,1:nstates))
  write(*,*) 'sum(HX) = ', sum(HQ(1:nbasis,1:nstates))
  write(*,*) 'sum(temp1) = ', sum(temp1)

  ! Calculate residual vectors: W <-- HX - XHX
  call zcopy(nbasis*nstates, HQ,1, Q(1,nstates+1), 1) ! W <-- HX
  call zgemm('N','N',nbasis,nstates,nstates,MZ_ONE,Q,nbasis,temp1,nstates, &
            Z_ONE,Q(1,nstates+1),nbasis)
  write(*,*) 'sum(Q(:,nstates+1:nstates2) = ', sum(Q(:,nstates+1:nstates2))
  ! Diagonalize
  write(*,*) 'before eig_zheevd_f90: sum(temp1) = ', sum(temp1)
  call eig_zheevd_f90(temp1,nstates,lambda,nstates)
  write(*,*) 'sum(lambda) = ', sum(lambda)
  
!
! Check convergence
!
  nconv = 0 ! Reset nconv
  nlock = 0
  do i=1,nstates
    ! TODO: use BLAS
    resnrm(i) = sqrt( zdotc(nbasis, Q(1,nstates+i),1, Q(1,nstates+i),1) )
    write(*,*) i, lambda(i), resnrm(i)
    if(resnrm(i) < tolerance) nconv = nconv + 1
    if(resnrm(i) < tolerance/TFUDGE) nlock = nlock + 1
  enddo

  if(nconv >= nstates-3) GOTO 10

  ! Apply preconditioner
  !do i=1,nstates
  !  Q(:,nstates+i) = prec*Q(:,nstates+i) ! element-wise multiplication
  !enddo
  
  if(nlock > 0) then
    write(*,*) 'WARNING: nlock=',nlock
  endif
  
!
! Apply Hamiltonian
!
  write(*,*) 'Before apply_Ham:'
  write(*,*) 'sum(Q(:,nstates+1:nstates2) = ', sum(Q(:,nstates+1:nstates2))
  write(*,*) 'sum(HQ(:,nstates+1:nstates2) = ', sum(HQ(:,nstates+1:nstates2))
  do i=1,nstates
   call apply_Ham(Q(1,nstates+i), HQ(1,nstates+i), nbasis, ik, rhoe, nnr) 
  enddo
  write(*,*) 'After apply_Ham:'
  write(*,*) 'sum(Q(:,nstates+1:nstates2) = ', sum(Q(:,nstates+1:nstates2))
  write(*,*) 'sum(HQ(:,nstates+1:nstates2) = ', sum(HQ(:,nstates+1:nstates2))

  ! C <-- W* W
  call zgemm('C','N',nstates,nstates,nbasis, Z_ONE,Q(1,nstates+1),nbasis, &
       Q(1,nstates+1),nbasis, Z_ZERO,temp1,nstates)
  temp1 = ( temp1 + conjg(temp1) )*0.5d0
  
  write(*,*) 'before chol: sum(temp1) = ', sum(temp1)
  ! Cholesky decomposition
  call zpotrf('U',nstates,temp1,nstates,info)
  if(info /= 0) then
    write(*,'(1x,A,I4)') 'ERROR calculating Cholesky decomposition : info ', info
    stop
  endif
  write(*,*) 'after chol: sum(temp1) = ', sum(temp1)

  ! Solve linear equations
  ! W <-- W/C
  call ztrsm('R', 'U', 'N', 'N', nbasis,nstates, Z_ONE,temp1,nstates, Q(1,nstates+1),nbasis)
  ! HW <-- HW/C
  call ztrsm('R', 'U', 'N', 'N', nbasis,nstates, Z_ONE,temp1,nstates, HQ(1,nstates+1),nbasis)
  write(*,*) 'After ztrsm:'
  write(*,*) 'sum(Q(:,nstates+1:nstates2) = ', sum(Q(:,nstates+1:nstates2))
  write(*,*) 'sum(HQ(:,nstates+1:nstates2) = ', sum(HQ(:,nstates+1:nstates2))

  ! T <-- Q* HQ
  write(*,*) 'sum(Q(:,1:nstates2) = ', sum(Q(:,1:nstates2))
  write(*,*) 'sum(HQ(:,1:nstates2) = ', sum(HQ(:,1:nstates2))
  call zgemm('C','N',nstates2,nstates2,nbasis, Z_ONE,Q,nbasis, HQ,nbasis, &
     Z_ZERO,T,nstates3)
  T = ( T + conjg(T) )*0.5d0
  ! G <-- Q* Q
  call zgemm('C','N',nstates2,nstates2,nbasis, Z_ONE,Q,nbasis, Q,nbasis, &
     Z_ZERO,G,nstates3)
  G = (G + conjg(G) )*0.5d0
  
  write(*,*) 'sum(T) = ', sum(T)
  write(*,*) 'sum(G) = ', sum(G)
  write(*,*) 'before: sum(U) = ', sum(U) 
  !call eig_zhegv_f90(T(1:nstates2,1:nstates2),nstates2, G(1:nstates2,1:nstates2),nstates2, &
  !   U(1:nstates2,1:nstates2),nstates2, nstates2)
  call eig_zhegv_f90(T,nstates3, G,nstates3, U,nstates3, nstates2)
  !call eig_zhegv(T(1:nstates2,1:nstates2),G(1:nstates2,1:nstates2), &
  !  U(1:nstates2,1:nstates2), nstates2)
  write(*,*) 'sum(U(1:nstates2,1:nstates2)) = ', sum(U(1:nstates2,1:nstates2))
  write(*,*) 'sum(U(1:nstates2,1:nstates)) = ', sum(U(1:nstates2,1:nstates))

  ! X <-- Q U
  call zgemm('N','N',nbasis,nstates,nstates2, Z_ONE,Q,nbasis, U,nstates3, &
     Z_ZERO,tempX,nbasis)
  call zcopy(nbasis*nstates, tempX,1, Q,1)
  write(*,*) 'sum(X) = ', sum(tempX)
  write(*,*) 'sum(Q(:,1:nstates) = ', sum(Q(:,1:nstates))
  ! HX <-- HQ U
  call zgemm('N','N',nbasis,nstates,nstates2, Z_ONE,HQ,nbasis, U,nstates3, &
     Z_ZERO,tempX,nbasis)
  call zcopy(nbasis*nstates, tempX,1, HQ,1)
  write(*,*) 'sum(X) = ', sum(tempX)
  write(*,*) 'sum(HQ(:,1:nstates)) = ', sum(HQ(:,1:nstates))
  ! P <-- W
  call zcopy(nbasis*nstates, Q(1,nstates+1),1, Q(1,nstates2+1),1)
  ! HP <-- HW
  call zcopy(nbasis*nstates, HQ(1,nstates+1),1, HQ(1,nstates2+1),1)
  write(*,*) 'sum(Q(:,nstates2+1:nstates3)) = ', sum(Q(:,nstates2+1:nstates3))
  write(*,*) 'sum(HQ(:,nstates2+1:nstates3)) = ', sum(HQ(:,nstates2+1:nstates3))

!-----------------------------------------
! Begin of LOBPCG main iteration
!-----------------------------------------

  do iter=2,maxIter
    ! XHX <-- X* HX
    call zgemm('C','N',nstates,nstates,nbasis, Z_ONE,Q,nbasis, HQ,nbasis, Z_ZERO,temp1,nstates)
    ! Calculate residual vectors
    call zcopy(nbasis*nstates, HQ,1, Q(1,nstates+1),1) ! W <-- HX
    call zgemm('N','N',nbasis,nstates,nstates, MZ_ONE,Q,nbasis, temp1,nstates, &
      Z_ONE,Q(1,nstates+1),nbasis)

    call eig_zheevd_f90(temp1,nstates,lambda,nstates)
       
    ! Check convergence
    nconv = 0 ! reset nconv
    nlock = 0
    do i=1,nstates
      resnrm(i) = sqrt( zdotc(nbasis, Q(1,nstates+i),1, Q(1,nstates+i),1) )
      if(resnrm(i) < tolerance) nconv = nconv + 1
      if(resnrm(i) < tolerance/TFUDGE) ilock = ilock + 1
    enddo
    write(*,*) 'iter = ', iter, 'nconv = ', nconv
  
    if(nconv >= nstates-3) GOTO 10

    ! Apply preconditioner
    !do i=1,nstates
    !  W(:,i) = prec*rV(:,i)
    !enddo
  
    if(nlock > 0) then
      write(*,*) 'Warning: nlock=',nlock
    endif
  
!
! Apply Hamiltonian
!
    do i=1,nstates
      call apply_Ham(Q(1,nstates+i), HQ(1,nstates+i), nbasis, ik, rhoe, nnr) 
    enddo

    ! C <-- W* W
    call zgemm('C','N',nstates,nstates,nbasis, Z_ONE,Q(1,nstates+1),nbasis, &
         Q(1,nstates+1),nbasis, Z_ZERO,temp1,nstates)
    temp1  = (temp1 + conjg(temp1))*0.5d0
    
    ! Cholesky decomposition
    call zpotrf('U',nstates,temp1,nstates,info)
    if(info /= 0) then
      write(*,'(1x,A,I4)') 'ERROR calculating Cholesky decomposition : info ', info
      stop
    endif

    ! Solve linear equations
    ! W <-- W/C
    call ztrsm('R', 'U', 'N', 'N', nbasis,nstates, Z_ONE,temp1,nstates, Q(1,nstates+1),nbasis)
    ! HW <-- HW/C
    call ztrsm('R', 'U', 'N', 'N', nbasis,nstates, Z_ONE,temp1,nstates, HQ(1,nstates+1),nbasis)

    ! T <-- Q* HQ
    call zgemm('C','N',nstates3,nstates3,nbasis, Z_ONE,Q,nbasis, HQ,nbasis, &
      Z_ZERO,T,nstates3)
    T = ( T + conjg(T) )*0.5d0
    ! G <-- Q* Q
    call zgemm('C','N',nstates3,nstates3,nbasis, Z_ONE,Q,nbasis, Q,nbasis, &
      Z_ZERO,G,nstates3)
    G = ( G + conjg(G) )*0.5d0

    call eig_zhegv_f90(T,nstates3, G,nstates3, U,nstates3, nstates3)
    !call eig_zhegv(T,G,U,nstates3)

    ! X <-- Q U
    call zgemm('N','N',nbasis,nstates,nstates3, Z_ONE,Q,nbasis, U,nstates3, &
       Z_ZERO,tempX,nbasis)
    call zcopy(nbasis*nstates, tempX,1, Q,1)
    ! HX <-- HQ U
    call zgemm('N','N',nbasis,nstates,nstates3, Z_ONE,HQ,nbasis, U,nstates3, &
      Z_ZERO,tempX,nbasis)
    call zcopy(nbasis*nstates, tempX,1, HQ,1)

    ! P
    call zgemm('N','N',nbasis,nstates,nstates, Z_ONE,Q(1,nstates2+1),nbasis, &
      U(nstates2+1,1),nstates3, Z_ZERO,tempX,nbasis)
    call zgemm('N','N',nbasis,nstates,nstates, Z_ONE,Q(1,nstates+1),nbasis, &
      U(nstates+1,1),nstates3, Z_ZERO,Q(1,nstates2+1),nbasis)
    call zaxpy(nbasis*nstates, Z_ONE,tempX,1, Q(1,nstates2+1),1)

    ! HP
    call zgemm('N','N',nbasis,nstates,nstates, Z_ONE,HQ(1,nstates2+1),nbasis, &
      U(nstates2+1,1),nstates3, Z_ZERO,tempX,nbasis)
    call zgemm('N','N',nbasis,nstates,nstates, Z_ONE,HQ(1,nstates+1),nbasis, &
      U(nstates+1,1),nstates3, Z_ZERO,HQ(1,nstates2+1),nbasis)
    call zaxpy(nbasis*nstates, Z_ONE,tempX,1, HQ(1,nstates2+1),1)

    ! C = P* P
    call zgemm('C','N',nstates,nstates,nbasis, &
      Z_ONE,Q(1,nstates2+1),nbasis, Q(1,nstates2+1),nbasis, Z_ZERO,temp1,nstates)
    temp1 = ( temp1 + conjg(temp1) )*0.5d0
    !
    ! Cholesky decomposition
    call zpotrf('U',nstates,temp1,nstates,info)
    if(info /= 0) then
      write(*,'(1x,A,I4)') 'ERROR calculating Cholesky decomposition : info = ', info
      stop
    endif
    ! P = P/C
    call ztrsm('R','U','N','N', nbasis,nstates, Z_ONE,temp1,nstates, &
      Q(1,nstates2+1),nbasis)
    ! HP = HP/C
    call ztrsm('R','U','N','N', nbasis,nstates, Z_ONE,temp1,nstates, &
      HQ(1,nstates2+1),nbasis)
  enddo ! end of iteration loop

10 continue

  ! XHX = X* HX
  !call zgemm('C','N',nstates,nstates,nbasis,ONE,X,nbasis,HX,nbasis,ZERO,XHX_temp,nstates)
  ! XHX = (XHX + XHX*)/2
  !call mkl_zomatadd('Col','N','C',nstates,nstates,HALF,XHX_temp,nstates,&
  !    HALF,XHX_temp,nstates, XHX,nstates)
  ! Calculate the eigenvalues and eigenvectors
  !call eig_zheev(XHX,lambda,nstates)
  !W = X ! save X to W, W must not be used again ...
  !call zgemm('N','N',nbasis,nstates,nstates, ONE,W,nbasis, XHX,nbasis, ZERO,X,nbasis)

  call mytimer(time2)
  write(*,*) 'Number of converged eigenvalues:', nconv
  write(*,*) 'Time for LOBPCG algorithm: ', time2-time1
  do i=1,nstates
    write(*,*) i, lambda(i), resnrm(i)
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



!-------------------------------------------
subroutine ortho_qr_f90(X, nbasis, nstates)
!-------------------------------------------
  implicit none
  ! Argument
  integer :: nbasis, nstates
  complex(8) :: X(nbasis,nstates)
  ! Local
  integer :: m, n, k, lwork, info
  complex(8), allocatable :: tau(:), work(:)

  m = nbasis
  n = nstates
  k = nstates
  lwork = n

  allocate(tau(k))
  allocate(work(lwork))

  ! Compute QR factorization
  call zgeqrf(m,n,X,m,tau,work,lwork,info)
  if(info /= 0) then
    write(*,'(1x,A,I4)') 'ERROR calling ZGEQRF : info = ', info
    return
  endif

  ! Construct Q explicitly
  call zungqr(m,k,k,X,m,tau,work,lwork,info)
  if(info /= 0) then
    write(*,'(1x,A,I4)') 'ERROR calling ZUNGQR : info = ', info
    return
  endif

  deallocate(tau)
  deallocate(work)
end subroutine

!----------------------------------------
subroutine eig_zheevd_f90(A,lda, lambda,N)
!----------------------------------------
  implicit none
  ! Arguments
  integer :: ldA,N
  real(8) :: lambda(N)
  complex(8) :: A(ldA,N)
  ! Local
  integer :: lwork, lrwork, liwork
  complex(8), allocatable :: work(:)
  real(8), allocatable :: rwork(:)
  integer, allocatable :: iwork(:)
  integer :: info

  lwork = N*N + 2*N
  lrwork = 2*N*N + 5*N + 1
  liwork = 5*N + 3
  
  allocate(work(lwork))
  allocate(rwork(lrwork))
  allocate(iwork(liwork))

  write(*,*) 'sum(A) = ', sum(A)
  call zheevd('V','U', N, A,ldA, lambda, work,lwork, rwork,lrwork, iwork,liwork, info)
  if(info /= 0) then
    write(*,'(1x,A,I4)') 'ERROR calling ZHEEVD : info = ', info
    stop
  endif

  deallocate(work)
  deallocate(rwork)
  deallocate(iwork)
end subroutine


!---------------------------------------------------
subroutine eig_zhegv_f90(A,ldA, B,ldB, evec,lde, N)
!---------------------------------------------------
  implicit none
  ! Arguments
  integer :: ldA,ldB,lde,N
  complex(8) :: A(ldA,N)
  complex(8) :: B(ldB,N)
  complex(8) :: evec(ldE,N)
  ! Local
  real(8), parameter :: SMALL=1.d-10
  complex(8), parameter :: Z_ZERO=(0.d0,0.d0), Z_ONE=(1.d0,0.d0)
  integer :: lwork, lrwork, liwork, info, i, nn
  complex(8), allocatable :: work(:)
  real(8), allocatable :: rwork(:)
  integer, allocatable :: iwork(:)
  real(8), allocatable :: eval(:)
  real(8) :: scal

  lwork = N*N + 2*N
  lrwork = 2*N*N + 5*N + 1
  liwork = 5*N + 3

  allocate(work(lwork))
  allocate(rwork(lrwork))
  allocate(iwork(liwork))
  allocate(eval(N))

  ! Diagonalize B
  call zheevd('V','U',N, B,ldB, eval, work,lwork, rwork,lrwork, iwork,liwork, info)
  if(info /= 0) then
    write(*,'(1x,A,I4)') 'ERROR calling ZHEEVD in EIG_ZHEGV_F90 : info = ', info
    stop
  endif

  nn = 0
  do i=1,N
    if( eval(i) > SMALL ) then
      nn = nn + 1
      scal = 1.d0/sqrt(eval(i))
      call zdscal(N, scal, B(:,i),1)
    endif
  enddo
  if(nn < N) then
    write(*,'(1x,A,I4)') 'ERROR: Number of linearly independent vectors = ', nn
    write(*,'(1x,A,I4)') '       while size of the problem = ', N
    stop
  endif

  ! Transform A:
  ! A <-- evec(B)* A evec(B)
  call zgemm('N','N', N,N,N, Z_ONE,A,ldA, B,ldB, Z_ZERO,evec,lde)
  call zgemm('C','N', N,N,N, Z_ONE,B,ldB, evec,lde, Z_ZERO,A,ldA)

  ! Diagonalize transformed A
  call zheevd('V','U',N, A,ldA, eval, work,lwork, rwork,lrwork, iwork,liwork, info)
  if(info /= 0) then
    write(*,'(1x,A,I4)') 'ERROR calling ZHEEVD in EIG_ZHEEVD_F90 : info = ', info
    stop
  endif

  ! Back transform eigenvectors
  call zgemm('N','N', N,N,N, Z_ONE,B,ldB, A,ldA, Z_ZERO,evec,lde)

  write(*,*) 'sum(eval) = ', sum(eval)

  deallocate(work)
  deallocate(rwork)
  deallocate(iwork)
  deallocate(eval)
end subroutine

! Generalized eigenvalue problem, Hermitian A, symmetric, positive definite B
!----------------------------------------------
subroutine eig_zhegv(A,B,Rvec,N)
!----------------------------------------------
  implicit none
  ! Arguments
  integer :: N
  complex(8) :: A(N,N)
  complex(8) :: B(N,N)
  complex(8) :: Rvec(N,N) ! eigenvectors
  ! Parameter
  integer, parameter :: NB=64
  ! Local variables
  complex(8) :: alpha(N), beta(N)
  integer :: lwork
  complex(8), allocatable :: work(:)
  real(8), allocatable :: rwork(:)
  real(8), allocatable :: eigval(:)
  integer :: info
  integer :: i

  lwork = (NB+1)*N
  allocate(eigval(N))
  allocate(work(lwork))
  allocate(rwork(8*N))
  Rvec = A ! Save A so that it will not be modified

  ! Solve the generalized eigenvalue problem
  ! Calculate right eigenvectors
  call zhegv(1,'Vectors','Upper',N,Rvec,N,B,N,eigval,work,LWORK,rwork,info)
  if(info /= 0) then
    write(*,*) 'Error in calling zhegv: info = ', info
    stop
  endif

  write(*,*) 'sum(eigval) = ', sum(eigval)

  deallocate(eigval)
  deallocate(work)
  deallocate(rwork)
end subroutine

