module LinearAlgebra

! Interface blocks
interface chol
  module procedure chol_dpotrf, chol_zpotrf
end interface

interface chol_linsolve_t
  module procedure c8_chol_linsolve_t, c8_r8_chol_linsolve_t
end interface

interface forsubst_t
  module procedure c8_r8_forsubs_t, c8_c8_forsubs_t
end interface

!============================
contains
!============================

! For emulating:
!   [A,eigval] = eig(A)
!   eigval = diag(eigval)
!----------------------------------------
subroutine eig_zheev(A,eigval,dimA)
!----------------------------------------
  implicit none
  ! Arguments
  integer :: dimA
  complex(8) :: A(dimA,dimA)
  real(8) :: eigval(dimA)
  ! Workspace array for ZHEEV
  integer :: lwork
  complex(8), allocatable :: work(:)
  complex(8), allocatable :: rwork(:)
  ! Local variables
  integer :: info

  lwork = 2*dimA-1
  allocate(work(lwork))
  allocate(rwork(3*dimA-2))

  call zheev('v', 'u', dimA, A, dimA, eigval, work, lwork, rwork, info)
  if(info /= 0) then
    write(*,*) 'Error on calling zheev: info=',info
    stop
  endif

  ! Free memory
  deallocate(work)
  deallocate(rwork)
end subroutine

! Only diagonalize matrix A
!----------------------------------------
subroutine diagon_zheev(A,dimA)
!----------------------------------------
  implicit none
  ! Arguments
  integer :: dimA
  complex(8) :: A(dimA,dimA)
  ! Workspace array for ZHEEV
  integer :: lwork
  complex(8), allocatable :: work(:)
  complex(8), allocatable :: rwork(:)
  real(8) :: eigval(dimA) ! use allocatable?
  ! Local variables
  integer :: info

  lwork = 2*dimA-1
  allocate(work(lwork))
  allocate(rwork(3*dimA-2))

  call zheev('v', 'u', dimA, A, dimA, eigval, work, lwork, rwork, info)
  if(info /= 0) then
    write(*,*) 'Error on calling zheev: info=',info
    stop
  endif

  ! Free memory
  deallocate(work)
  deallocate(rwork)
end subroutine


! Cholesky decomposition, real(8) version
!------------------------------------------
subroutine chol_dpotrf(A,dimA,cholFlag)
!------------------------------------------
  implicit none
  ! Argument
  integer :: dimA
  real(8) :: A(dimA,dimA)
  integer :: cholFlag
  ! Iterators
  integer :: i,j

  ! Use upper tridiagonal
  write(*,*) 'Calling dpotrf'
  call dpotrf('U',dimA,A,dimA,cholFlag)
  if(cholFlag /= 0) then
    write(*,*) 'Error in calling dpotrf: info = ', cholFlag
    stop
  endif
  ! Zero out lower diagonal
  do j=1,dimA
  do i=j+1,dimA
    A(i,j) = 0.d0
  enddo
  enddo
end subroutine

! Cholesky decomposition, complex(8) version
!--------------------------------------
subroutine chol_zpotrf(A,dimA,cholFlag)
!--------------------------------------
  implicit none
  ! Argument
  integer :: dimA
  complex(8) :: A(dimA,dimA)
  integer :: cholFlag
  ! Iterators
  integer :: i,j

  ! Use upper diagonal
  call zpotrf('U',dimA,A,dimA,cholFlag)
  if(cholFlag /= 0) then
    write(*,*) 'Error in calling zpotrf: info = ', cholFlag
    stop
  endif
  ! Zero out lower diagonal
  do j=1,dimA
  do i=j+1,dimA
    A(i,j) = cmplx(0.d0,0.d0)
  enddo
  enddo
end subroutine

! From Burden, et. al.
subroutine r8_chol_linsolve(U,n,nrhs,b,x)
  implicit none
  ! arguments
  integer :: n
  integer :: nrhs
  real(8) :: U(n,n)
  real(8) :: b(n,nrhs)
  real(8) :: x(n,nrhs)
  ! local variables
  integer :: i,j,k
  real(8) :: sum1(nrhs)

  x(1,:) = b(1,:)/U(1,1)
  do i=2,n
    sum1(:) = 0.d0
    do j=1,i-1
      sum1(:) = sum1(:) + U(j,i)*x(j,:)
    enddo
    x(i,:) = (b(i,:) - sum1(:))/U(i,i)
  enddo
  x(n,:) = x(n,:)/U(n,n)

  do i=n-1,1,-1
    sum1(:) = 0.d0
    do j=i+1,n
      sum1(:) = sum1(:) + U(i,j)*x(j,:)
    enddo
    x(i,:) = (x(i,:) - sum1(:))/U(i,i)
  enddo
end subroutine

! From Burden, et. al.
subroutine c8_chol_linsolve(U,n,nrhs,b,x)
  implicit none
  ! arguments
  integer :: n
  integer :: nrhs
  complex(8) :: U(n,n)
  complex(8) :: b(n,nrhs)
  complex(8) :: x(n,nrhs)
  ! local variables
  integer :: i,j,k
  complex(8) :: sum1(nrhs)

  x(1,:) = b(1,:)/U(1,1)
  do i=2,n
    sum1(:) = 0.d0
    do j=1,i-1
      sum1(:) = sum1(:) + conjg(U(j,i))*x(j,:)
    enddo
    x(i,:) = (b(i,:) - sum1(:))/U(i,i)
  enddo
  x(n,:) = x(n,:)/U(n,n)

  do i=n-1,1,-1
    sum1(:) = 0.d0
    do j=i+1,n
      sum1(:) = sum1(:) + U(i,j)*x(j,:)
    enddo
    x(i,:) = (x(i,:) - sum1(:))/U(i,i)
  enddo
end subroutine

! b = A/b
subroutine c8_chol_linsolve_t(U,n,nrhs,b)
  implicit none
  ! arguments
  integer :: n
  integer :: nrhs
  complex(8) :: U(n,n)
  complex(8) :: b(nrhs,n)
  ! local variables
  integer :: i,j,k
  complex(8) :: sum1(nrhs)

  write(*,*) 'Calling c8_chol_linsolve_t'

  b(:,1) = b(:,1)/U(1,1)
  do i=2,n
    sum1(:) = 0.d0
    do j=1,i-1
      sum1(:) = sum1(:) + conjg(U(j,i))*b(:,j)
    enddo
    b(:,i) = (b(:,i) - sum1(:))/U(i,i)
  enddo
  b(:,n) = b(:,n)/U(n,n)

  do i=n-1,1,-1
    sum1(:) = 0.d0
    do j=i+1,n
      sum1(:) = sum1(:) + U(j,i)*b(:,j)
    enddo
    b(:,i) = (b(:,i) - sum1(:))/U(i,i)
  enddo
end subroutine

! A = A/b
! L is similar to b
subroutine c8_r8_chol_linsolve_t(U,n,nrhs,b)
  implicit none
  ! arguments
  integer :: n
  integer :: nrhs
  real(8) :: U(n,n)
  complex(8) :: b(nrhs,n)
  ! local variables
  integer :: i,j,k
  complex(8) :: sum1(nrhs)

  write(*,*) 'Calling c8_r8_chol_linsolve_t'
  b(:,1) = b(:,1)/U(1,1)
  do i=2,n
    sum1(:) = 0.d0
    do j=1,i-1
      sum1(:) = sum1(:) + U(j,i)*b(:,j)
    enddo
    b(:,i) = (b(:,i) - sum1(:))/U(i,i)
  enddo
  b(:,n) = b(:,n)/U(n,n)

  do i=n-1,1,-1
    sum1(:) = 0.d0
    do j=i+1,n
      sum1(:) = sum1(:) + U(i,j)*b(:,j)
    enddo
    b(:,i) = (b(:,i) - sum1(:))/U(i,i)
  enddo
end subroutine

! Backsubtitution of lower triangular matrix L
! L is obtained from transpose(U)
!   L(i,j) = U(j,i)
! B wil be overwritten with the solution
subroutine c8_r8_forsubs_t(L,dimL,nrhs,B)
  implicit none
  ! Arguments
  integer :: dimL
  integer :: nrhs
  real(8) :: L(dimL,dimL)
  complex(8) :: B(nrhs,dimL)
  ! Local variables
  integer :: i,j
  complex(8) :: sum1(nrhs)

  B(:,1) = B(:,1)/L(1,1) ! From row 1
  ! Loop over row, start from row 2 till N
  do j=2,dimL
    sum1(:) = cmplx(0.d0,0.d0)
    ! Loop over column 1 to m-1
    do i=1,j-1
      sum1(:) = sum1(:) + L(i,j)*B(:,i)
    enddo
    B(:,j) = ( B(:,j) - sum1(:) )/L(j,j)
  enddo
end subroutine

!------------------------------------------
subroutine c8_c8_forsubs_t(L,dimL,nrhs,B)
!------------------------------------------
  implicit none
  ! Arguments
  integer :: dimL
  integer :: nrhs
  complex(8) :: L(dimL,dimL)
  complex(8) :: B(nrhs,dimL)
  ! Local variables
  integer :: i,j
  complex(8) :: sum1(nrhs)

  B(:,1) = B(:,1)/L(1,1) ! From row 1
  ! Loop over row, start from row 2 till N
  do j=2,dimL
    sum1(:) = cmplx(0.d0,0.d0)
    ! Loop over column 1 to j-1
    do i=1,j-1
      sum1(:) = sum1(:) + L(i,j)*B(:,i) ! use conjg or not?
    enddo
    B(:,j) = ( B(:,j) - sum1(:) )/L(j,j)
  enddo
end subroutine

!--------------------------
subroutine qr_decomp(A,m,n)
!--------------------------
  implicit none
  ! Arguments
  integer :: m, n
  complex(8) :: A(m,n)
  ! Local variables
  integer :: n2, k, lwork
  complex(8), allocatable :: A2(:,:)
  complex(8), allocatable :: tau(:), work(:)
  integer :: i, info, lda

  n2 = m
  k = min(m,n2)
  lwork = n2

  allocate(tau(k))
  allocate(work(lwork))

  ! Compute the QR factorization
  lda = m
  call zgeqrf(m,n2,A,lda,tau,work,lwork,info)
  if(info /= 0) then
    write(*,*)
    write(*,*) 'zgeqrf returned nonzero info = ', info
    return
  endif

  ! Construct Q=A explicitly
  call zungqr(m,k,k,A,lda,tau,work,lwork,info)
  if(info /= 0) then
    write(*,*)
    write(*,*) 'zungqr returned nonzero info = ', info
    return
  endif

  deallocate(tau)
  deallocate(work)
end subroutine qr_decomp

! 'Economy size QR decomposition
!---------------------------
subroutine qr_decomp0(A,Q,R)
!---------------------------
  implicit none
  ! Arguments
  complex(8) :: A(:,:) ! size m x n
  complex(8) :: Q(:,:) ! size m x k, with k = min(m,n)
  complex(8) :: R(:,:) ! size k x n
  ! Local variables
  integer :: m, n, k
  integer :: lwork, i, info, lda
  complex(8), allocatable :: tau(:), work(:)
  
  m = size(A,1)
  n = size(A,2)
  k = min(m,n)
  lwork = n

  ! Allocate arrays
  allocate(tau(k))
  allocate(work(lwork))
  
  ! Compute QR factorization
  lda = m
  call zgeqrf(m,n,A,lda,tau,work,lwork,info)
  if(info /= 0) then
    write(*,*) 'zgeqrf returned nonzero info = ', info
    return
  endif

  R(1:k,1:n) = 0d0
  do i=1,k
    R(i,i:n) = A(i,i:n)
  enddo

  ! Construct Q explicitly
  Q(1:m,1:k) = A(1:m,1:k)
  call zungqr(m,k,k,Q,lda,tau,work,lwork,info)
  if(info /= 0) then
    write(*,*) 'zungqr returned nonzero info = ', info
    return
  endif

  deallocate(tau)
  deallocate(work)
end subroutine qr_decomp0

!------------------------------
subroutine r8_inverse1(A)
!------------------------------
  implicit none
  ! Local variables
  real(8) :: A(:,:)
  ! Local variables
  integer :: dim1,dim2
  integer :: lwork
  real(8), allocatable :: work(:)
  integer, allocatable :: ipiv(:)
  integer :: lda
  integer :: info

  ! Check dimension
  dim1 = size(A,1); dim2 = size(A,2)
  if(dim1 <= 0 .OR. dim2 <= 0) then
    write(*,*) 'error, dim1=',dim1,' dim2=',dim2
    stop
  endif
  if(dim1 /= dim2) then
    write(*,*) 'error in r8_inverse1: dim1 /= dim2 ', dim1,dim2
    stop
  endif
  
  ! Factorize A
  lda = dim1
  allocate(ipiv(dim1)) ! allocate ipiv first
  call dgetrf(dim1,dim2,A,lda,ipiv,info)
  if(info /= 0) then
    write(*,*) 'Error in r8_inverse: dgetrf returned info = ', info
    stop
  endif

  ! Calculate inverse of A
  lwork = 2*dim1
  allocate(work(lwork))
  call dgetri(dim1,A,lda,ipiv,work,lwork,info)
  if(info /= 0) then
    write(*,*) 'Error in r8_inverse: dgetri returned info = ', info
    stop
  endif

  ! Free memory
  deallocate(ipiv)
  deallocate(work)
end subroutine

!------------------------------
subroutine c8_inverse(A)
!------------------------------
  implicit none
  ! Local variables
  complex(8) :: A(:,:)
  ! Local variables
  integer :: dim1,dim2
  integer :: lwork
  complex(8), allocatable :: work(:)
  integer, allocatable :: ipiv(:)
  integer :: lda
  integer :: info

  ! Check dimension
  dim1 = size(A,1); dim2 = size(A,2)
  if(dim1 <= 0 .OR. dim2 <= 0) then
    write(*,*) 'error, dim1=',dim1,' dim2=',dim2
    stop
  endif
  if(dim1 /= dim2) then
    write(*,*) 'error in c8_inverse: dim1 /= dim2 ', dim1,dim2
    stop
  endif
  
  ! Factorize A
  lda = dim1
  allocate(ipiv(dim1)) ! allocate ipiv first
  call zgetrf(dim1,dim2,A,lda,ipiv,info)
  if(info /= 0) then
    write(*,*) 'Error in c8_inverse: zgetrf returned info = ', info
    stop
  endif

  ! Calculate inverse of A
  lwork = 2*dim1
  allocate(work(lwork))
  call zgetri(dim1,A,lda,ipiv,work,lwork,info)
  if(info /= 0) then
    write(*,*) 'Error in c8_inverse: zgetri returned info = ', info
    stop
  endif

  ! Free memory
  deallocate(ipiv)
  deallocate(work)
end subroutine

! Generalized eigenvalue problem, nonsymmetric
!----------------------------------------------
subroutine eig_zggev(A,B,eigval,Rvec,N)
!----------------------------------------------
  implicit none
  ! Arguments
  integer :: N
  complex(8) :: A(N,N)
  complex(8) :: B(N,N)
  complex(8) :: eigval(N) ! eigenvalues
  complex(8) :: Rvec(N,N) ! eigenvectors
  ! Local variables
  complex(8) :: alpha(N), beta(N)
  integer :: lwork
  complex(8) :: dummy(1,1) 
  complex(8), allocatable :: work(:)
  real(8) :: rwork(8*N)
  integer :: info
  integer :: i

  lwork = N + 64*N
  allocate(work(lwork))

  ! Solve the generalized eigenvalue problem
  ! Calculate right eigenvectors
  call zggev('N','V',N,A,N,B,N,alpha,beta,dummy,1,Rvec,N,work,lwork,rwork,info)
  if(info /= 0) then
    write(*,*) 'Error in calling zggev: info = ', info
    stop
  endif
  ! Calculate eigenvalues
  do i=1,N
    eigval(i) = alpha(i)/beta(i)
  enddo

  deallocate(work)
end subroutine

! Generalized eigenvalue problem, Hermitian A, symmetric, positive definite B
!----------------------------------------------
subroutine eig_zhegv(A,B,eigval,Rvec,N)
!----------------------------------------------
  implicit none
  ! Arguments
  integer :: N
  complex(8) :: A(N,N)
  complex(8) :: B(N,N)
  real(8) :: eigval(N) ! eigenvalues
  complex(8) :: Rvec(N,N) ! eigenvectors
  ! Parameter
  integer, parameter :: NB=64
  ! Local variables
  complex(8) :: alpha(N), beta(N)
  integer :: lwork
  complex(8), allocatable :: work(:)
  real(8) :: rwork(8*N) ! not allocatable?
  integer :: info
  integer :: i

  lwork = (NB+1)*N
  allocate(work(lwork))
  Rvec = A ! Save A so that it will not be modified

  ! Solve the generalized eigenvalue problem
  ! Calculate right eigenvectors
  call zhegv(1,'Vectors','Upper',N,Rvec,N,B,N,eigval,work,LWORK,rwork,info)
  if(info /= 0) then
    write(*,*) 'Error in calling zhegv: info = ', info
    stop
  endif

  deallocate(work)
end subroutine

end module

