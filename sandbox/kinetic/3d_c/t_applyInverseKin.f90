! efefer, 20 Jan 2016
! test preconditioning matrix as described in JCP-119-8842

!------------------------------------------------------------------------------
PROGRAM t_precond
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF3d
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: N = 3
  REAL(8), PARAMETER :: L = 2.d0
  !
  TYPE(LF3d_t) :: LF
  !
  INTEGER :: ip1, ip2, i1,i2,j1,j2,k1,k2
  REAL(8) :: PrecInv(N**3,N**3)
  REAL(8) :: Prec(N**3,N**3), vecIn(N**3), vecOut(N**3)
  REAL(8) :: vecOutIter(N**3)
  
  ! Initialize the basis functions
  CALL init_LF3d_c( LF, (/N,N,N/), (/0.d0,0.d0,0.d0/), (/L,L,L/) )
 
  ! Model Hamiltonian Kin + 1
  DO ip2 = 1, N**3
    DO ip1 = ip2, N**3
      i1 = LF%lin2xyz(1,ip1)
      i2 = LF%lin2xyz(1,ip2)
      !
      j1 = LF%lin2xyz(2,ip1)
      j2 = LF%lin2xyz(2,ip2)
      !
      k1 = LF%lin2xyz(3,ip1)
      k2 = LF%lin2xyz(3,ip2)
      !
      PrecInv(ip1,ip2) = 0.d0
      !
      IF( j1 == j2 .AND. k1 == k2 ) THEN
        PrecInv(ip1,ip2) = PrecInv(ip1,ip2) + LF%LFx%D2jl(i1,i2)
      ENDIF
      !
      IF( i1 == i2 .AND. k1 == k2 ) THEN
        PrecInv(ip1,ip2) = PrecInv(ip1,ip2) + LF%LFy%D2jl(j1,j2)
      ENDIF
      !
      IF( i1 == i2 .AND. j1 == j2 ) THEN
        PrecInv(ip1,ip2) = PrecInv(ip1,ip2) + LF%LFz%D2jl(k1,k2)
      ENDIF
      !
      PrecInv(ip1,ip2) = -0.5d0*PrecInv(ip1,ip2)  ! Kinetic
      PrecInv(ip2,ip1) = PrecInv(ip1,ip2)
      !WRITE(*,'(2I8,F18.10)') ip1, ip2, Hamiltonian(ip1,ip2)
    ENDDO
    PrecInv(ip2,ip2) = PrecInv(ip2,ip2) + 1.d0 ! diagonal only
  ENDDO
  !

  Prec = PrecInv
  CALL r8_inverse1(Prec)

  ! Test apply_invPrec
  !vecIn(:) = 1.d0
  !vecOut = matmul( PrecInv, vecIn )
  !CALL apply_invPrec( N**3, LF, vecIn, vecOutIter )
  ! This value should be very small
  !WRITE(*,*) sum( (vecOut - vecOutIter)**2 )

  ! Test apply_Prec or CG
  vecIn(:) = 2.d0
  vecOut = matmul( Prec, vecIn )
  CALL apply_PrecCG( N**3, LF, vecIn, vecOutIter )
  !CALL apply_invPrec( N**3, LF, vecOutIter, vecOut )
  ! This value should be very small
  WRITE(*,*) sum( (vecOut - vecOutIter)**2 )
  !WRITE(*,*) sum( (vecOut - vecIn)**2 )
  WRITE(*,*) sum(vecOut)

CONTAINS


SUBROUTINE apply_invPrec( N, LF, veci, veco )
  IMPLICIT NONE
  TYPE(LF3d_t) :: LF
  INTEGER :: N
  REAL(8) :: veci(N), veco(N)
  !
  INTEGER :: i, j, k, ip, ii, jj, kk, Nx, Ny, Nz

  Nx = LF%LFx%N
  Ny = LF%LFy%N
  Nz = LF%LFz%N
  !
  DO ip=1,N
    i = LF%lin2xyz(1,ip)
    j = LF%lin2xyz(2,ip)
    k = LF%lin2xyz(3,ip)
    !
    veco(ip) = veci(ip)  ! due to unit matrix
    !
    DO ii=1,Nx
      veco(ip) = veco(ip) + -0.5d0*LF%LFx%D2jl(ii,i)*veci(LF%xyz2lin(ii,j,k))
    ENDDO
    !
    DO jj=1,Ny
      veco(ip) = veco(ip) + -0.5d0*LF%LFy%D2jl(jj,j)*veci(LF%xyz2lin(i,jj,k))
    ENDDO
    !
    DO kk=1,Nz
      veco(ip) = veco(ip) + -0.5d0*LF%LFz%D2jl(kk,k)*veci(LF%xyz2lin(i,j,kk))
    ENDDO
  ENDDO

END SUBROUTINE


SUBROUTINE apply_PrecCG( N, LF, rho, phi )
  IMPLICIT NONE
  INTEGER :: N
  TYPE(LF3d_t) :: LF
  REAL(8) :: rho(N), phi(N)
  !
  REAL(8), ALLOCATABLE :: r(:), p(:) ! residual
  REAL(8), ALLOCATABLE :: nabla2(:)
  INTEGER :: ip, iter
  REAL(8) :: rsold, rsnew, alpha

  ALLOCATE( r(N), p(N), nabla2(N) )

  !
  DO ip=1,N
    vecOut(ip) = 0.d0
  ENDDO

  CALL apply_invPrec( N, LF, phi, nabla2 )
  r(:) = rho(:) - nabla2(:)
  p(:) = r(:)

  !
  rsold = dot_product(r,r)
  WRITE(*,*) 'rsold = ', rsold

  DO iter=1,N
    CALL apply_invPrec( N, LF, p, nabla2 )
    !
    alpha = rsold/dot_product(p,nabla2)
    !
    phi(:) = phi(:) + alpha*p(:)
    !
    r(:) = r(:) - alpha*nabla2(:)
    !
    rsnew = dot_product(r,r)
    WRITE(*,*) 'rsnew = ', sqrt(rsnew)
    !
    IF(sqrt(rsnew) < 1.d-10) THEN
    !IF(rsnew < 1.d-10) THEN
      WRITE(*,*) 'Convergence in apply_PrecCG in iter:', iter
      EXIT
    ENDIF
    p(:) = r(:) + (rsnew/rsold)*p(:)
    rsold = rsnew
  ENDDO

  DEALLOCATE( r, p, nabla2 )
END SUBROUTINE



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


END PROGRAM


