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
  INTEGER :: ip, ip1, ip2, i1,i2,j1,j2,k1,k2, i,j,k
  REAL(8) :: evecsT(N**3,N**3), evalsT(N**3), Hamiltonian(N**3,N**3)
  REAL(8) :: invHam(N**3,N**3), tMat(N**3,N**3)
  
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
      Hamiltonian(ip1,ip2) = 0.d0
      !
      IF( j1 == j2 .AND. k1 == k2 ) THEN
        Hamiltonian(ip1,ip2) = Hamiltonian(ip1,ip2) + LF%LFx%D2jl(i1,i2)
      ENDIF
      !
      IF( i1 == i2 .AND. k1 == k2 ) THEN
        Hamiltonian(ip1,ip2) = Hamiltonian(ip1,ip2) + LF%LFy%D2jl(j1,j2)
      ENDIF
      !
      IF( i1 == i2 .AND. j1 == j2 ) THEN
        Hamiltonian(ip1,ip2) = Hamiltonian(ip1,ip2) + LF%LFz%D2jl(k1,k2)
      ENDIF
      !
      Hamiltonian(ip1,ip2) = -0.5d0*Hamiltonian(ip1,ip2)  ! Kinetic
      Hamiltonian(ip2,ip1) = Hamiltonian(ip1,ip2)
      !WRITE(*,'(2I8,F18.10)') ip1, ip2, Hamiltonian(ip1,ip2)
    ENDDO
    Hamiltonian(ip2,ip2) = Hamiltonian(ip2,ip2) + 1.d0 ! diagonal only
  ENDDO
  !
  !evecsT(:,:) = Hamiltonian(:,:)

  !CALL write_matrix_hdf5( Hamiltonian, N**3, N**3, 'Precond.h5', 'precond' )

!  CALL eig_dsyev(evecsT,evalsT,N**3)
!
!  WRITE(*,*)
!  WRITE(*,*) 'Eigenvalues of T:'
!  DO ip=1,N**3
!    i = LF%lin2xyz(1,ip)
!    j = LF%lin2xyz(2,ip)
!    k = LF%lin2xyz(3,ip)
!    WRITE(*,'(1x,4I5,F18.10)') ip, i, j, k, evalsT(ip)
!  ENDDO

  invHam = Hamiltonian
  CALL r8_inverse1(invHam)
  tMat = matmul( Hamiltonian, invHam )

  CALL write_matrix_hdf5( invHam, N**3, N**3, 'invHam.h5', 'inverse' ) 
  CALL write_matrix_hdf5( tMat, N**3, N**3, 'tMat.h5', 'test_inverse' ) 


CONTAINS

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


