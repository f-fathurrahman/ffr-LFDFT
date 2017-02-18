!------------------------------------------------------------------------------
PROGRAM t_fun
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
  REAL(8) :: evalsK(N**3), Kx,Ky,Kz
  INTEGER :: idxdum(N**3)
  
  ! Initialize the basis functions
  CALL init_LF3d_c( LF, (/N,N,N/), (/0.d0,0.d0,0.d0/), (/L,L,L/) )
 
  ! Kinetic part
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
      Hamiltonian(ip2,ip1) = Hamiltonian(ip1,ip2)
      WRITE(*,'(2I8,F18.10)') ip1, ip2, Hamiltonian(ip1,ip2)
    ENDDO
  ENDDO
  !
  Hamiltonian(:,:) = -0.5d0*Hamiltonian(:,:)
  evecsT(:,:) = Hamiltonian(:,:)

  CALL eig_dsyev(evecsT,evalsT,N**3)

  WRITE(*,*)
  WRITE(*,*) 'Eigenvalues of T:'
  DO ip=1,N**3
    i = LF%lin2xyz(1,ip)
    j = LF%lin2xyz(2,ip)
    k = LF%lin2xyz(3,ip)
    WRITE(*,'(1x,4I5,F18.10)') ip, i, j, k, evalsT(ip)
  ENDDO
  
  WRITE(*,*)
  WRITE(*,*) 'Eigenvalues based on Kvec:'
  ip = 0
  DO k = 1,N
    DO j = 1,N
      DO i = 1,N
        Kx = i*PI/L
        Ky = j*PI/L
        Kz = k*PI/L
        ip = ip + 1
        evalsK(ip) = 0.5d0*( Kx**2 + Ky**2 + Kz**2 )
        idxdum(ip) = ip
      ENDDO
    ENDDO
  ENDDO
  CALL hpsort( N**3, evalsT, idxdum )

  DO ip=1,N**3
    WRITE(*,'(1x,I5,F18.10)') ip, evalsT(ip) 
  ENDDO
  !

END PROGRAM


!----------------------------------
SUBROUTINE eig_dsyev(A,eigval,dimA)
!----------------------------------
  IMPLICIT NONE
  ! Arguments
  INTEGER :: dimA
  REAL(8) :: A(dimA,dimA)
  REAL(8) :: eigval(dimA)
  ! Workspace array for ZHEEV
  INTEGER :: lwork
  REAL(8), ALLOCATABLE :: work(:)
  ! Local variables
  INTEGER :: info

  lwork = 3*dimA-1
  ALLOCATE(work(lwork))

  CALL dsyev('v', 'u', dimA, A, dimA, eigval, work, lwork, info)
  IF(info /= 0) THEN
    WRITE(*,*) 'Error on calling dsyev: info=',info
    STOP
  ENDIF

  ! Free memory
  DEALLOCATE(work)
END SUBROUTINE
