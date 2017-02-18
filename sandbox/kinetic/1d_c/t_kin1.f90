!------------------------------------------------------------------------------
PROGRAM t_fun
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF1d
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: N = 7
  REAL(8), PARAMETER :: L = 2.d0
  !
  TYPE(LF1d_t) :: LF
  !
  INTEGER :: ip, i
  REAL(8) :: evecsT(N,N), evalsT(N)
  REAL(8) :: evalsK(N), Kvec
  
  ! Initialize the basis functions
  CALL init_LF1d_c( LF, N, -L/2, L/2 )
  CALL info_LF1d( LF, .TRUE. )

  evecsT(:,:) = -0.5d0*LF%D2jl(:,:)
  CALL eig_dsyev(evecsT,evalsT,N)

  WRITE(*,*)
  WRITE(*,*) 'Eigenvalues of T:'
  DO ip=1,N
    WRITE(*,'(1x,I5,F18.10)') ip, evalsT(ip)
  ENDDO
  
  DO i=1,N
    Kvec = i*PI/L
    evalsK(i) = 0.5d0*Kvec**2
  ENDDO
  !
  WRITE(*,*)
  WRITE(*,*) 'Eigenvalues based on Kvec:'
  DO ip=1,N
    WRITE(*,'(1x,I5,F18.10)') ip, evalsK(ip)
  ENDDO
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
