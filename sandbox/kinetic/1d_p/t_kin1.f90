!------------------------------------------------------------------------------
PROGRAM t_fun
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF1d
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: N = 5
  REAL(8), PARAMETER :: L = 2.d0
  !
  TYPE(LF1d_t) :: LF
  !
  INTEGER :: ip, i
  REAL(8) :: evecsT(N,N), evalsT(N)
  REAL(8) :: evalsG(N), Gvec
  
  ! Initialize the basis functions
  CALL init_LF1d_p( LF, N, -L/2, L/2 )
  CALL info_LF1d( LF, .TRUE. )

  evecsT(:,:) = -0.5d0*LF%D2jl(:,:)
  CALL eig_dsyev(evecsT,evalsT,N)

  WRITE(*,*)
  WRITE(*,*) 'Eigenvalues of T:'
  DO ip=1,N
    WRITE(*,'(1x,I5,F18.10)') ip, evalsT(ip)
  ENDDO
  
  ip = 1
  DO i=0,N/2
    Gvec = 2.d0*PI/L * i
    IF( i == 0 ) THEN
      evalsG(i) = 0.d0
      ip = ip + 1 ! ip will be 2
    ELSE
      evalsG(ip) = 0.5d0*Gvec**2  ! remember the (imag)**2 factor in 2nd deriv of PW
      evalsG(ip+1) = 0.5d0*Gvec**2
      ip = ip + 2
    ENDIF
  ENDDO
  !
  WRITE(*,*)
  WRITE(*,*) 'Eigenvalues based on Gvec:'
  DO ip=1,N
    WRITE(*,'(1x,I5,F18.10)') ip, evalsG(ip)
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
