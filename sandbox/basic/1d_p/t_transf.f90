PROGRAM t_transf
  USE m_LF1d
  IMPLICIT NONE
  TYPE(LF1d_t) :: LF
  INTEGER, PARAMETER :: N = 9
  REAL(8), PARAMETER :: L = 2.d0
  INTEGER :: i, j, ip
  COMPLEX(8), ALLOCATABLE :: invTal(:,:)
  COMPLEX(8), ALLOCATABLE :: TT(:,:)
  COMPLEX(8), ALLOCATABLE :: KinT(:,:)
  REAL(8), ALLOCATABLE :: vecIn(:), vecOut(:)
  COMPLEX(8), ALLOCATABLE :: vecTmp(:)
  REAL(8) :: Gvec
  REAL(8), ALLOCATABLE :: evecsT(:,:), evalsT(:)
  INTEGER, ALLOCATABLE :: idxG(:)

  ALLOCATE( invTal(N,N), TT(N,N), KinT(N,N), vecTmp(N) )

  ALLOCATE( vecIn(N), vecOut(N), evecsT(N,N), evalsT(N) )

  ALLOCATE( idxG(N) )

  idxG(1) = 0
  ip = 1
  DO i=2,N,2
    !WRITE(*,*) i
    idxG(i) = ip
    idxG(i+1) = ip
    ip = ip + 1
  ENDDO

  DO i=1,N
    WRITE(*,*) idxG(i)
  ENDDO

  CALL init_LF1d_p( LF, N, -L/2.d0, L/2.d0 )
  CALL info_LF1d( LF, .TRUE. )

  evecsT = -0.5d0*LF%D2jl
  CALL eig_dsyev( evecsT, evalsT, N )

!  WRITE(*,*)
!  WRITE(*,*) 'Matrix Tal'
!  DO i = 1, N
!    WRITE(*,*)
!    DO j = 1, N
!      WRITE(*,'(A,2F10.5,A)', advance='no') '(',LF%Tal(i,j),')'
!    ENDDO
!  ENDDO
!  WRITE(*,*)

  invTal(:,:) = LF%Tal(:,:)
  !CALL c8_inverse(N,invTal)  ! not working?
 
!  WRITE(*,*)
!  WRITE(*,*) 'Inverse of matrix Tal'
!  DO i = 1, N
!    WRITE(*,*)
!    DO j = 1, N
!      WRITE(*,'(A,2F10.5,A)', advance='no') '(',invTal(i,j),')'
!    ENDDO
!  ENDDO
!  WRITE(*,*)

  CALL calc_invTal( N, LF, invTal )

  TT = matmul( invTal, LF%Tal )

!  WRITE(*,*)
!  WRITE(*,*) 'invTal * Tal'
!  DO i = 1, N
!    WRITE(*,*)
!    DO j = 1, N
!      WRITE(*,'(A,2F10.5,A)', advance='no') '(',TT(i,j),')'
!    ENDDO
!  ENDDO
!  WRITE(*,*)

  KinT = matmul( matmul( invTal, -0.5 * LF%D2jl ), LF%Tal )

!  WRITE(*,*)
!  WRITE(*,*) 'KinT'
!  DO i = 1, N
!    WRITE(*,*)
!    DO j = 1, N
!      WRITE(*,'(A,2F10.5,A)', advance='no') '(',KinT(i,j),')'
!    ENDDO
!  ENDDO
!  WRITE(*,*)


  vecIn(:) = 1.d0
  !vecIn(1) = 5.d0
  vecTmp = matmul( invTal, vecIn )



  WRITE(*,*) 'Before'
  DO i=1,N
    WRITE(*,'(3F18.10)') vecIn(i), vecTmp(i)
  ENDDO
!  ip = 1
!  DO i=-N/2,N,2
!   !WRITE(*,*) vecTmp(i)
!    IF( i /= 0 ) THEN
!      Gvec = 2.d0*PI/L * i 
!      vecTmp(ip) = vecTmp(ip)/( 0.5*Gvec**2 + 1 )
!    ENDIF
!    ip = ip + 1
!  ENDDO
  WRITE(*,*) 'After'
  DO i=1,N
    WRITE(*,'(2F18.10)') vecTmp(i)
  ENDDO
  
!  DO i=1,N
!    IF( idxG(i) /= 0 ) THEN
!      Gvec = 2.d0*PI/L * idxG(i)
!      vecTmp(i) = vecTmp(ip)/(0.5d0*Gvec**2 )
!    ENDIF
!  ENDDO

  vecOut = real( matmul( LF%Tal, vecTmp ) )

  DO i = 1, N 
    WRITE(*,'(I5,2F18.10)',advance='no') i, evalsT(i), vecOut(i)
    IF( evalsT(i) >= 1d-10 ) THEN
      WRITE(*,'(F18.10)', advance='no') vecIn(i)/( evalsT(i) )
    ENDIF
    WRITE(*,*)
  ENDDO


  DEALLOCATE( invTal, TT, KinT, vecTmp )
  DEALLOCATE( vecIn, vecOut, evecsT, evalsT )

CONTAINS


! FIXME Not working?
SUBROUTINE calc_invTal( N, LF, invTal )
  IMPLICIT NONE
  INTEGER :: N
  TYPE(LF1d_t) :: LF
  COMPLEX(8) :: invTal(N,N)
  INTEGER :: kl, alpha, l

  kl = -(LF%N-1)/2
  DO l = 1, LF%N
    DO alpha = 1, LF%N
      !WRITE(*,*) 'kl = ', kl
      invTal(alpha,l) = exp( -2.d0*PI*(0.d0,1.d0)*kl*LF%grid(alpha)/LF%L ) / sqrt(dble(LF%N))
    ENDDO
    kl = kl + 1
  ENDDO

END SUBROUTINE


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

END PROGRAM



