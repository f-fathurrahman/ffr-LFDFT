
FUNCTION deltaMN(m,n) RESULT( ff )
  IMPLICIT NONE
  INTEGER :: m, n
  REAL(8) :: ff
  IF( m /= n ) THEN
    ff = 0.d0
  ELSE
    ff = 1.d0
  ENDIF
END FUNCTION



!----------------------------------
SUBROUTINE eig_dsyev(A,eigval,dimA)
!----------------------------------
  IMPLICIT NONE
  ! Arguments
  INTEGER :: dimA
  REAL(8) :: A(dimA,dimA)
  REAL(8) :: eigval(dimA)
  ! Workspace array for DSYEV
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


! real(8) version
!------------------------------------------------
SUBROUTINE ortho_gram_schmidt(v, ldv, nrow, ncol)
!------------------------------------------------
  IMPLICIT NONE
  INTEGER :: ldv, nrow, ncol
  REAL(8) :: v(ldv,ncol)
  !
  INTEGER :: ii, jj
  REAL(8) :: zz, puv
  !
  REAL(8) :: ddot

  DO ii = 1, ncol
    zz = ddot( nrow, v(1:nrow,ii),1, v(1:nrow,ii),1 )
    v(1:nrow,ii) = v(1:nrow,ii)/sqrt( zz )
    DO jj = ii+1, ncol
      puv = prj( nrow, v(1:nrow,ii), v(1:nrow,jj) )
      v(1:nrow,jj) = v(1:nrow,jj) - puv*v(1:nrow,ii)
    ENDDO
  ENDDO

  CONTAINS

    ! compute prj = <v|u>/<u|u>
    FUNCTION prj(N,u,v)
      IMPLICIT NONE
      !
      REAL(8) :: prj
      INTEGER :: N
      REAL(8) :: u(N), v(N)
      !
      REAL(8) :: vu, uu
      REAL(8) :: ddot
      !
      ! FIXME: I got the vectors to be orthogonal when I reverse the arguments
      ! for zdotc
      vu = ddot( N, u,1, v,1 )
      uu = ddot( N, u,1, u,1 )
      prj = vu/uu
    END FUNCTION

END SUBROUTINE


!---------------------------
SUBROUTINE r8_rand_vec(N, v)
!---------------------------
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: v(N)
  !
  INTEGER :: i

  DO i=1,N
    CALL random_number( v(i) )
  ENDDO

END SUBROUTINE
