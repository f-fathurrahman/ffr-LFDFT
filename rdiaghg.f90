!----------------------------------------------------------------------------
SUBROUTINE rdiaghg( n, m, h, s, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates eigenvalues and eigenvectors of the generalized problem
  ! ... Hv=eSv, with H symmetric matrix, S overlap matrix.
  ! ... On output both matrix are unchanged
  !
  ! ... LAPACK version - uses both DSYGV and DSYGVX
  !
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: DP=8
  INTEGER, INTENT(IN) :: n, m, ldh
    ! dimension of the matrix to be diagonalized
    ! number of eigenstates to be calculated
    ! leading dimension of h, as declared in the calling pgm unit
  REAL(DP), INTENT(INOUT) :: h(ldh,n), s(ldh,n)
    ! matrix to be diagonalized
    ! overlap matrix
  !
  REAL(DP), INTENT(OUT) :: e(n)
    ! eigenvalues
  REAL(DP), INTENT(OUT) :: v(ldh,m)
    ! eigenvectors (column-wise)
  !
  INTEGER               :: i, j, lwork, nb, mm, info
    ! mm = number of calculated eigenvectors
  REAL(DP)              :: abstol
  REAL(DP), PARAMETER   :: one = 1_DP
  REAL(DP), PARAMETER   :: zero = 0_DP
  INTEGER,  ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP), ALLOCATABLE :: work(:), sdiag(:), hdiag(:)
  LOGICAL               :: all_eigenvalues
  INTEGER,  EXTERNAL    :: ILAENV
  ! ILAENV returns optimal block size "nb"

  !
  ! ... save the diagonal of input S (it will be overwritten)
  !
  ALLOCATE( sdiag( n ) )
  DO i = 1, n
     sdiag(i) = s(i,i)
  END DO
  !
  all_eigenvalues = ( m == n )
  !
  ! ... check for optimal block size
  !
  nb = ILAENV( 1, 'DSYTRD', 'U', n, -1, -1, -1 )
  !
  IF ( nb < 5 .OR. nb >= n ) THEN
     !
     lwork = 8*n
     !
  ELSE
     !
     lwork = ( nb + 3 )*n
     !
  END IF
  !
  ALLOCATE( work( lwork ) )
  !
  IF ( all_eigenvalues ) THEN
     !
     ! ... calculate all eigenvalues
     !
     v(:,:) = h(:,:)
     !
     CALL DSYGV( 1, 'V', 'U', n, v, ldh, s, ldh, e, work, lwork, info )
     !
  ELSE
     !
     ! ... calculate only m lowest eigenvalues
     !
     ALLOCATE( iwork( 5*n ) )
     ALLOCATE( ifail( n ) )
     !
     ! ... save the diagonal of input H (it will be overwritten)
     !
     ALLOCATE( hdiag( n ) )
     DO i = 1, n
        hdiag(i) = h(i,i)
     END DO
     !
     abstol = 0.D0
     ! abstol = 2.D0*DLAMCH( 'S' )
     !
     CALL DSYGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
                  0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
                  work, lwork, iwork, ifail, info )
     !
     DEALLOCATE( ifail )
     DEALLOCATE( iwork )
     !
     ! ... restore input H matrix from saved diagonal and lower triangle
     !
     DO i = 1, n
        h(i,i) = hdiag(i)
        DO j = i + 1, n
           h(i,j) = h(j,i)
        END DO
        DO j = n + 1, ldh
           h(j,i) = 0.0_DP
        END DO
     END DO
     !
     DEALLOCATE( hdiag )
     !
  END IF
  !
  DEALLOCATE( work )
  !
  IF ( info > n ) THEN
    WRITE(*,*) 'ERROR in rdiaghg: S matrix not positive definite', ABS( info )
    STOP
  ENDIF
  IF ( info > 0 ) THEN 
    WRITE(*,*) 'ERROR in rdiaghg: eigenvectors failed to converge', ABS( info )
    STOP 
  ENDIF 
  IF ( info < 0 ) THEN
    WRITE(*,*) 'ERROR in rdiaghg: incorrect call to DSYGV*', ABS( info )
    STOP 
  ENDIF 
  
  ! ... restore input S matrix from saved diagonal and lower triangle
  !
  DO i = 1, n
     s(i,i) = sdiag(i)
     DO j = i + 1, n
        s(i,j) = s(j,i)
     END DO
     DO j = n + 1, ldh
        s(j,i) = 0.0_DP
     END DO
  END DO
  !
  DEALLOCATE( sdiag )
  !
  RETURN
  !
END SUBROUTINE

