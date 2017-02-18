! from pwscf of QuantumEspresso distribution
! modified by ffr

!----------------------------------------------------------------------------
SUBROUTINE regterg( npw, npwx, nvec, nvecx, evc, ethr, &
                    gstart, e, btype, notcnv, lrot, dav_iter )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an uspp matrix, evc is a complex vector
  ! ... (real wavefunctions with only half plane waves stored)
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: DP=8
  INTEGER, PARAMETER :: stdout=6
  REAL(DP), PARAMETER :: ZERO=0.d0
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx, gstart
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set
    !    (the basis set is refreshed when its dimension would exceed nvecx)
  REAL(DP), INTENT(INOUT) :: evc(npwx,nvec)
    !  evc   contains the  refined estimates of the eigenvectors
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence: root improvement is stopped,
    ! when two consecutive estimates of the root differ by less than ethr.
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  LOGICAL, INTENT(IN) :: lrot
    ! .TRUE. if the wfc have already been rotated
  REAL(DP), INTENT(OUT) :: e(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer  number of iterations performed
    ! number of unconverged roots
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 500
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, n, m, nb1, ibnd
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
    ! counter on the bands
  INTEGER :: ierr
  REAL(DP), ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:), ew(:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! eigenvectors of the Hamiltonian
    ! eigenvalues of the reduced hamiltonian
  REAL(DP), ALLOCATABLE :: psi(:,:), hpsi(:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  LOGICAL, ALLOCATABLE :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr 
    ! threshold for empty bands
  INTEGER :: npw2, npwx2
  !
  REAL(DP), EXTERNAL :: ddot
  !
  INTEGER :: i,j
  !
  ! EXTERNAL  h_psi, s_psi, g_psi
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi> 
    ! s_psi(npwx,npw,nvec,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,nvec)
    ! g_psi(npwx,npw,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  IF ( nvec > nvecx / 2 ) THEN
    WRITE(*,*) 'Error in regterg: nvec is too small'
    STOP
  ENDIF
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  ALLOCATE( psi(  npwx, nvecx ), STAT=ierr )
  ALLOCATE( hpsi( npwx, nvecx ), STAT=ierr )
  !
  ALLOCATE( sr( nvecx, nvecx ), STAT=ierr )
  ALLOCATE( hr( nvecx, nvecx ), STAT=ierr )
  ALLOCATE( vr( nvecx, nvecx ), STAT=ierr )
  ALLOCATE( ew( nvecx ), STAT=ierr )
  ALLOCATE( conv( nvec ), STAT=ierr )
  !
  npw2  = npw
  npwx2  = npwx
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
  hpsi = ZERO
  psi  = ZERO
  psi(:,1:nvec) = evc(:,1:nvec)
  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL h_psi( npwx, npw, nvec, psi, hpsi )
  !
  !
  ! ... hr contains the projection of the hamiltonian onto the reduced
  ! ... space vr contains the eigenvectors of hr
  !
  hr(:,:) = 0.D0
  sr(:,:) = 0.D0
  vr(:,:) = 0.D0
  !
  CALL DGEMM( 'T', 'N', nbase, nbase, npw, 2.D0 , &
              psi, npwx, hpsi, npwx, 0.D0, hr, nvecx )
  !
  IF ( gstart == 2 ) CALL DGER( nbase, nbase, -1.D0, psi, npwx, hpsi, npwx, hr, nvecx ) 
  !
  CALL DGEMM( 'T', 'N', nbase, nbase, npw, 2.D0, &
              psi, npwx, psi, npwx, 0.D0, sr, nvecx )
  !
  IF ( gstart == 2 ) CALL DGER( nbase, nbase, -1.D0, psi, npwx, psi, npwx, sr, nvecx )
  !
  IF ( lrot ) THEN
     !
     DO n = 1, nbase
        !
        e(n) = hr(n,n)
        vr(n,n) = 1.D0
        !
     END DO
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL rdiaghg( nbase, nvec, hr, sr, nvecx, ew, vr )
     !
     e(1:nvec) = ew(1:nvec)
     !
  END IF
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
     dav_iter = kter
     !     !
     np = 0
     !
     DO n = 1, nvec
        !
        IF ( .NOT. conv(n) ) THEN
           !
           ! ... this root not yet converged ... 
           !
           np = np + 1
           !
           ! ... reorder eigenvectors so that coefficients for unconverged
           ! ... roots come first. This allows to use quick matrix-matrix 
           ! ... multiplications to set a new basis vector (see below)
           !
           IF ( np /= n ) vr(:,np) = vr(:,n)
           !
           ! ... for use in g_psi
           !
           ew(nbase+np) = e(n)
           !   
        END IF
        !
     END DO
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL DGEMM( 'N', 'N', npw2, notcnv, nbase, 1.D0, &
                  psi, npwx2, vr, nvecx, 0.D0, psi(1,nb1), npwx2 )
     !
     DO np = 1, notcnv
        !
        psi(:,nbase+np) = - ew(nbase+np) * psi(:,nbase+np)
        !
     END DO
     !
     CALL DGEMM( 'N', 'N', npw2, notcnv, nbase, 1.D0, &
                 hpsi, npwx2, vr, nvecx, 1.D0, psi(1,nb1), npwx2 )
     !
     !
     ! ... approximate inverse iteration
     !
     !CALL g_psi( npwx, npw, notcnv, psi(1,nb1), ew(nb1) )
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in 
     ! ... order to improve numerical stability of subspace diagonalization 
     ! ... (rdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     DO n = 1, notcnv
        !
        ew(n) = 2.D0 * ddot( npw2, psi(1,nbase+n), 1, psi(1,nbase+n), 1 )
        IF ( gstart == 2 ) ew(n) = ew(n) - psi(1,nbase+n) * psi(1,nbase+n)
        !
     END DO
     !
     DO n = 1, notcnv
        !
        psi(:,nbase+n) = psi(:,nbase+n) / SQRT( ew(n) )
        !
     END DO
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     CALL h_psi( npwx, npw, notcnv, psi(1,nb1), hpsi(1,nb1) )
     !
     ! ... update the reduced hamiltonian
     !
     CALL DGEMM( 'T', 'N', nbase+notcnv, notcnv, npw2, 2.D0, psi, &
                 npwx2, hpsi(1,nb1), npwx2, 0.D0, hr(1,nb1), nvecx )
     IF ( gstart == 2 ) CALL DGER( nbase+notcnv, notcnv, -1.D0, psi, &
                 npwx2, hpsi(1,nb1), npwx2, hr(1,nb1), nvecx )
     !
     !
     CALL DGEMM( 'T', 'N', nbase+notcnv, notcnv, npw2, 2.D0, psi, &
                 npwx2, psi(1,nb1), npwx2, 0.D0, sr(1,nb1) , nvecx )
     !
     IF ( gstart == 2 ) CALL DGER( nbase+notcnv, notcnv, -1.D0, psi, &
                 npwx2, psi(1,nb1), npwx2, sr(1,nb1), nvecx )
    
     !
     nbase = nbase + notcnv
     !
     DO n = 1, nbase
        !
        DO m = n + 1, nbase
           !
           hr(m,n) = hr(n,m)
           sr(m,n) = sr(n,m)
           !
        END DO
        !
     END DO
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL rdiaghg( nbase, nvec, hr, sr, nvecx, ew, vr )
     !
     ! ... test for convergence
     !
     WHERE( btype(1:nvec) == 1 )
        !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < ethr ) )
        !
     ELSEWHERE
        !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < empty_ethr ) )
        !
     END WHERE
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
     e(1:nvec) = ew(1:nvec)
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     IF ( notcnv == 0 .OR. &
          nbase+notcnv > nvecx .OR. dav_iter == maxter ) THEN
        !        !
        CALL DGEMM( 'N', 'N', npw2, nvec, nbase, 1.D0, &
                    psi, npwx2, vr, nvecx, 0.D0, evc, npwx2 )
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           WRITE(*,*) 'All root converged ...'
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           WRITE( stdout, '(5X,"WARNING: ",I5, &
                &   " eigenvalues not converged in regterg")' ) notcnv
           !           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        psi(:,1:nvec) = evc(:,1:nvec)
        !
        !
        CALL DGEMM( 'N', 'N', npw2, nvec, nbase, 1.D0, hpsi, &
                    npwx2, vr, nvecx, 0.D0, psi(1,nvec+1), npwx2 )
        !
        hpsi(:,1:nvec) = psi(:,nvec+1:nvec+nvec)
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        hr(:,1:nbase) = 0.D0
        sr(:,1:nbase) = 0.D0
        vr(:,1:nbase) = 0.D0
        !
        DO n = 1, nbase
           !
           hr(n,n) = e(n)
           sr(n,n) = 1.D0
           vr(n,n) = 1.D0
           !
        END DO
        !        !
     END IF
     !
  END DO iterate
  !
  DEALLOCATE( conv )
  DEALLOCATE( ew )
  DEALLOCATE( vr )
  DEALLOCATE( hr )
  DEALLOCATE( sr )
  !
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )  
  !  !
  RETURN
  !
END SUBROUTINE regterg



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
#if defined (__ESSL)
     !
     ! ... there is a name conflict between essl and lapack ...
     !
     CALL DSYGV( 1, v, ldh, s, ldh, e, v, ldh, n, work, lwork )
     !
     info = 0
#else
     CALL DSYGV( 1, 'V', 'U', n, v, ldh, s, ldh, e, work, lwork, info )
#endif
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
  IF ( info > n ) WRITE(*,*) 'rdiaghg: S matrix not positive definite', ABS( info )
  IF ( info > 0 ) WRITE(*,*) 'rdiaghg: eigenvectors failed to converge', ABS( info )
  IF ( info < 0 ) WRITE(*,*) 'rdiaghg: incorrect call to DSYGV*', ABS( info )
  
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
END SUBROUTINE rdiaghg



