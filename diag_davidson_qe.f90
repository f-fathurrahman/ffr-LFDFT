
! gstart = 2
! lrot = .false.
SUBROUTINE diag_davidson_qe( Nbasis, nvec, nvecx, evc, ethr, &
                             e, btype, notcnv, dav_iter )
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: DP=8
  INTEGER, PARAMETER :: stdout=6
  REAL(DP), PARAMETER :: ZERO=0.d0
  INTEGER, INTENT(IN) :: Nbasis, nvec, nvecx
    ! dimension of the matrix to be diagonalized
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set
    !    (the basis set is refreshed when its dimension would exceed nvecx)
  REAL(DP), INTENT(INOUT) :: evc(Nbasis,nvec)
    !  evc   contains the  refined estimates of the eigenvectors
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence: root improvement is stopped,
    ! when two consecutive estimates of the root differ by less than ethr.
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  REAL(DP), INTENT(OUT) :: e(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer  number of iterations performed
    ! number of unconverged roots
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 100
    ! maximum number of iterations
  !
  INTEGER :: kter, Nred, np, n, m, nb1
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
    ! counter on the bands
  REAL(DP), ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:), ew(:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! eigenvectors of the Hamiltonian
    ! eigenvalues of the reduced hamiltonian
  REAL(DP), ALLOCATABLE :: psi(:,:), hpsi(:,:)
    ! work space, contains psi
    ! the product of H and psi
  LOGICAL, ALLOCATABLE :: conv(:) ! true if the root is converged
  REAL(DP) :: empty_ethr ! threshold for empty bands
  REAL(DP), EXTERNAL :: ddot
  INTEGER :: ib
  !
  IF ( nvec > nvecx / 2 ) THEN
    WRITE(*,*) 'Error in diag_davidson_qe: nvec is too small'
    STOP
  ENDIF
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  ALLOCATE( psi(  Nbasis, nvecx ) )
  ALLOCATE( hpsi( Nbasis, nvecx ) )
  !
  ALLOCATE( sr( nvecx, nvecx ) )
  ALLOCATE( hr( nvecx, nvecx ) )
  ALLOCATE( vr( nvecx, nvecx ) )
  ALLOCATE( ew( nvecx ) )
  ALLOCATE( conv( nvec ) )
  !
  notcnv = nvec
  Nred  = nvec
  conv   = .FALSE.
  !
  hpsi = ZERO
  psi  = ZERO
  psi(:,1:nvec) = evc(:,1:nvec)
  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL op_H( nvec, psi, hpsi )
  !
  !
  ! ... hr contains the projection of the hamiltonian onto the reduced
  ! ... space vr contains the eigenvectors of hr
  !
  hr(:,:) = 0.D0
  sr(:,:) = 0.D0
  vr(:,:) = 0.D0
  ! hr = psi^T * hpsi
  CALL DGEMM( 'T', 'N', Nred, Nred, Nbasis, 1.D0 , &  ! 2.0 -> 1.0
              psi, Nbasis, hpsi, Nbasis, 0.D0, hr, nvecx )
  !
  CALL DGER( Nred, Nred, -1.D0, psi, Nbasis, hpsi, Nbasis, hr, nvecx )
  !
  CALL DGEMM( 'T', 'N', Nred, Nred, Nbasis, 1.D0, &  ! 2.0 -> 1.0
              psi, Nbasis, psi, Nbasis, 0.D0, sr, nvecx )
  !
  CALL DGER( Nred, Nred, -1.D0, psi, Nbasis, psi, Nbasis, sr, nvecx )
  !

  !
  ! ... diagonalize the reduced hamiltonian
  !
  CALL rdiaghg( Nred, nvec, hr, sr, nvecx, ew, vr )
  !
  e(1:nvec) = ew(1:nvec)
  !
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     WRITE(*,*) 'kter = ', kter
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
           ew(Nred+np) = e(n)
           !
        END IF
        !
     END DO
     !
     nb1 = Nred + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL DGEMM( 'N', 'N', Nbasis, notcnv, Nred, 1.D0, &
                  psi, Nbasis, vr, nvecx, 0.D0, psi(1,nb1), Nbasis )
     !
     DO np = 1, notcnv
        !
        psi(:,Nred+np) = - ew(Nred+np) * psi(:,Nred+np)
        !
     END DO
     !
     CALL DGEMM( 'N', 'N', Nbasis, notcnv, Nred, 1.D0, &
                 hpsi, Nbasis, vr, nvecx, 1.D0, psi(1,nb1), Nbasis )
     !
     !
     ! ... approximate inverse iteration
     !
     !CALL g_psi( Nbasis, Nbasis, notcnv, 1, psi(1,nb1), ew(nb1) )
     !WRITE(*,*) 'nb1 = ', nb1
     DO ib = 1,notcnv
       CALL prec_ilu0_inplace( psi(:,nb1+ib-1) )
     ENDDO
     !
     !
     ! ... "normalize" correction vectors psi(:,nb1:Nred+notcnv) in
     ! ... order to improve numerical stability of subspace diagonalization
     ! ... (rdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = Nred + 1, Nred + notcnv
     !
     DO n = 1, notcnv
        !
        ew(n) = 2.D0 * ddot( Nbasis, psi(1,Nred+n), 1, psi(1,Nred+n), 1 )
        ew(n) = ew(n) - psi(1,Nred+n) * psi(1,Nred+n)
        !
     END DO
     !
     DO n = 1, notcnv
        !
        psi(:,Nred+n) = psi(:,Nred+n) / SQRT( ew(n) )
        !
     END DO
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     CALL op_H( notcnv, psi(1,nb1), hpsi(1,nb1) )
     !
     ! ... update the reduced hamiltonian
     !
     CALL DGEMM( 'T', 'N', Nred+notcnv, notcnv, Nbasis, 1.D0, psi, &  ! 2.0 -> 1.0
                 Nbasis, hpsi(1,nb1), Nbasis, 0.D0, hr(1,nb1), nvecx )
     CALL DGER( Nred+notcnv, notcnv, -1.D0, psi, &
                 Nbasis, hpsi(1,nb1), Nbasis, hr(1,nb1), nvecx )
     !
     !
     CALL DGEMM( 'T', 'N', Nred+notcnv, notcnv, Nbasis, 1.D0, psi, &  ! 2 -> 1.0
                 Nbasis, psi(1,nb1), Nbasis, 0.D0, sr(1,nb1) , nvecx )
     !
     CALL DGER( Nred+notcnv, notcnv, -1.D0, psi, &
                 Nbasis, psi(1,nb1), Nbasis, sr(1,nb1), nvecx )

     !
     Nred = Nred + notcnv
     !
     DO n = 1, Nred
        !
        DO m = n + 1, Nred
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
     CALL rdiaghg( Nred, nvec, hr, sr, nvecx, ew, vr )
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
     !DO n = 1,nvec
     !  WRITE(*,*) n, ABS( ew(1:nvec) - e(1:nvec) )
     !ENDDO
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
          Nred+notcnv > nvecx .OR. dav_iter == maxter ) THEN
        !        !
        CALL DGEMM( 'N', 'N', Nbasis, nvec, Nred, 1.D0, &
                    psi, Nbasis, vr, nvecx, 0.D0, evc, Nbasis )
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           !WRITE(*,*) 'All root converged ...'
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
        CALL DGEMM( 'N', 'N', Nbasis, nvec, Nred, 1.D0, hpsi, &
                    Nbasis, vr, nvecx, 0.D0, psi(1,nvec+1), Nbasis )
        !
        hpsi(:,1:nvec) = psi(:,nvec+1:nvec+nvec)
        !
        ! ... refresh the reduced hamiltonian
        !
        Nred = nvec
        !
        hr(:,1:Nred) = 0.D0
        sr(:,1:Nred) = 0.D0
        vr(:,1:Nred) = 0.D0
        !
        DO n = 1, Nred
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
END SUBROUTINE
