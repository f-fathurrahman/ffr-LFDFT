SUBROUTINE diag_Ham( Nbasis, Nstates, evecs, evals )
  IMPLICIT NONE
  ! arguments
  INTEGER :: Nstates
  INTEGER :: Nbasis
  !
  REAL(8) :: evecs(Nbasis,Nstates)
  REAL(8) :: evals(Nstates)
  ! Local
  INTEGER :: gstart
  REAL(8) :: ethr
  INTEGER, ALLOCATABLE :: btype(:)
  LOGICAL :: lrot
  INTEGER :: dav_iter, notcnv
  INTEGER :: ist
  
  gstart = 2
  ethr = 1.d-5
  ALLOCATE( btype(Nstates) )
  btype(:) = 1 ! all bands are occupied
  lrot = .FALSE.
 
  ! random initial guess of eigenvectors
  DO ist = 1, Nstates
    CALL r8_rand_vec( Nbasis, evecs(:,ist) )
  ENDDO
  CALL ortho_gram_schmidt( evecs, Nbasis, Nbasis, Nstates )

  !WRITE(*,*) 'Calling regterg ...'
  CALL regterg( Nbasis, Nbasis, Nstates, 4*Nstates, evecs, ethr, &
                    gstart, evals, btype, notcnv, lrot, dav_iter )
  !WRITE(*,*) 'dav_iter = ', dav_iter
  !WRITE(*,*) 'notcnv   = ', notcnv

  ! renormalize eigenstates
  evecs(:,:) = sqrt(2.d0)*evecs(:,:)
  !DO ist = 1, Nstates
  !  WRITE(*,*) evals(ist)
  !ENDDO

  DEALLOCATE( btype )


END SUBROUTINE

