
SUBROUTINE solve_poisson_cg( Nbasis, rho, phi )
  IMPLICIT NONE
  !
  INTEGER :: Nbasis
  REAL(8) :: rho(Nbasis), phi(Nbasis)
  !
  REAL(8), ALLOCATABLE :: r(:), p(:) ! residual
  REAL(8), ALLOCATABLE :: nabla2_phi(:)
  INTEGER :: ip, iter
  REAL(8) :: rsold, rsnew, alpha

  ALLOCATE( r(Nbasis), p(Nbasis), nabla2_phi(Nbasis) )

  !
  DO ip=1,Nbasis
    phi(ip) = 0.d0
  ENDDO

  CALL apply_laplacian( phi, nabla2_phi )
  r(:) = rho(:) - nabla2_phi(:)
  p(:) = r(:)

  !
  rsold = dot_product(r,r)
  !WRITE(*,*) 'rsold = ', rsold

  DO iter=1,1000
    CALL apply_laplacian( p, nabla2_phi )
    !
    alpha = rsold/dot_product(p,nabla2_phi)
    !
    phi(:) = phi(:) + alpha*p(:)
    !
    r(:) = r(:) - alpha*nabla2_phi(:)
    !
    rsnew = dot_product(r,r)
    !WRITE(*,*) 'cg-poisson = ', sqrt(rsnew)
    !
    IF(sqrt(rsnew) < 1.d-9) THEN
      WRITE(*,'(/1x,A,I5,E18.10)') 'Convergence in solve_poisson_cg in iter:', iter, sqrt(rsnew)
      EXIT
    ENDIF
    p(:) = r(:) + (rsnew/rsold)*p(:)
    rsold = rsnew
  ENDDO

  DEALLOCATE( r, p, nabla2_phi )
END SUBROUTINE


