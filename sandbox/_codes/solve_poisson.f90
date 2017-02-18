
SUBROUTINE solve_poisson_cg( rho, phi )
  USE m_globals, ONLY : N
  IMPLICIT NONE
  !
  REAL(8) :: rho(N**3), phi(N**3)
  !
  REAL(8), ALLOCATABLE :: r(:), p(:) ! residual
  REAL(8), ALLOCATABLE :: nabla2_phi(:)
  INTEGER :: ip, iter
  REAL(8) :: rsold, rsnew, alpha

  ALLOCATE( r(N**3), p(N**3), nabla2_phi(N**3) )

  !
  DO ip=1,N
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
    !WRITE(*,*) 'rsnew = ', sqrt(rsnew)
    !
    IF(sqrt(rsnew) < 1.d-10) THEN
    !IF(rsnew < 1.d-10) THEN
      WRITE(*,*) 'Convergence in solve_poisson_cg in iter:', iter
      EXIT
    ENDIF
    p(:) = r(:) + (rsnew/rsold)*p(:)
    rsold = rsnew
  ENDDO

  DEALLOCATE( r, p, nabla2_phi )
END SUBROUTINE


