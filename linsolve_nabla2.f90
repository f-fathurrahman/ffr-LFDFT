
!
! solve A*x = b, where A = nabla2 (Laplacian matrix
!
SUBROUTINE linsolve_nabla2( b, x )
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  IMPLICIT NONE
  !
  REAL(8) :: b(Npoints), x(Npoints)
  !
  REAL(8), ALLOCATABLE :: r(:), p(:) ! residual
  REAL(8), ALLOCATABLE :: nabla2_x(:)
  INTEGER :: ip, iter
  REAL(8) :: rsold, rsnew, alpha

  ALLOCATE( r(Npoints), p(Npoints), nabla2_x(Npoints) )

  ! TODO: Add option to use x as input
  DO ip=1,Npoints
    x(ip) = 0.d0
  ENDDO

  CALL op_nabla2( x, nabla2_x )
  r(:) = b(:) - nabla2_x(:)
  p(:) = r(:)

  !
  rsold = dot_product(r,r)
  !WRITE(*,*) 'rsold = ', rsold

  DO iter=1,1000
    CALL apply_laplacian( p, nabla2_x )
    !
    alpha = rsold/dot_product(p,nabla2_x)
    !
    x(:) = x(:) + alpha*p(:)
    !
    r(:) = r(:) - alpha*nabla2_x(:)
    !
    rsnew = dot_product(r,r)
    !WRITE(*,*) 'cg-poisson = ', sqrt(rsnew)
    !
    IF(sqrt(rsnew) < 1.d-10) THEN
    !IF(rsnew < 1.d-10) THEN
      WRITE(*,*) 'Convergence in solve_poisson_cg in iter:', iter
      EXIT
    ENDIF
    p(:) = r(:) + (rsnew/rsold)*p(:)
    rsold = rsnew
  ENDDO

  DEALLOCATE( r, p, nabla2_x )
END SUBROUTINE


