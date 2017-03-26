
!
! solve A*x = b, where A = H (Hamiltonian matrix)
!
! NOTE: x will be used as initial guess.
! Explicitly set x to zeros before calling this subroutine if needed.
SUBROUTINE linsolve_H( b, x )
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  IMPLICIT NONE
  !
  REAL(8) :: b(Npoints), x(Npoints)
  !
  REAL(8), ALLOCATABLE :: r(:), p(:) ! residual
  REAL(8), ALLOCATABLE :: Hx(:)
  INTEGER :: iter
  REAL(8) :: rsold, rsnew, alpha

  ALLOCATE( r(Npoints), p(Npoints), Hx(Npoints) )

  CALL op_H( x, Hx )
  r(:) = b(:) - Hx(:)
  p(:) = r(:)

  !
  rsold = dot_product(r,r)
  !WRITE(*,*) 'rsold = ', rsold

  DO iter=1,1000
    CALL op_H( p, Hx )
    !
    alpha = rsold/dot_product(p,Hx)
    !
    x(:) = x(:) + alpha*p(:)
    !
    r(:) = r(:) - alpha*Hx(:)
    !
    rsnew = dot_product(r,r)
    !WRITE(*,*) 'cg-poisson = ', sqrt(rsnew)
    !
    IF(sqrt(rsnew) < 1.d-10) THEN
    !IF(rsnew < 1.d-10) THEN
      WRITE(*,*) 'linsolve_H converged in iter:', iter
      EXIT
    ENDIF
    p(:) = r(:) + (rsnew/rsold)*p(:)
    rsold = rsnew
  ENDDO

  DEALLOCATE( r, p, Hx )
END SUBROUTINE


