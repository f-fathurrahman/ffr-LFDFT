!
! solve A*x = b, where A = H (Hamiltonian matrix)
!
! NOTE: x will be used as initial guess.
! Explicitly set x to zeros before calling this subroutine if needed.
SUBROUTINE prec_linsolve_cg_H( b, x, NiterMax )

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  IMPLICIT NONE
  !
  REAL(8) :: b(Npoints), x(Npoints)
  INTEGER :: NiterMax
  !
  REAL(8), ALLOCATABLE :: r(:), p(:) ! residual
  REAL(8), ALLOCATABLE :: Hx(:)
  INTEGER :: iter
  REAL(8) :: rsold, rsnew, alpha
  REAL(8) :: ddot

  ALLOCATE( r(Npoints), p(Npoints), Hx(Npoints) )

  CALL op_H( 1, x, Hx )
  r(:) = b(:) - Hx(:)
  p(:) = r(:)

  !
  rsold = ddot( Npoints, r,1, r,1 )

  DO iter=1,NiterMax

    CALL op_H( 1, p, Hx )
    !
    alpha = rsold/ddot( Npoints, p,1, Hx,1 )
    !
    x(:) = x(:) + alpha*p(:)
    !
    r(:) = r(:) - alpha*Hx(:)
    !
    rsnew = ddot( Npoints, r,1, r,1 )
    !
    WRITE(*,*) 'iterCG, sqrt(rsnew) = ', iter, sqrt(rsnew)
    IF(sqrt(rsnew) < 5.d-9) THEN
      WRITE(*,*) 'prec_linsolve_cg_H converged in iter:', iter
      GOTO 100
    ENDIF
    p(:) = r(:) + (rsnew/rsold)*p(:)
    rsold = rsnew
  ENDDO

  WRITE(*,*) 'prec_linsolve_cg_H is not converged after iter:', iter
  WRITE(*,*) 'Last sqrt(rsnew) = ', sqrt(rsnew)

  100 DEALLOCATE( r, p, Hx )
END SUBROUTINE


