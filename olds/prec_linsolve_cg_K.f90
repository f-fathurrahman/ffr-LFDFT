
!
! solve A*x = b, where A = K (Kinetic operator matrix)
!
! NOTE: x will be used as initial guess.
! Explicitly set x to zeros before calling this subroutine if needed.
SUBROUTINE prec_linsolve_cg_K( b, x, NiterMax )

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  IMPLICIT NONE
  !
  REAL(8) :: b(Npoints), x(Npoints)
  INTEGER :: NiterMax
  !
  REAL(8), ALLOCATABLE :: r(:), p(:) ! residual
  REAL(8), ALLOCATABLE :: Kx(:)
  INTEGER :: iter
  REAL(8) :: rsold, rsnew, alpha
  REAL(8) :: ddot

  ALLOCATE( r(Npoints), p(Npoints), Kx(Npoints) )

  CALL op_K( x, Kx )
  r(:) = b(:) - Kx(:)
  p(:) = r(:)

  !
  rsold = ddot( Npoints, r,1, r,1 )

  DO iter=1,NiterMax

    CALL op_K( p, Kx )
    !
    alpha = rsold/ddot( Npoints, p,1, Kx,1 )
    !
    x(:) = x(:) + alpha*p(:)
    !
    r(:) = r(:) - alpha*Kx(:)
    !
    rsnew = ddot( Npoints, r,1, r,1 )
    !
    !WRITE(*,*) 'iterCG, sqrt(rsnew) = ', iter, sqrt(rsnew)
    IF(sqrt(rsnew) < 1.d-10) THEN
      WRITE(*,*) 'prec_linsolve_cg_K converged in iter:', iter
      GOTO 100
    ENDIF
    p(:) = r(:) + (rsnew/rsold)*p(:)
    rsold = rsnew
  ENDDO

  WRITE(*,*) 'prec_linsolve_cg_K is not converged after iter:', iter
  WRITE(*,*) 'Last sqrt(rsnew) = ', sqrt(rsnew)

  100 DEALLOCATE( r, p, Kx )
END SUBROUTINE


