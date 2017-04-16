
!
! solve A*x = b, where A = H (Hamiltonian matrix)
!
! NOTE: x will be used as initial guess.
! Explicitly set x to zeros before calling this subroutine if needed.
SUBROUTINE prec_linsolve_pcg_H( b, x, NiterMax )

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  IMPLICIT NONE
  !
  REAL(8) :: b(Npoints), x(Npoints)
  INTEGER :: NiterMax
  !
  REAL(8), ALLOCATABLE :: r(:), p(:), z(:), r_old(:)
  REAL(8), ALLOCATABLE :: Hx(:)
  INTEGER :: iter
  REAL(8) :: rsold, rsnew, alpha
  REAL(8) :: ddot

  ALLOCATE( r(Npoints), p(Npoints), Hx(Npoints), z(Npoints), r_old(Npoints) )

  CALL op_H( 1, x, Hx )
  r(:) = b(:) - Hx(:)
  CALL prec_H_diag( r, z )
  p(:) = z(:)

  rsold = ddot( Npoints, r,1, z,1 )
  rsnew = 0.d0
  r_old(:) = r(:)

  DO iter=1,NiterMax

    CALL op_H( 1, p, Hx )
    !
    alpha = rsold/ddot( Npoints, p,1, Hx,1 )
    !
    x(:) = x(:) + alpha*p(:)
    !
    r(:) = r(:) - alpha*Hx(:)
    CALL prec_H_diag( r, z )
    rsnew = ddot( Npoints, z,1, r,1 )
    !
    WRITE(*,*) 'iterCG, sqrt(rsnew) = ', iter, sqrt(abs(rsnew))
    IF(sqrt(abs(rsnew)) < 1.d-9) THEN
      WRITE(*,*) 'prec_linsolve_pcg_H converged in iter:', iter
      GOTO 100
    ENDIF
    !p(:) = z(:) + (rsnew/rsold)*p(:)
    p(:) = z(:) + ddot(Npoints, z,1, r(:) - r_old(:), 1)/rsold*p(:)
    rsold = rsnew
    r_old(:) = r(:)
  ENDDO

  WRITE(*,*) 'prec_linsolve_cg_H is not converged after iter:', iter
  WRITE(*,*) 'Last sqrt(rsnew) = ', sqrt(rsnew)

  100 DEALLOCATE( r, p, Hx, z, r_old )
END SUBROUTINE


