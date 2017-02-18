! multicolumn version
! v are asssumed to be orthogonal
SUBROUTINE get_grad( Ncol, v, grad, do_ortho)
  USE m_globals, ONLY : N
  IMPLICIT NONE
  !
  INTEGER :: Ncol
  REAL(8) :: v(N**3,Ncol)
  REAL(8) :: grad(N**3,Ncol)
  LOGICAL, OPTIONAL :: do_ortho  ! FIXME The keyword OPTIONAL is not working (?)
  !
  REAL(8) :: Hv(N**3)
  INTEGER :: ic, icc
  !
  REAL(8) :: ddot

  IF( .NOT. present( do_ortho ) ) do_ortho = .FALSE.

  IF( do_ortho ) THEN
    CALL ortho_gram_schmidt( v, N**3, N**3, ncol )
  ENDIF

  !WRITE(*,*) 'Calling get_grad'
  !
  DO ic = 1, Ncol
    CALL apply_Ham( v(:,ic), Hv(:) )
    grad(:,ic) = Hv(:)
    DO icc = 1, Ncol
      grad(:,ic) = grad(:,ic) - ddot( N**3, v(:,icc),1, Hv(:),1 )*v(:,icc)
    ENDDO
  ENDDO
END SUBROUTINE
