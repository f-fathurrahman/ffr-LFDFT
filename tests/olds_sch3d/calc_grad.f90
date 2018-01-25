! multicolumn version
! v are asssumed to be orthogonal
SUBROUTINE calc_grad( Ncols, v, grad )
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  IMPLICIT NONE
  !
  INTEGER :: Ncols
  REAL(8) :: v(Npoints,Ncols)
  REAL(8) :: grad(Npoints,Ncols)

  !
  REAL(8), ALLOCATABLE :: Hv(:)
  INTEGER :: ic, icc
  !
  REAL(8) :: ddot

  ALLOCATE( Hv(Npoints) )

  DO ic = 1, Ncols
    CALL op_H_1col( v(:,ic), Hv(:) )
    grad(:,ic) = Hv(:)
    DO icc = 1, Ncols
      grad(:,ic) = grad(:,ic) - ddot( Npoints, v(:,icc),1, Hv(:),1 )*v(:,icc)*dVol
    ENDDO
  ENDDO

  DEALLOCATE( Hv )
END SUBROUTINE
