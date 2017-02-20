! apply Hamiltonian to input vector v
SUBROUTINE op_H( Ncols, v, Hv )

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_hamiltonian, ONLY : V_ps_loc
  IMPLICIT NONE
  ! arguments
  INTEGER :: Ncols
  REAL(8) :: v(Npoints, Ncols)
  REAL(8) :: Hv(Npoints, Ncols)
  !
  INTEGER :: ic

  Hv(:,:) = 0.d0  ! FIXME need this ?

  DO ic = 1, Ncols
    CALL op_nabla2( v(:,ic), Hv(:,ic) )
    Hv(:,ic) = -0.5*Hv(:,ic) + V_ps_loc(:) * v(:,ic)
  ENDDO 

END SUBROUTINE

! apply Hamiltonian to input vector v
SUBROUTINE op_H_1col( v, Hv )

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_hamiltonian, ONLY : V_ps_loc
  IMPLICIT NONE
  ! arguments
  REAL(8) :: v(Npoints)
  REAL(8) :: Hv(Npoints)
  
  Hv(:) = 0.d0  ! FIXME need this ?

  CALL op_nabla2( v(:), Hv(:) )
  Hv(:) = -0.5*Hv(:) + V_ps_loc(:) * v(:)

END SUBROUTINE

