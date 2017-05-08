! apply Hamiltonian to input vector v
SUBROUTINE op_H( Ncols, v, Hv )

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_hamiltonian, ONLY : V_ps_loc, V_Hartree, V_xc
  USE m_nabla2_sparse, ONLY : nzval => nabla2_nzval, &
                              rowval => nabla2_rowval, &
                              colptr => nabla2_colptr
  USE m_options, ONLY : FREE_NABLA2
  USE m_PsPot, ONLY : NbetaNL
  IMPLICIT NONE
  ! arguments
  INTEGER :: Ncols
  REAL(8) :: v(Npoints, Ncols)
  REAL(8) :: Hv(Npoints, Ncols)
  REAL(8), ALLOCATABLE :: V_ps_NL_psi(:,:)
  !
  INTEGER :: ic

  Hv(:,:) = 0.d0  ! FIXME need this ?

  IF( FREE_NABLA2 ) THEN 
    DO ic = 1, Ncols
      CALL op_nabla2( v(:,ic), Hv(:,ic) )
      Hv(:,ic) = -0.5*Hv(:,ic) + ( V_ps_loc(:) + V_Hartree(:) + V_xc(:) )* v(:,ic)
    ENDDO 
  ELSE 
    DO ic = 1, Ncols
      CALL amux( Npoints, v(:,ic), Hv(:,ic), nzval, rowval, colptr )
      Hv(:,ic) = -0.5*Hv(:,ic) + ( V_ps_loc(:) + V_Hartree(:) + V_xc(:) )* v(:,ic)
    ENDDO 
  ENDIF 

  ! add non-local contrib if any
  IF( NbetaNL > 0 ) THEN 
    ALLOCATE( V_ps_NL_psi(Npoints,Ncols) )
    CALL op_V_ps_NL( Ncols, V_ps_NL_psi(:,:) )
    Hv(:,:) = Hv(:,:) + V_ps_NL_psi(:,:)
    DEALLOCATE( V_ps_NL_psi )
  ENDIF 

END SUBROUTINE



! apply Hamiltonian to input vector v
SUBROUTINE op_H_1col( ist, v, Hv )

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_hamiltonian, ONLY : V_ps_loc, V_Hartree, V_xc
  USE m_nabla2_sparse, ONLY : nzval => nabla2_nzval, &
                              rowval => nabla2_rowval, &
                              colptr => nabla2_colptr
  USE m_options, ONLY : FREE_NABLA2
  USE m_PsPot, ONLY : NbetaNL
  IMPLICIT NONE
  ! arguments
  INTEGER :: ist
  REAL(8) :: v(Npoints)
  REAL(8) :: Hv(Npoints)
  REAL(8), ALLOCATABLE :: V_ps_NL_psi(:)
  
  Hv(:) = 0.d0  ! FIXME need this ?

  IF( FREE_NABLA2 ) THEN 
    CALL op_nabla2( v(:), Hv(:) )
    Hv(:) = -0.5*Hv(:) + ( V_ps_loc(:) + V_Hartree(:) + V_xc(:) ) * v(:)
  ELSE 
    CALL amux( Npoints, v(:), Hv(:), nzval, rowval, colptr )
    Hv(:) = -0.5*Hv(:) + ( V_ps_loc(:) + V_Hartree(:) + V_xc(:) )* v(:)
  ENDIF 

  ! add non-local contrib if any
  IF( NbetaNL > 0 ) THEN 
    ALLOCATE( V_ps_NL_psi(Npoints) )

    CALL op_V_ps_NL_1col( ist, V_ps_NL_psi(:) )
    Hv(:) = Hv(:) + V_ps_NL_psi(:)

    DEALLOCATE( V_ps_NL_psi )
  ENDIF 

END SUBROUTINE

