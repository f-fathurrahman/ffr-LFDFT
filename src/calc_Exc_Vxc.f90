SUBROUTINE calc_Exc_Vxc()
  USE m_hamiltonian, ONLY : Rhoe, V_xc
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_xc
  IMPLICIT NONE 
  REAL(8), ALLOCATABLE :: gRhoe(:)

  ! these calls should only be done for LDA
  IF( XC_NAME == 'VWN' ) THEN 
    CALL excVWN( Npoints, Rhoe, EPS_XC )
    CALL excpVWN( Npoints, Rhoe, d_EPS_XC_RHO )
    ! calculate potential
    V_xc(:) = EPS_XC(:) + Rhoe(:)*d_EPS_XC_RHO(:)
  ELSEIF( XC_NAME == 'PBE' ) THEN 
    !
    ALLOCATE( gRhoe(Npoints) )
    !
    CALL excPBE( Npoints, Rhoe, gRhoe, EPS_XC )
    CALL excpPBE( Npoints, Rhoe, gRhoe, d_EPS_XC_RHO, d_EPS_XC_GRHO )
    !
    ! calculate potential
    V_xc(:) = EPS_XC(:) + Rhoe(:)*d_EPS_XC_RHO(:) + gRhoe(:)*d_EPS_XC_GRHO(:)
    !
    DEALLOCATE( gRhoe )
  ENDIF 

!  WRITE(*,*) 'sum(epsxc) = ', sum(epsxc)
!  WRITE(*,*) 'sum(V_xc) = ', sum(V_xc)

END SUBROUTINE 


