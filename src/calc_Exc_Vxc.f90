SUBROUTINE calc_Exc_Vxc()
  USE m_hamiltonian, ONLY : Rhoe, V_xc
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_xc
  USE xc_f90_types_m
  USE xc_f90_lib_m
  IMPLICIT NONE 
  REAL(8), ALLOCATABLE :: gRhoe(:)
  REAL(8), ALLOCATABLE :: eps_x(:), eps_c(:)
  REAL(8), ALLOCATABLE :: vrho_x(:), vrho_c(:), vgrho_x(:), vgrho_c(:)
  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info

  ! these calls should only be done for LDA
  IF( XC_NAME == 'VWN' ) THEN 
    !
!    CALL excVWN( Npoints, Rhoe, EPS_XC )
!    CALL excpVWN( Npoints, Rhoe, d_EPS_XC_RHO )

    ALLOCATE( eps_x(Npoints), eps_c(Npoints) )
    ALLOCATE( vrho_x(Npoints), vrho_c(Npoints) )

    eps_x(:)  = 0.d0
    eps_c(:)  = 0.d0
    vrho_x(:) = 0.d0
    vrho_c(:) = 0.d0

    ! LDA exchange 
    CALL xc_f90_func_init(xc_func, xc_info, 1, XC_UNPOLARIZED)
    CALL xc_f90_lda_exc_vxc(xc_func, Npoints, Rhoe(1), eps_x(1), vrho_x(1))
    CALL xc_f90_func_end(xc_func)

    ! VWN correlation
    ! LDA_C_VWN_1 = 28
    ! LDA_C_VWN   = 7
    CALL xc_f90_func_init(xc_func, xc_info, 28, XC_UNPOLARIZED)
    CALL xc_f90_lda_exc_vxc(xc_func, Npoints, Rhoe(1), eps_c(1), vrho_c(1))
    CALL xc_f90_func_end(xc_func)

    EPS_XC(:) = eps_x(:) + eps_c(:)

    ! calculate potential
    V_xc(:) = EPS_XC(:) + vrho_x(:) + vrho_c(:)
    !
    DEALLOCATE( eps_x, eps_c )
    DEALLOCATE( vrho_x, vrho_c )

  ELSEIF( XC_NAME == 'PBE' ) THEN 
    !
    ALLOCATE( gRhoe(Npoints) )
    ALLOCATE( eps_x(Npoints), eps_c(Npoints) )
    ALLOCATE( vrho_x(Npoints), vrho_c(Npoints) )
    ALLOCATE( vgrho_x(Npoints), vgrho_c(Npoints) )

    ! PBE exchange 
    CALL xc_f90_func_init(xc_func, xc_info, 101, XC_UNPOLARIZED)
    CALL xc_f90_gga_exc_vxc(xc_func, Npoints, Rhoe(1), gRhoe(1), eps_x(1), vrho_x(1), vgrho_x(1))
    CALL xc_f90_func_end(xc_func)

    ! PBE correlation
    CALL xc_f90_func_init(xc_func, xc_info, 130, XC_UNPOLARIZED)
    CALL xc_f90_gga_exc_vxc(xc_func, Npoints, Rhoe(1), gRhoe(1), eps_c(1), vrho_c(1), vgrho_c(1))
    CALL xc_f90_func_end(xc_func)

    EPS_XC(:) = eps_x(:) + eps_c(:)
    d_EPS_XC_RHO(:) = vrho_x(:) + vrho_c(:)
    d_EPS_XC_GRHO(:) = vgrho_x(:) + vgrho_c(:)
    !
    ! calculate potential
    V_xc(:) = EPS_XC(:) + Rhoe(:)*d_EPS_XC_RHO(:) + gRhoe(:)*d_EPS_XC_GRHO(:)
    !
    DEALLOCATE( gRhoe )
    DEALLOCATE( eps_x, eps_c )
    DEALLOCATE( vrho_x, vrho_c )
    DEALLOCATE( vgrho_x, vgrho_c )
  ENDIF 

!  WRITE(*,*) 'sum(epsxc) = ', sum(epsxc)
!  WRITE(*,*) 'sum(V_xc) = ', sum(V_xc)

END SUBROUTINE 


