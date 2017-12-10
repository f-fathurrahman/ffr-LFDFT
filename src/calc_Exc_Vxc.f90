SUBROUTINE calc_Exc_Vxc()
  USE m_hamiltonian, ONLY : Rhoe, V_xc
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_xc
  USE xc_f90_types_m
  USE xc_f90_lib_m
  IMPLICIT NONE 
  REAL(8), ALLOCATABLE :: gRhoe(:,:), gRhoe2(:)
  REAL(8), ALLOCATABLE :: eps_x(:), eps_c(:)
  REAL(8), ALLOCATABLE :: vrho_x(:), vrho_c(:), vgrho_x(:), vgrho_c(:)
  REAL(8), ALLOCATABLE :: h(:,:), divh(:)
  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  INTEGER :: ip

  ! these calls should only be done for LDA
  IF( XC_NAME == 'VWN' ) THEN 
    !
    !CALL excVWN( Npoints, Rhoe, EPS_XC )
    !CALL excpVWN( Npoints, Rhoe, d_EPS_XC_RHO )

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
    ! FIXME: need to optimize memory uses here
    !
    ALLOCATE( gRhoe(3,Npoints) )
    ALLOCATE( gRhoe2(Npoints) )
    ALLOCATE( eps_x(Npoints), eps_c(Npoints) )
    ALLOCATE( vrho_x(Npoints), vrho_c(Npoints) )
    ALLOCATE( vgrho_x(Npoints), vgrho_c(Npoints) )
    ALLOCATE( h(3,Npoints) )
    ALLOCATE( divh(Npoints) )

    CALL op_nabla( Rhoe, gRhoe )

    DO ip = 1,Npoints
      gRhoe2(ip) = gRhoe(1,ip)**2 + gRhoe(2,ip)**2 + gRhoe(3,ip)**2
    ENDDO 

    ! PBE exchange 
    CALL xc_f90_func_init(xc_func, xc_info, 101, XC_UNPOLARIZED)
    CALL xc_f90_gga_exc_vxc(xc_func, Npoints, Rhoe(1), gRhoe2(1), eps_x(1), vrho_x(1), vgrho_x(1))
    CALL xc_f90_func_end(xc_func)

    ! PBE correlation
    CALL xc_f90_func_init(xc_func, xc_info, 130, XC_UNPOLARIZED)
    CALL xc_f90_gga_exc_vxc(xc_func, Npoints, Rhoe(1), gRhoe2(1), eps_c(1), vrho_c(1), vgrho_c(1))
    CALL xc_f90_func_end(xc_func)

    EPS_XC(:) = eps_x(:) + eps_c(:)
    
    ! vgrho * gRhoe
    DO ip = 1,Npoints
      h(:,ip) = ( vgrho_x(ip) + vgrho_c(ip) ) * gRhoe(:,ip)
    ENDDO 

    ! div ( vgrho * gRhoe )
    CALL op_nabla_dot( h, divh )

    ! need factor 2 for divh ???
    DO ip = 1,Npoints 
!      V_xc(ip) = EPS_XC(ip) + vrho_x(ip) + vrho_c(ip) - 2.d0*divh(ip)
      V_xc(ip) = vrho_x(ip) + vrho_c(ip) - 2d0*divh(ip)
    ENDDO 

    !
    DEALLOCATE( gRhoe )
    DEALLOCATE( gRhoe2 )
    DEALLOCATE( eps_x, eps_c )
    DEALLOCATE( vrho_x, vrho_c )
    DEALLOCATE( vgrho_x, vgrho_c )
    DEALLOCATE( h )
    DEALLOCATE( divh )
  ENDIF 

END SUBROUTINE 


