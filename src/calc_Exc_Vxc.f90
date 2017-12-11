SUBROUTINE calc_Exc_Vxc()
  USE m_hamiltonian, ONLY : Rhoe, V_xc
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_xc
  USE xc_f90_types_m
  USE xc_f90_lib_m
  IMPLICIT NONE 
  REAL(8), ALLOCATABLE :: gRhoe(:,:), gRhoe2(:)
  REAL(8), ALLOCATABLE :: vrho_x(:), vrho_c(:), vgrho_x(:), vgrho_c(:)
  REAL(8), ALLOCATABLE :: h(:,:), divh(:)
  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  INTEGER :: ip
  REAL(8), ALLOCATABLE :: d_EPS_XC_RHO(:), EPS_XC(:)

  ! these calls should only be done for LDA
  IF( XC_NAME == 'VWN' ) THEN 
    
    IF( USE_ARIAS_VWN ) THEN 
      ALLOCATE( d_EPS_XC_RHO(Npoints) )
      ALLOCATE( EPS_XC(Npoints) )
      CALL excVWN( Npoints, Rhoe, EPS_XC )
      CALL excpVWN( Npoints, Rhoe, d_EPS_XC_RHO )
      ! calculate potential
      V_xc(:) = EPS_XC(:) + d_EPS_XC_RHO(:)*Rhoe(:)
      DEALLOCATE( d_EPS_XC_RHO )
      DEALLOCATE( EPS_XC )
    ELSE 

    ALLOCATE( vrho_x(Npoints), vrho_c(Npoints) )

    vrho_x(1:Npoints) = 0.d0
    vrho_c(1:Npoints) = 0.d0

    ! LDA exchange 
    CALL xc_f90_func_init(xc_func, xc_info, 1, XC_UNPOLARIZED)
    CALL xc_f90_lda_vxc(xc_func, Npoints, Rhoe(1), vrho_x(1))
    CALL xc_f90_func_end(xc_func)

    ! VWN correlation
    ! LDA_C_VWN_1 = 28
    ! LDA_C_VWN   = 7
    CALL xc_f90_func_init(xc_func, xc_info, 28, XC_UNPOLARIZED)
    CALL xc_f90_lda_vxc(xc_func, Npoints, Rhoe(1), vrho_c(1))
    CALL xc_f90_func_end(xc_func)

    ! calculate potential
    DO ip = 1,Npoints
!      V_xc(ip) = eps_x(ip) + eps_c(ip) + (vrho_x(ip) + vrho_c(ip))*Rhoe(ip)
      V_xc(ip) = vrho_x(ip) + vrho_c(ip)
    ENDDO 
    
    DEALLOCATE( vrho_x, vrho_c )
    
    ENDIF ! USE_ARIAS_VWN

  ELSEIF( XC_NAME == 'PBE' ) THEN 
    !
    ! FIXME: need to optimize memory uses here
    !
    ALLOCATE( gRhoe(3,Npoints) )
    ALLOCATE( gRhoe2(Npoints) )
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
    CALL xc_f90_gga_vxc(xc_func, Npoints, Rhoe(1), gRhoe2(1), vrho_x(1), vgrho_x(1))
    CALL xc_f90_func_end(xc_func)

    ! PBE correlation
    CALL xc_f90_func_init(xc_func, xc_info, 130, XC_UNPOLARIZED)
    CALL xc_f90_gga_vxc(xc_func, Npoints, Rhoe(1), gRhoe2(1), vrho_c(1), vgrho_c(1))
    CALL xc_f90_func_end(xc_func)

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
    DEALLOCATE( vrho_x, vrho_c )
    DEALLOCATE( vgrho_x, vgrho_c )
    DEALLOCATE( h )
    DEALLOCATE( divh )
  ENDIF 

END SUBROUTINE 


