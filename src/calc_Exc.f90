SUBROUTINE calc_Exc()
  USE m_energies, ONLY : E_xc
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, dVol => LF3d_dVol
  USE m_hamiltonian, ONLY : Rhoe 
  USE m_xc
  USE xc_f90_types_m
  USE xc_f90_lib_m
  USE m_options, ONLY : USE_ARIAS_VWN
  IMPLICIT NONE 
  REAL(8), ALLOCATABLE :: gRhoe(:,:), gRhoe2(:)
  REAL(8), ALLOCATABLE :: eps_x(:), eps_c(:)
  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  INTEGER :: ip

  ! these calls should only be done for LDA
  IF( XC_NAME == 'VWN' ) THEN 

    IF( USE_ARIAS_VWN ) THEN 
      CALL excVWN( Npoints, Rhoe, EPS_XC )
    ELSE

    ALLOCATE( eps_x(Npoints), eps_c(Npoints) )
    eps_x(1:Npoints) = 0.d0
    eps_c(1:Npoints) = 0.d0

    ! LDA exchange 
    CALL xc_f90_func_init(xc_func, xc_info, 1, XC_UNPOLARIZED)
    CALL xc_f90_lda_exc(xc_func, Npoints, Rhoe(1), eps_x(1))
    CALL xc_f90_func_end(xc_func)

    ! VWN correlation
    ! LDA_C_VWN_1 = 28
    ! LDA_C_VWN   = 7
    CALL xc_f90_func_init(xc_func, xc_info, 28, XC_UNPOLARIZED)
    CALL xc_f90_lda_exc(xc_func, Npoints, Rhoe(1), eps_c(1))
    CALL xc_f90_func_end(xc_func)

    EPS_XC(:) = eps_x(:) + eps_c(:)
    ENDIF ! USE_ARIAS_VWN

    E_xc = sum( Rhoe(:) * EPS_XC(:) )*dVol

  ELSEIF( XC_NAME == 'PBE' ) THEN 
    !
    ALLOCATE( gRhoe(3,Npoints) )
    ALLOCATE( gRhoe2(Npoints) )
    ALLOCATE( eps_x(Npoints), eps_c(Npoints) )

    CALL op_nabla( Rhoe, gRhoe )
    DO ip = 1,Npoints
      gRhoe2(ip) = gRhoe(1,ip)**2 + gRhoe(2,ip)**2 + gRhoe(3,ip)**2
    ENDDO 

    ! PBE exchange 
    CALL xc_f90_func_init(xc_func, xc_info, 101, XC_UNPOLARIZED)
    CALL xc_f90_gga_exc(xc_func, Npoints, Rhoe(1), gRhoe2(1), eps_x(1))
    CALL xc_f90_func_end(xc_func)

    ! PBE correlation
    CALL xc_f90_func_init(xc_func, xc_info, 130, XC_UNPOLARIZED)
    CALL xc_f90_gga_exc(xc_func, Npoints, Rhoe(1), gRhoe2(1), eps_c(1))
    CALL xc_f90_func_end(xc_func)

    EPS_XC(:) = eps_x(:) + eps_c(:)

    E_xc = sum( Rhoe(:) * EPS_XC(:) )*dVol

!    WRITE(*,*) 'E_xc PBE = ', E_xc

    !
    DEALLOCATE( gRhoe )
    DEALLOCATE( eps_x, eps_c )

  ENDIF 

END SUBROUTINE 

