PROGRAM compare
  IMPLICIT NONE 
  INTEGER :: Npoints
  REAL(8), ALLOCATABLE :: Rhoe(:)

  Npoints = 5
  ALLOCATE( Rhoe(Npoints) )
  Rhoe(:) = (/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0 /)

  CALL test1()
  CALL test2()

CONTAINS 


SUBROUTINE test1()
  REAL(8), ALLOCATABLE :: EPS_XC(:), d_EPS_XC_RHO(:), V_xc(:)
  INTEGER :: i

  ALLOCATE( V_xc(Npoints) )
  ALLOCATE( EPS_XC(Npoints) )
  ALLOCATE( d_EPS_XC_RHO(Npoints) )

  CALL excVWN( Npoints, Rhoe, EPS_XC )
  CALL excpVWN( Npoints, Rhoe, d_EPS_XC_RHO )

  ! calculate potential
  V_xc(:) = EPS_XC(:) + Rhoe(:)*d_EPS_XC_RHO(:)

  DO i = 1, Npoints
    WRITE(*,'(1x,4F18.10)') Rhoe(i), EPS_XC(i), d_EPS_XC_RHO(i), V_xc(i)
  ENDDO 

END SUBROUTINE 


SUBROUTINE test2()

  USE xc_f90_types_m
  USE xc_f90_lib_m

  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  REAL(8), ALLOCATABLE :: EPS_XC(:), d_EPS_XC_RHO(:)
  INTEGER :: i
  REAL(8), ALLOCATABLE :: eps_x(:), eps_c(:)
  REAL(8), ALLOCATABLE :: vrho_x(:), vrho_c(:)

  ALLOCATE( EPS_XC(Npoints) )
  ALLOCATE( d_EPS_XC_RHO(Npoints) )

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
  CALL xc_f90_func_init(xc_func, xc_info, 7, XC_UNPOLARIZED)
  CALL xc_f90_lda_exc_vxc(xc_func, Npoints, Rhoe(1), eps_c(1), vrho_c(1))
  CALL xc_f90_func_end(xc_func)

  EPS_XC(:) = eps_x(:) + eps_c(:)
  d_EPS_XC_RHO(:) = vrho_x(:) + vrho_c(:)

  WRITE(*,*)
  DO i = 1, Npoints
    WRITE(*,'(1x,2F18.10,18x,F18.10)') Rhoe(i), EPS_XC(i), d_EPS_XC_RHO(i)
  ENDDO 

  !
  DEALLOCATE( eps_x, eps_c )
  DEALLOCATE( vrho_x, vrho_c )


END SUBROUTINE 

END PROGRAM 

