SUBROUTINE excVWN_libxc( Npts, rho, epsxc )
  USE m_constants, ONLY : PI
  USE xc_f90_types_m
  USE xc_f90_lib_m
  IMPLICIT NONE
  INTEGER :: Npts
  REAL(8) :: rho(Npts)
  REAL(8) :: epsxc(Npts)
  !
  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  INTEGER :: func_id

  func_id = 1

  CALL xc_f90_func_init( xc_func, xc_info, func_id, XC_UNPOLARIZED )

  CALL xc_f90_lda_exc( xc_func, Npts, rho(1), epsxc(1) )

  CALL xc_f90_func_end( xc_func )
END



SUBROUTINE excpVWN_libxc( Npts, rho, depsxc )
  USE m_constants, ONLY : PI
  USE xc_f90_types_m
  USE xc_f90_lib_m
  IMPLICIT NONE
  INTEGER :: Npts
  REAL(8) :: rho(Npts)
  REAL(8) :: depsxc(Npts)
  !
  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  INTEGER :: func_id

  func_id = 1

  CALL xc_f90_func_init( xc_func, xc_info, func_id, XC_UNPOLARIZED )

  CALL xc_f90_lda_vxc( xc_func, Npts, rho(1), depsxc(1) )

  CALL xc_f90_func_end( xc_func )

END SUBROUTINE

