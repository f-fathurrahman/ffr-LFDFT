PROGRAM lxctest

  USE xc_f90_types_m
  USE xc_f90_lib_m

  IMPLICIT NONE

  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  REAL(8) :: rho(5) = (/0.1, 0.2, 0.3, 0.4, 0.5/)
  REAL(8) :: sigma(5) = (/0.2, 0.3, 0.4, 0.5, 0.6/)
  REAL(8) :: exc(5)
  INTEGER :: i, vmajor, vminor, func_id = 1

  CALL xc_f90_version(vmajor, vminor)
  WRITE(*,'("Libxc version: ",I1,".",I1)') vmajor, vminor

  CALL xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)

  SELECT CASE (xc_f90_info_family(xc_info))
  CASE(XC_FAMILY_LDA)
    CALL xc_f90_lda_exc(xc_func, 5, rho(1), exc(1))
  CASE(XC_eFAMILY_GGA, XC_FAMILY_HYB_GGA)
    CALL xc_f90_gga_exc(xc_func, 5, rho(1), sigma(1), exc(1))
  END SELECT

  DO i = 1, 5
    WRITE(*,"(F8.6,1X,F9.6)") rho(i), exc(i)
  ENDDO

  CALL xc_f90_func_end(xc_func)

END PROGRAM lxctest

