PROGRAM lxctest
  USE xc_f90_types_m
  USE xc_f90_lib_m

  IMPLICIT NONE 

  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  !
  CHARACTER(120) :: str
  !
  REAL(8) :: rho(5) = (/0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0/)
  REAL(8) :: sigma(5) = (/0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0/)
  !
  REAL(8) :: vrho(5)
  REAL(8) :: vsigma(5)
  REAL(8) :: exc(5)
  !
  INTEGER :: i, vmajor, vminor, vmicro
  !INTEGER :: func_id = 1  ! LDA
  !INTEGER :: func_id = 101  ! X_PBE
  INTEGER :: func_id = 130  ! C_PBE

  CALL xc_f90_version(vmajor, vminor, vmicro)
  WRITE(*,'("Libxc version: ",I1,".",I1,".",I1)') vmajor, vminor, vmicro

  CALL xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)

  i = 0
  CALL xc_f90_info_refs(xc_info, i, str)
  DO WHILE( i >= 0 )
    WRITE(*,'(a,I1,2A)') '[',i,'] ', trim(str)
    i = i + 1
    CALL xc_f90_info_refs(xc_info, i, str)
  ENDDO 

  !
  ! evaluate exc, vrho and/or vsigma
  !
  SELECT CASE (xc_f90_info_family(xc_info))
  CASE(XC_FAMILY_LDA)
    WRITE(*,*) 'XC_FAMILY_LDA'
    !CALL xc_f90_lda_exc(xc_func, 5, rho(1), exc(1))
    CALL xc_f90_lda_exc_vxc(xc_func, 5, rho(1), exc(1), vrho(1))  
  CASE(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
    WRITE(*,*) 'XC_FAMILY_GGA'
    !CALL xc_f90_gga_exc(xc_func, 5, rho(1), sigma(1), exc(1))
    CALL xc_f90_gga_exc_vxc(xc_func, 5, rho(1), sigma(1), exc(1), vrho(1), vsigma(1))
  END SELECT 

  DO i = 1, 5
    WRITE(*,'(1x,5F18.10)') rho(i), sigma(i), exc(i), vrho(i), vsigma(i)
  ENDDO

  CALL xc_f90_func_end(xc_func)

END PROGRAM 



