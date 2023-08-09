program lxctest
  use xc_f03_lib_m
  implicit none
  TYPE(xc_f03_func_t) :: xc_func
  TYPE(xc_f03_func_info_t) :: xc_info
  integer(8), parameter :: npoints = 5
  real(8) :: rho(npoints) = (/0.1, 0.2, 0.3, 0.4, 0.5/)
  real(8) :: sigma(npoints) = (/0.2, 0.3, 0.4, 0.5, 0.6/)
  real(8) :: exc(npoints)
  integer :: i, vmajor, vminor, vmicro, func_id = 1
  ! Print out the version  
  call xc_f03_version(vmajor, vminor, vmicro)
  write(*,'("Libxc version: ",I1,".",I1,".",I1)') vmajor, vminor, vmicro
  ! Initialize the functional
  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
  ! Get the information
  xc_info = xc_f03_func_get_info(xc_func)
  ! Evaluate energy depending on the family
  select case (xc_f03_func_info_get_family(xc_info))
    case(XC_FAMILY_LDA)
      call xc_f03_lda_exc(xc_func, npoints, rho(1), exc(1))
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call xc_f03_gga_exc(xc_func, npoints, rho(1), sigma(1), exc(1))
  end select
  ! Print out density and energy density per particle
  do i = 1, npoints
    write(*,"(F8.6,1X,F9.6)") rho(i), exc(i)
  end do
  ! Deallocate memory
  call xc_f03_func_end(xc_func)
 end program lxctest
