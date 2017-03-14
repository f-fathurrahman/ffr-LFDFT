SUBROUTINE init_V_ps_loc_H_hgh( Npoints, r, V )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  INTEGER :: Npoints
  REAL(8) :: r(Npoints)
  REAL(8) :: V(Npoints)
  INTEGER :: ip
  REAL(8) :: r1, r2
  !
  REAL(8) :: z_val
  REAL(8) :: rlocal
  REAL(8) :: c(2)

  z_val = 1.d0
  rlocal = 0.2d0
  c(1) = -4.180237d0
  c(2) = 0.725075d0

  DO ip = 1, Npoints

    r1 = r(ip)/rlocal
    r2 = r1**2

    IF(r(ip) < 1.0d-7) THEN
      V(ip) = - (2.d0 * z_val)/(sqrt(2.d0*PI)*rlocal) + c(1)
      WRITE(*,*) 'Small r, using limiting value'
    ELSE
      ! using erf from intrinsic function
      V(ip) = - (z_val/r(ip))*erf(r1/sqrt(2.d0)) + exp( -0.5d0*r2 ) * ( c(1) + c(2)*r2 )
    ENDIF

  ENDDO


END SUBROUTINE

