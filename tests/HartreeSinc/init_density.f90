SUBROUTINE init_density( num_gaussian, positions, coefs, exponents, density, potential )

  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : hh => LF3d_hh, &
                     Npoints => LF3d_Npoints, &
                     lingrid => LF3d_lingrid
  IMPLICIT NONE 
  !
  INTEGER :: num_gaussian
  REAL(8) :: positions(3,num_gaussian)
  REAL(8) :: exponents(num_gaussian), coefs(num_gaussian)
  REAL(8) :: density(Npoints)
  REAL(8) :: potential(Npoints)
  !
  REAL(8) :: weight, exponent2, r2
  REAL(8) :: r(3)
  INTEGER :: igauss, ip

  weight = hh(1)*hh(2)*hh(3)  ! FIXME use dVol ??

  density(:)   = 0.d0
  potential(:) = 0.d0
  !
  DO igauss = 1, num_gaussian
    exponent2 = exponents(igauss)**2
    DO ip = 1,Npoints
      r(:) = lingrid(:,ip) - positions(:,igauss)
      r2 = r(1)**2 + r(2)**2 + r(3)**2
      density(ip) = density(ip) + coefs(igauss)*(exponent2/pi)**1.5d0 * exp(-exponent2*r2)*sqrt(weight)
      ! the potential (analytic) due to the Gaussian density
      potential(ip) = potential(ip) + coefs(igauss)*erf(exponents(igauss)*sqrt(r2))/sqrt(r2)
    ENDDO 
  ENDDO 
END SUBROUTINE 

