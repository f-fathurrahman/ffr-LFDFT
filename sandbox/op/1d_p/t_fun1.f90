FUNCTION gaussian( sigma, mu, x )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  REAL(8) :: gaussian
  REAL(8) ::sigma, mu, x

  gaussian = 1.d0/(sigma*sqrt(2.d0*PI))*exp( -0.5d0*( (x-mu)/sigma )**2 )
END FUNCTION

!
FUNCTION dgaussian( sigma, mu, x )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  REAL(8) :: dgaussian
  REAL(8) :: sigma, mu, x

  dgaussian = -sqrt(2.d0)*(-2.d0*mu + 2.d0*x)*exp(-(-mu + x)**2/sigma**2)/(2*sqrt(pi)*sigma**3)
END FUNCTION

FUNCTION d2gaussian( sigma, mu,x )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  REAL(8) :: d2gaussian
  REAL(8) :: sigma, mu, x

 d2gaussian = -sqrt(2.d0)*exp(-(-mu + x)**2/sigma**2)/(sqrt(PI)*sigma**3) + &
   sqrt(2.d0)*(-2.d0*mu + 2.d0*x)**2*exp(-(-mu + x)**2/sigma**2)/(2*sqrt(PI)*sigma**5)
END FUNCTION


!------------------------------------------------------------------------------
PROGRAM t_fun
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF1d
  IMPLICIT NONE
  !
  INTEGER :: N = 201
  REAL(8) :: L = 2.d0*PI
  !
  TYPE(LF1d_t) :: LF
  REAL(8), ALLOCATABLE :: coefs(:), dcoefs(:), d2coefs(:)
  !
  INTEGER :: ii, jj
  ! Functions
  REAL(8) :: gaussian, dgaussian, d2gaussian
  !
  REAL(8) :: sigma, mu, mu2
  INTEGER :: NPTS_PLOT=201
  REAL(8) :: xx, yy

  sigma = 0.12d0
  mu  = L/3.d0
  mu2 = 2.d0*L/3.d0

  ! Initialize the basis functions
  CALL init_LF1d_p( LF, N, 0.d0, L )

  ! Calculate coefficients of expansion for the function
  ALLOCATE( coefs(N) )
  ALLOCATE( dcoefs(N) )
  ALLOCATE( d2coefs(N) )

  DO ii=1,N
    coefs(ii) = sqrt(L/N)*( gaussian( sigma, mu, LF%grid(ii) ) + &
       gaussian( sigma, mu2, LF%grid(ii) ) )
    WRITE(11,*) LF%grid(ii), sqrt(N/L)*coefs(ii)
  ENDDO

  dcoefs = matmul( LF%D1jl, coefs )
  d2coefs = matmul( LF%D2jl, coefs )

  ! Write to file for plotting
  DO ii=1,NPTS_PLOT
    !
    !xx = (ii-1)*L/NPTS_PLOT
    xx = LF%grid(ii)  ! only works for NPTS_PLOT = N
    yy = 0.d0
    DO jj=1,N
      yy = yy + eval_LF1d_p( LF, jj, xx )*coefs(jj)
    ENDDO
    WRITE(12,*) xx
    WRITE(13,*) yy
    WRITE(14,*) gaussian( sigma, mu, xx ) + gaussian( sigma, mu2, xx )
    !
    ! First derivative
    !
    yy = 0.d0
    DO jj=1,N
      !yy = yy + eval_LF1d_p_d1( LF, jj, xx )*coefs(jj)
      yy = yy + eval_LF1d_p( LF, jj, xx )*dcoefs(jj)
    ENDDO
    WRITE(15,*) yy
    WRITE(16,*) dgaussian( sigma, mu, xx ) + dgaussian( sigma, mu2, xx )
    !WRITE(*,*) yy/dgaussian( sigma, mu, xx )
    !
    ! Second derivative
    !
    yy = 0.d0
    DO jj=1,N
      yy = yy + eval_LF1d_p( LF, jj, xx )*d2coefs(jj)
    ENDDO
    WRITE(17,*) yy
    WRITE(18,*) d2gaussian( sigma, mu, xx ) + d2gaussian( sigma, mu2, xx )
    !WRITE(*,*) yy/d2gaussian( sigma, mu, xx )
  ENDDO

  ! Free memory
  CALL dealloc_LF1d( LF )
  DEALLOCATE( coefs, dcoefs, d2coefs )


END PROGRAM

