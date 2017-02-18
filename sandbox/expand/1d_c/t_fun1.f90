FUNCTION gaussian( sigma, mu, x )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  REAL(8) :: gaussian
  REAL(8) ::sigma, mu, x

  gaussian = 1.d0/(sigma*sqrt(2.d0*PI))*exp( -0.5d0*( (x-mu)/sigma )**2 )
END FUNCTION


!------------------------------------------------------------------------------
PROGRAM t_fun
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF1d
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: A, B, L
  !
  TYPE(LF1d_t) :: LF
  REAL(8), ALLOCATABLE :: coefs(:)
  !
  INTEGER :: ii, jj
  ! Functions
  REAL(8) :: gaussian
  !
  REAL(8) :: sigma, mu, mu2, h
  INTEGER :: NPTS_PLOT=201
  REAL(8) :: xx, yy

  A = -5.d0
  B = 5.d0
  L = B - A
  N = 30

  ! Gaussian parameters
  sigma = 0.3d0
  mu    = -2.d0
  mu2   = 2.d0

  ! Initialize the basis functions
  CALL init_LF1d_c( LF, N, A, B )
  CALL info_LF1d( LF, .TRUE. )
  h = LF%h
  
  ! Calculate coefficients of expansion for the function
  ALLOCATE( coefs(N) )

  DO ii=1,N
    coefs(ii) =  sqrt(h)*( gaussian( sigma, mu, LF%grid(ii) ) + &
       gaussian( sigma, mu2, LF%grid(ii) ) )
    WRITE(11,*) LF%grid(ii), coefs(ii)/sqrt(h)
  ENDDO
  
  WRITE(*,*) 'integration test: ', sum( coefs(:)/sqrt(h) )*h

  ! Write to file for plotting
  DO ii=1,NPTS_PLOT
    !
    xx = A + (ii-1)*L/NPTS_PLOT
    yy = 0.d0
    DO jj=1,N
      yy = yy + eval_LF1d_c( LF, jj, xx )*coefs(jj)
    ENDDO
    WRITE(12,*) xx
    WRITE(13,*) yy
    WRITE(14,*) gaussian( sigma, mu, xx ) + gaussian( sigma, mu2, xx )
  ENDDO

  ! Free memory
  CALL dealloc_LF1d( LF )
  DEALLOCATE( coefs )


END PROGRAM

