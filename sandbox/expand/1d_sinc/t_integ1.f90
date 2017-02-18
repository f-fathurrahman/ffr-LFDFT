! Test the accuracy of integration using LF

FUNCTION gaussian( sigma, mu, x )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  REAL(8) :: gaussian
  REAL(8) ::sigma, mu, x

  gaussian = 1.d0/(sigma*sqrt(2.d0*PI))*exp( -0.5d0*( (x-mu)/sigma )**2 )
END FUNCTION


!------------------------------------------------------------------------------
PROGRAM t_integ
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
  INTEGER :: ii
  CHARACTER(20) :: buffer
  ! Functions
  REAL(8) :: gaussian
  !
  REAL(8) :: sigma, mu, mu2, h
  
  CALL getarg(1,buffer)
  READ(buffer,*) N
  WRITE(*,*) 'N = ', N

  CALL getarg(2,buffer)
  READ(buffer,*) h
  WRITE(*,*) 'h = ', h

  sigma = 0.3d0
  mu  = -2.d0
  mu2 = 2.d0

  ! Initialize the basis functions
  CALL init_LF1d_sinc( LF, N, h )
  A = LF%A
  B = LF%B
  L = LF%L
  
  ! Calculate coefficients of expansion for the function
  ALLOCATE( coefs(N) )

  DO ii=1,N
    coefs(ii) =  sqrt(h)*( gaussian( sigma, mu, LF%grid(ii) ) + &
       gaussian( sigma, mu2, LF%grid(ii) ) )
  ENDDO
  
  WRITE(*,*) 'N, integration test: ', N, sum( coefs(:)/sqrt(h) )*h
  WRITE(*,*) 'N, error: ', N, abs( 2.d0 - sum(coefs(:)/sqrt(h))*h )

  DEALLOCATE( coefs )

END PROGRAM

