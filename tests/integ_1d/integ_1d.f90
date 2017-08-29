! test convergence of integrals
! TODO: also test egg-box effect in 1d ???

! A periodic function
FUNCTION funcx( c1, L, x )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  REAL(8) :: funcx
  REAL(8) :: c1, L, x, dx

  dx = abs(x - c1)
  funcx = 1.d-6*exp( 25.d0*cos(2.d0*PI/L*dx) )
END FUNCTION


! non-periodic function
FUNCTION funcx2( A, alpha, c1, x )
  IMPLICIT NONE 
  REAL(8) :: c1, x, funcx2, A, alpha
  funcx2 = A*exp( -alpha*(x-c1)**2 )
END FUNCTION 


! periodic function
FUNCTION funcx3( A, alpha, c1, L, x )
  IMPLICIT NONE 
  REAL(8) :: c1, x, funcx3, A, alpha, dx1, dx2, dx3, dx, L
  dx1 = abs( x - c1 )
  dx2 = abs( L + x - c1 )
  dx3 = abs( x - c1 - L )
  IF( dx1 < dx2 ) THEN 
    dx = dx1
  ELSE 
    dx = dx2
  ENDIF 
  IF( dx3 < dx ) THEN 
    dx = dx3
  ENDIF 
  funcx3 = A*exp( -alpha*dx**2 )
END FUNCTION 


PROGRAM test_integral
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : grid_x => LF3d_grid_x
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: L
  !
  REAL(8), ALLOCATABLE :: coefs(:)
  !
  INTEGER :: ii, jj
  ! Functions
  REAL(8) :: funcx, funcx2, funcx3
  !
  REAL(8) :: A, alpha
  REAL(8) :: xx, yy, h, c1
  !
  REAL(8) :: eval_LF1d_p
  !
  CHARACTER(8) :: chars_args
  REAL(8) :: scal
  INTEGER :: iargc
  
  IF( iargc() /= 4 ) THEN 
    WRITE(*,*) 'ERROR: exactly two arguments are needed'
    STOP 
  ENDIF 

  CALL getarg( 1, chars_args )
  READ( chars_args, * ) N

  CALL getarg( 2, chars_args )
  READ( chars_args, * ) scal

  CALL getarg( 3, chars_args )
  READ( chars_args, * ) A

  CALL getarg( 4, chars_args )
  READ( chars_args, * ) alpha

  IF( scal > 1.d0 .or. scal < 0.d0 ) THEN 
    WRITE(*,*) 'scal must be 0.0 < scal < 1.0'
    STOP 
  ENDIF 

  L = 10.d0

  c1 = L*scal  ! position

  ! Initialize the basis functions
  ALLOCATE( grid_x(N) )
  CALL init_grid_1d_p( N, 0.d0, L, grid_x )

  h = L/dble(N)  ! manually calculate h

  ! Calculate coefficients of expansion for the function
  ALLOCATE( coefs(N) )

  DO ii=1,N
!    coefs(ii) = funcx( c1, L, grid_x(ii) ) 
!    coefs(ii) = funcx2( 1.d0, 15.d0, c1, grid_x(ii) ) 
    coefs(ii) = funcx3( A, alpha, c1, L, grid_x(ii) ) 
    WRITE(11,'(1x,2F20.10)') grid_x(ii), coefs(ii)
  ENDDO

  WRITE(*,'(1x,A,F18.10)') 'Integ = ', sum(coefs)*h
  
END PROGRAM 

