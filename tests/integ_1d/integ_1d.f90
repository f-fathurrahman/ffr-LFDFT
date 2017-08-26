! test convergence of integrals
! TODO: also test egg-box effect in 1d ???

! A periodic function
FUNCTION funcx( c1, L, x )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  REAL(8) :: funcx
  REAL(8) :: c1, L, x

  funcx = 1.d-6*exp( 25.d0*cos(2.d0*PI/L *(x - c1) ) )
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
  REAL(8) :: funcx
  !
  INTEGER :: NPTS_PLOT=401
  REAL(8) :: xx, yy, h, c1
  !
  REAL(8) :: eval_LF1d_p
  !
  CHARACTER(8) :: chars_args
  REAL(8) :: scal
  INTEGER :: iargc
  
  IF( iargc() /= 2 ) THEN 
    WRITE(*,*) 'ERROR: exactly one argument is needed'
    STOP 
  ENDIF 

  CALL getarg( 1, chars_args )
  READ( chars_args, * ) N

  CALL getarg( 2, chars_args )
  READ( chars_args, * ) scal

  IF( scal > 1.d0 .or. scal < 0.d0 ) THEN 
    WRITE(*,*) 'scal must be 0.0 < scal < 1.0'
    STOP 
  ENDIF 

  L = 10.d0

  c1 = L*0.5d0  ! position

  ! Initialize the basis functions
  ALLOCATE( grid_x(N) )
  CALL init_grid_1d_p( N, 0.d0, L, grid_x )

  h = L/dble(N)  ! manually calculate h

  ! Calculate coefficients of expansion for the function
  ALLOCATE( coefs(N) )

  DO ii=1,N
    coefs(ii) = funcx( c1, L, grid_x(ii) ) 
    WRITE(11,'(1x,2ES18.10)') grid_x(ii), coefs(ii)
  ENDDO

  WRITE(*,'(1x,A,F18.10)') 'Integ = ', sum(coefs)*h
  
END PROGRAM 

