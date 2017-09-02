! eggbox effect for HGH pseudopot in 1d

! periodic function
FUNCTION funcx3( psp, c1, L, x )
  USE m_Ps_HGH
  IMPLICIT NONE 
  TYPE(Ps_HGH_Params_T) :: psp
  REAL(8) :: c1, x, funcx3, dx1, dx2, dx3, dx, L

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
  funcx3 = hgh_eval_Vloc_R_short( psp, dx )
END FUNCTION 


PROGRAM test_integral
  USE m_constants, ONLY : PI
  USE m_Ps_HGH
  USE m_LF3d, ONLY : grid_x => LF3d_grid_x
  IMPLICIT NONE
  !
  TYPE(Ps_HGH_Params_T) :: ps
  INTEGER :: N
  REAL(8) :: L
  !
  REAL(8), ALLOCATABLE :: coefs(:)
  !
  INTEGER :: ii, jj
  ! Functions
  REAL(8) :: funcx, funcx2, funcx3
  !
  CHARACTER(56) :: filpsp
  REAL(8) :: xx, yy, h, c1
  !
  REAL(8) :: eval_LF1d_p
  !
  CHARACTER(8) :: chars_args
  REAL(8) :: scal
  INTEGER :: iargc
  
  IF( iargc() /= 3 ) THEN 
    WRITE(*,*) 'ERROR: exactly three arguments are needed'
    STOP 
  ENDIF 

  CALL getarg( 1, chars_args )
  READ( chars_args, * ) N

  CALL getarg( 2, chars_args )
  READ( chars_args, * ) scal

  CALL getarg(3,filpsp)
  CALL init_Ps_HGH_Params( ps, filpsp )

  IF( scal > 1.d0 .or. scal < 0.d0 ) THEN 
    WRITE(*,*) 'scal must be 0.0 < scal < 1.0'
    STOP 
  ENDIF 

  L = 16.d0

  c1 = L*scal  ! position

  ! Initialize the basis functions
  ALLOCATE( grid_x(N) )
  CALL init_grid_1d_p( N, 0.d0, L, grid_x )

  h = L/dble(N)  ! manually calculate h

  ! Calculate coefficients of expansion for the function
  ALLOCATE( coefs(N) )

  DO ii=1,N
    coefs(ii) = funcx3( ps, c1, L, grid_x(ii) ) 
    WRITE(11,'(1x,2F20.10)') grid_x(ii), coefs(ii)
  ENDDO

  WRITE(*,'(1x,2F18.10)') scal, sum(coefs)*h
  
END PROGRAM 

