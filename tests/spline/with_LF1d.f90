PROGRAM test
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : grid_x => LF3d_grid_x, LL => LF3d_LL
  IMPLICIT NONE 
  REAL(8) :: x0, x1
  INTEGER :: num_x  ! no of grid points
  INTEGER :: x0_code, x1_code ! BC types
  REAL(8) :: x0_val, x1_val ! BC vals
  REAL(8), ALLOCATABLE :: dat(:), work(:)
  INTEGER(8) :: spline
  !
  INTEGER :: ip, ix
  REAL(8) :: AA(3), BB(3)
  INTEGER :: NN(3)
  REAL(8) :: x, val, delta
  REAL(8) :: func

  AA(:) = 0.d0
  BB(:) = 4.d0
  NN(:) = 75

  CALL init_LF3d_p( NN, AA, BB )

  num_x = NN(1)
  ALLOCATE( dat(num_x+1) )
  ALLOCATE( work(num_x+1) )

  x0 = grid_x(1)
  x1 = grid_x(1) + LL(1)
  
  x0_code = 0  ! periodic
  x1_code = 0

  delta = grid_x(2)-grid_x(1)

  work(1:num_x) = grid_x(:)
  DO ix = 1,num_x
    dat(ix) = func(grid_x(ix))
    WRITE(201,'(2F22.12)') grid_x(ix), dat(ix)
  ENDDO 
  dat(num_x+1) = dat(1)  ! value at the end point is the same as the value at start point
  work(num_x+1) = grid_x(num_x) + delta  ! the end points
  WRITE(201,'(2F22.12)') work(num_x+1), dat(num_x+1)

  CALL FCREATE_UBSPLINE_1D_D(x0, x1, num_x, x0_code, x0_val, x1_code, x1_val, dat, spline)

  !x0_val and x1_val are not used ??

  ! Now evaluate the value at FFT / FD grid
  DO ix = 1,num_x
    x = ( grid_x(1) + delta*0.5d0 ) + (ix-1)*delta
    CALL feval_ubspline_1d_d( spline, x, val ) 
    WRITE(*,'(1x,2F18.10,ES18.10)') x, val, abs( val - func(x) )
    WRITE(202,'(2F22.12)') x, val
  ENDDO 

  CALL fdestroy_bspline( spline )

  DEALLOCATE( dat )

END PROGRAM 


FUNCTION func( x ) RESULT( res )

  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : LL => LF3d_LL
  REAL(8) :: res
  REAL(8) :: x
  INTEGER, PARAMETER :: N = 5
  INTEGER :: i

  res = 0.d0
  DO i = 1, N
    res = res + cos( i*2.d0*PI*x/LL(1) )
  ENDDO 

END FUNCTION 


