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
  NN(:) = 55

  CALL init_LF3d_p( NN, AA, BB )

  num_x = NN(1)
  ALLOCATE( dat(num_x+1) )
  ALLOCATE( work(num_x+1) )

  x0 = grid_x(1)
  x1 = grid_x(1) + LL(1)
  
  x0_code = 0  ! periodic
  x1_code = 0

  work(1:num_x) = grid_x(:)
  DO ix = 1,num_x
    dat(ix) = func(grid_x(ix))  !sin( 2.d0*PI*grid_x(ix)/LL(1) )
    WRITE(201,*) grid_x(ix), dat(ix)
  ENDDO 
  dat(num_x+1) = dat(1)
  work(num_x+1) = grid_x(num_x) + ( grid_x(2) - grid_x(1) )
  WRITE(201,*) work(num_x+1), dat(num_x+1)

  CALL FCREATE_UBSPLINE_1D_D(x0, x1, num_x, x0_code, x0_val, x1_code, x1_val, dat, spline)

  !x0_val ??
  !x1_val ??

  x = 3.9d0
  CALL feval_ubspline_1d_d( spline, x, val ) 
  WRITE(*,'(1x,2F18.10,ES18.10)') x, val, abs( val - func(x) )

  x = 2.9d0
  CALL feval_ubspline_1d_d( spline, x, val ) 
  WRITE(*,'(1x,2F18.10,ES18.10)') x, val, abs( val - func(x) )

  x = 1.9d0
  CALL feval_ubspline_1d_d( spline, x, val ) 
  WRITE(*,'(1x,2F18.10,ES18.10)') x, val, abs( val - func(x) )

  CALL fdestroy_bspline( spline )

  DEALLOCATE( dat )

END PROGRAM 


FUNCTION func( x ) RESULT( res )

  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : LL => LF3d_LL
  REAL(8) :: res
  REAL(8) :: x

  res = cos( 2.d0*PI*x/LL(1) )
  !WRITE(*,*) 'x = ', x, ' res = ', res

END FUNCTION 


