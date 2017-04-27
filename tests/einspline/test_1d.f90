PROGRAM test

  USE m_constants, ONLY : PI
  IMPLICIT NONE 
  REAL(8) :: x0, x1
  INTEGER :: num_x  ! no of grid points
  INTEGER :: x0_code, x1_code ! BC types
  REAL(8) :: x0_val, x1_val ! BC vals
  REAL(8), ALLOCATABLE :: dat(:), grid_x(:)
  INTEGER(8) :: spline
  !
  INTEGER :: ip, ix
  REAL(8) :: AA(3), BB(3)
  INTEGER :: NN(3)
  REAL(8) :: x, val, delta, L

  AA(:) = 0.d0
  BB(:) = 4.d0
  NN(:) = 45

  L = BB(1) - AA(1)

  num_x = NN(1)
  ALLOCATE( dat(num_x) )
  ALLOCATE( grid_x(num_x) )

  DO ix = 1,num_x
    grid_x(ix) = AA(1) + dble(ix-1)*( BB(1) - AA(1) )/dble(num_x)
    dat(ix) = sin( grid_x(ix)*2.d0*PI/L ) + sin( 5.d0*grid_x(ix)*2.d0*PI/L )
    WRITE(*,*) grid_x(ix), dat(ix)
  ENDDO 
  delta = ( BB(1) - AA(1) )/num_x

  x0 = grid_x(1)
  x1 = grid_x(num_x)
  
  x0_code = 0  ! periodic
  x1_code = 0

  x0_val = dat(1)
  x1_val = dat(num_x)

  CALL FCREATE_UBSPLINE_1D_D (x0, x1, num_x, x0_code, x0_val, x1_code, x1_val, dat, spline)

  DO ix = 1,num_x
    x = grid_x(ix) + 0.5d0*delta
    CALL FEVAL_UBSPLINE_1D_D( spline, x, val ) 
    WRITE(*,*) grid_x(ix), dat(ix)
    WRITE(*,*) x, val
  ENDDO 

  CALL FDESTROY_BSPLINE( spline )

  DEALLOCATE( dat )
END PROGRAM 

