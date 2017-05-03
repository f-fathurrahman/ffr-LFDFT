PROGRAM test_bspline
  USE bspline
  INTEGER :: Nx
  INTEGER :: i
  REAL(8), ALLOCATABLE :: x(:), fx(:)
  REAL(8) :: Lx
  REAL(8) :: xx, fxx
  !
  REAL(8), ALLOCATABLE :: tx(:)
  INTEGER :: kx, iknot, idx
  INTEGER :: iflag(6)
  REAL(8) :: val, err, errmax, tru
  !
  REAL(8) :: fun1d

  Nx = 25
  Lx = 1.d0

  ALLOCATE( x(Nx) )
  ALLOCATE( fx(Nx) )

  DO i = 1,Nx
    x(i) = dble(i-1)/dble(Nx-1) * Lx
    fx(i) = fun1d( x(i), Lx )
    WRITE(101,*) x(i), fx(i)
  ENDDO 

  kx = 4 ! order in x
  iknot = 0 ! automatically select the knots
  ALLOCATE( tx(Nx+kx) )

  CALL db1ink( x, Nx, fx, kx, iknot, tx, fx, iflag(1) )
  !WRITE(*,*) 'iflag(1) = ', iflag(1)
  !DO i = 1,Nx
  !  WRITE(*,*) x(i), fx(i)
  !ENDDO 
  !DO i = 1,Nx+kx
  !  WRITE(*,*) i, tx(i)
  !ENDDO 

  errmax = 0.d0
  idx = 0
  inbvx = 1
  DO i = 1,2*Nx
    xx = dble(i-1)/dble(2*Nx-1)*Lx
    CALL db1val( xx, idx, tx, nx, kx, fx, fxx, iflag(1), inbvx )
    tru = fun1d( xx, Lx )
    err = abs( tru - fxx )
    errmax = max( err, errmax )
    WRITE(102,*) xx, tru, err
  ENDDO 

  DEALLOCATE( x, fx )
END PROGRAM 

FUNCTION fun1d( x, L ) RESULT(f)
  USE m_constants, ONLY : PI
  REAL(8) :: f
  REAL(8) :: x, L

  f = cos(2.d0*PI/L * x) + sin(3*2.d0*PI/L * x)

END FUNCTION 

