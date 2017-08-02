PROGRAM test_gauleg
  IMPLICIT NONE 
  REAL(8) :: x1, x2
  INTEGER :: N
  REAL(8), ALLOCATABLE :: x(:), w(:)
  INTEGER :: i

  x1 = 0.d0
  x2 = 2.d0
  N = 11
  ALLOCATE( x(N) )
  ALLOCATE( w(N) )
  CALL gauleg(x1, x2, N, x, w)
  
  DO i = 1,N
    WRITE(*,'(I4,1x,F18.10,1x,F18.10)') i, x(i), w(i)
  ENDDO 

  DEALLOCATE( x, w )
END PROGRAM 
