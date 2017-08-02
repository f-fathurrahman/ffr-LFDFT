! Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
! arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss-
! Legendre n-point quadrature formula.
SUBROUTINE gauleg( x1, x2, N, x, w )
  USE m_constants, ONLY : PI
  IMPLICIT NONE 
  !
  REAL(8) :: x1, x2
  INTEGER :: N
  REAL(8) :: x(N), w(N)
  !
  REAL(8), PARAMETER :: EPS = 3.d-11 ! Relative PRECISION
  !
  INTEGER :: m, j, i
  REAL(8) :: z1, z, xm, xl, pp, p3, p2, p1  ! High precision is a good idea for this routine.

  m = (N + 1)/2        ! The roots are symmetric in the interval, so
  xm = 0.5d0*(x2 + x1) ! we only have to find half of them.
  xl = 0.5d0*(x2 - x1)

  ! FIXME: Need this ?
  z1 = 100.d0  ! some initial number to get the while loop enter the first time
  pp = 0.d0    ! make pp visible outside for loop

  DO i = 1, m ! Loop over the desired roots.
    z = cos( pi*(i-0.25)/(n+0.5) )
    ! Starting with the above approximation to the ith root, we enter the main loop of
    ! refinement by Newton’s method.
    DO 
      p1 = 1.d0
      p2 = 0.d0
      DO j = 1, N     ! Loop up the recurrence relation to get the
        p3 = p2       ! Legendre polynomial evaluated at z.
        p2 = p1
        p1 = ( (2.d0*j - 1.d0)*z*p2 - (j - 1.d0)*p3 )/j
      ENDDO 
      ! p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
      ! by a standard relation involving also p2, the polynomial of one lower order.
      pp = N*(z*p1-p2)/(z*z-1.d0)
      z1 = z
      z  = z1 - p1/pp ! Newton’s method.
      IF( abs(z-z1) < EPS ) EXIT 
    ENDDO
    !
    x(i)     = xm - xl*z                  ! Scale the root to the desired interval,
    x(n-i+1) = xm + xl*z                  ! and put in its symmetric counterpart.
    w(i)     = 2.d0*xl/((1.d0-z*z)*pp*pp) ! Compute the weight
    w(n-i+1) = w(i)                       ! and its symmetric counterpart.
  ENDDO

END SUBROUTINE 



