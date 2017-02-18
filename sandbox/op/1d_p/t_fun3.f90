FUNCTION soft_coul( Q, A, L, x )
  IMPLICIT NONE
  !
  REAL(8) :: soft_coul
  REAL(8) :: Q, A, L, x

  soft_coul = -Q**2 / sqrt( A**2 + (x - L/2.d0)**2 )
END FUNCTION

!
FUNCTION dsoft_coul( Q, A, L, x )
  IMPLICIT NONE
  !
  REAL(8) :: dsoft_coul
  REAL(8) :: Q, A, L, x

  dsoft_coul = -Q**2*(L/2.d0 - x)/(A**2 + (-L/2.d0 + x)**2)**(1.5d0)
END FUNCTION

FUNCTION d2soft_coul( Q, A, L,x )
  IMPLICIT NONE
  !
  REAL(8) :: d2soft_coul
  REAL(8) :: Q, A, L, x

 d2soft_coul = Q**2/(A**2 + (-L/2 + x)**2)**(1.5d0) - &
     Q**2*(L/2.d0 - x)*(3.d0*L/2.d0 - 3*x)/(A**2 + (-L/2.d0 + x)**2)**(2.5d0)
END FUNCTION


!------------------------------------------------------------------------------
PROGRAM t_fun
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF1d
  IMPLICIT NONE
  !
  INTEGER :: N = 201
  REAL(8) :: L = 5.d0
  !
  TYPE(LF1d_t) :: LF
  REAL(8), ALLOCATABLE :: coefs(:), dcoefs(:), d2coefs(:)
  !
  INTEGER :: ii, jj
  ! Functions
  REAL(8) :: soft_coul, dsoft_coul, d2soft_coul
  !
  REAL(8) :: Q, A
  INTEGER :: NPTS_PLOT=201
  REAL(8) :: xx, yy

  Q = 1.d0
  A = 1.d0

  ! Initialize the basis functions
  CALL init_LF1d_p( LF, N, 0.d0, L )

  ! Calculate coefficients of expansion for the function
  ALLOCATE( coefs(N) )
  ALLOCATE( dcoefs(N) )
  ALLOCATE( d2coefs(N) )

  DO ii=1,N
    coefs(ii) = sqrt(L/N)* soft_coul( Q, A, L, LF%grid(ii) ) 
    WRITE(11,*) LF%grid(ii), sqrt(N/L)*coefs(ii)
  ENDDO

  dcoefs = matmul( LF%D1jl, coefs )
  d2coefs = matmul( LF%D2jl, coefs )

  ! Write to file for plotting
  DO ii=1,NPTS_PLOT
    !
    !xx = (ii-1)*L/NPTS_PLOT
    xx = LF%grid(ii)  ! only works for NPTS_PLOT = N
    yy = 0.d0
    DO jj=1,N
      yy = yy + eval_LF1d_p( LF, jj, xx )*coefs(jj)
    ENDDO
    WRITE(12,*) xx
    WRITE(13,*) yy
    WRITE(14,*) soft_coul( Q, A, L, xx )
    !
    ! First derivative
    !
    yy = 0.d0
    DO jj=1,N
      !yy = yy + eval_LF1d_p_d1( LF, jj, xx )*coefs(jj)
      yy = yy + eval_LF1d_p( LF, jj, xx )*dcoefs(jj)
    ENDDO
    WRITE(15,*) yy
    WRITE(16,*) dsoft_coul( Q, A, L, xx )
    !
    ! Second derivative
    !
    yy = 0.d0
    DO jj=1,N
      yy = yy + eval_LF1d_p( LF, jj, xx )*d2coefs(jj)
    ENDDO
    WRITE(17,*) yy
    WRITE(18,*) d2soft_coul( Q, A, L, xx )
  ENDDO

  ! Free memory
  CALL dealloc_LF1d( LF )
  DEALLOCATE( coefs, dcoefs, d2coefs )


END PROGRAM

