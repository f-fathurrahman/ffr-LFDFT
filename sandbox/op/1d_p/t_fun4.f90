FUNCTION pw_Gx( Gx, x )
  IMPLICIT NONE
  !
  REAL(8) :: pw_Gx
  REAL(8) :: Gx, x

  pw_Gx = cos( Gx*x )
END FUNCTION

!
FUNCTION dpw_Gx( Gx, x )
  IMPLICIT NONE
  !
  REAL(8) :: dpw_Gx
  REAL(8) :: Gx, x

  dpw_Gx = -Gx * sin( Gx*x )
END FUNCTION

FUNCTION d2pw_Gx( Gx, x )
  IMPLICIT NONE
  !
  REAL(8) :: d2pw_Gx
  REAL(8) :: Gx, x

 d2pw_Gx = -Gx**2 * cos( Gx*x )
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
  REAL(8) :: pw_Gx, dpw_Gx, d2pw_Gx
  !
  REAL(8) :: Gx
  INTEGER :: NPTS_PLOT=201
  REAL(8) :: xx, yy

  Gx = 2*PI*1.d0/L

  ! Initialize the basis functions
  CALL init_LF1d_p( LF, N, 0.d0, L )

  ! Calculate coefficients of expansion for the function
  ALLOCATE( coefs(N) )
  ALLOCATE( dcoefs(N) )
  ALLOCATE( d2coefs(N) )

  DO ii=1,N
    coefs(ii) = sqrt(L/N)* pw_Gx( Gx, LF%grid(ii) ) 
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
    WRITE(14,*) pw_Gx( Gx, xx )
    !
    ! First derivative
    !
    yy = 0.d0
    DO jj=1,N
      yy = yy + eval_LF1d_p( LF, jj, xx )*dcoefs(jj)
    ENDDO
    WRITE(15,*) yy
    WRITE(16,*) dpw_Gx( Gx, xx )
    !
    ! Second derivative
    !
    yy = 0.d0
    DO jj=1,N
      yy = yy + eval_LF1d_p( LF, jj, xx )*d2coefs(jj)
    ENDDO
    WRITE(17,*) yy
    WRITE(18,*) d2pw_Gx( Gx, xx )
  ENDDO

  ! Free memory
  CALL dealloc_LF1d( LF )
  DEALLOCATE( coefs, dcoefs, d2coefs )


END PROGRAM

