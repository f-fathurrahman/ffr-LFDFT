FUNCTION funcx( c1, L, x )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  REAL(8) :: funcx
  REAL(8) :: c1, L, x

  funcx = exp( 20.5*cos(2.d0*PI/L *(x - c1) ) )
END FUNCTION


!------------------------------------------------------------------------------
PROGRAM t_expand
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF1d
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: L
  !
  TYPE(LF1d_t) :: LF
  REAL(8), ALLOCATABLE :: coefs(:)
  !
  INTEGER :: ii, jj
  ! Functions
  REAL(8) :: funcx
  !
  INTEGER :: NPTS_PLOT=401
  REAL(8) :: xx, yy, h, c1

  N = 15
  L = 10.d0

  c1 = L/2.d0

  ! Initialize the basis functions
  CALL init_LF1d_p( LF, N, 0.d0, L )

  h = LF%h

  ! Calculate coefficients of expansion for the function
  ALLOCATE( coefs(N) )

  DO ii=1,N
    coefs(ii) = sqrt(h)* funcx( c1, L, LF%grid(ii) ) 
    WRITE(11,*) LF%grid(ii), coefs(ii)/sqrt(h)
  ENDDO

  ! Write to file for plotting
  DO ii=1,NPTS_PLOT
    !
    xx = (ii-1)*2*L/NPTS_PLOT ! plot for double periods
    yy = 0.d0
    DO jj=1,N
      yy = yy + eval_LF1d_p( LF, jj, xx )*coefs(jj)
    ENDDO
    WRITE(12,*) xx
    WRITE(13,*) yy
    WRITE(14,*) funcx( c1, L, xx )
  ENDDO

  ! Free memory
  CALL dealloc_LF1d( LF )
  DEALLOCATE( coefs )


END PROGRAM

