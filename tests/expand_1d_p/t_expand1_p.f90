FUNCTION funcx( c1, L, x )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  REAL(8) :: funcx
  REAL(8) :: c1, L, x

  funcx = exp( 1.5*cos(2.d0*PI/L *(x - c1) ) )
END FUNCTION


!------------------------------------------------------------------------------
PROGRAM t_expand
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : grid_x => LF3d_grid_x
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: L
  !
  REAL(8), ALLOCATABLE :: coefs(:)
  !
  INTEGER :: ii, jj
  ! Functions
  REAL(8) :: funcx
  !
  INTEGER :: NPTS_PLOT=401
  REAL(8) :: xx, yy, h, c1
  !
  REAL(8) :: eval_LF1d_p

  N = 15
  L = 10.d0

  c1 = L/2.d0

  ! Initialize the basis functions
  ALLOCATE( grid_x(N) )
  CALL init_grid_1d_p( N, 0.d0, L, grid_x )

  h = L/dble(N)  ! manually calculate h

  ! Calculate coefficients of expansion for the function
  ALLOCATE( coefs(N) )

  DO ii=1,N
    coefs(ii) = sqrt(h)* funcx( c1, L, grid_x(ii) ) 
    WRITE(11,*) grid_x(ii), coefs(ii)/sqrt(h)
  ENDDO

  ! Write to file for plotting
  DO ii=1,NPTS_PLOT
    !
    xx = (ii-1)*2*L/NPTS_PLOT ! plot for double periods
    yy = 0.d0
    DO jj=1,N
      yy = yy + eval_LF1d_p( N, L, grid_x, jj, xx )*coefs(jj)
    ENDDO
    WRITE(12,*) xx
    WRITE(13,*) yy
    WRITE(14,*) funcx( c1, L, xx )
  ENDDO

  ! Free memory
  DEALLOCATE( grid_x )
  DEALLOCATE( coefs )


END PROGRAM

