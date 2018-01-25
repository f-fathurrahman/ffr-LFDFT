!------------------------------------------------------------------------------
PROGRAM t_eval1d
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : grid_x => LF3d_grid_x
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: L
  !
  INTEGER :: ii, jj
  !
  INTEGER :: NPTS_PLOT=401
  REAL(8) :: xx, yy, h, c1
  !
  REAL(8) :: eval_LF1d_p

  N = 5
  L = 1.d0

  c1 = L/2.d0

  ! Initialize the basis functions
  ALLOCATE( grid_x(N) )
  CALL init_grid_1d_p( N, 0.d0, L, grid_x )

  h = L/dble(N)  ! manually calculate h

  ! Write to file for plotting
  DO jj = 1,N
    DO ii=1,NPTS_PLOT
      !
      xx = (ii-1)*L/NPTS_PLOT ! plot for double periods
      yy = eval_LF1d_p( N, L, grid_x, jj, xx )
      WRITE(10+jj,*) xx, yy
    ENDDO
  ENDDO 

  ! Free memory
  DEALLOCATE( grid_x )

END PROGRAM

