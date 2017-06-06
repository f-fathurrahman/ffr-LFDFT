!------------------------------------------------------------------------------
PROGRAM test_V_ps_loc
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : grid_x => LF3d_grid_x
  USE m_Ps_HGH, ONLY : Ps_HGH_Params_T, hgh_eval_Vloc_R_short, init_Ps_HGH_Params
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
  REAL(8) :: xx, yy, h, A
  REAL(8) :: center
  !
  REAL(8) :: eval_LF1d_c, eval_LF1d_p
  TYPE(Ps_HGH_Params_T) :: Ps

  A = 0.d0
  N = 45
  L = 16.d0

  ! Initialize the basis functions
  ALLOCATE( grid_x(N) )
  CALL init_grid_1d_c( N, A, L, grid_x )

  h = L/dble(N)  ! manually calculate h

  ! Calculate coefficients of expansion for the function
  ALLOCATE( coefs(N) )

  ! 
  CALL init_Ps_HGH_Params( Ps, '../../HGH/C.hgh' )
  center = L/2.d0

  DO ii=1,N
    coefs(ii) = sqrt(h)* hgh_eval_Vloc_R_short( Ps, grid_x(ii)-center )
    WRITE(11,*) grid_x(ii), coefs(ii)/sqrt(h)
  ENDDO

  ! Write to file for plotting
  DO ii=1,NPTS_PLOT
    !
    xx = (ii-1)*L/NPTS_PLOT
    yy = 0.d0
    DO jj=1,N
      yy = yy + eval_LF1d_c( N, L, A, grid_x, jj, xx )*coefs(jj)
    ENDDO
    WRITE(12,*) xx
    WRITE(13,*) yy
    WRITE(14,*) hgh_eval_Vloc_R_short( Ps, xx-center )
  ENDDO

  ! Free memory
  DEALLOCATE( grid_x )
  DEALLOCATE( coefs )


END PROGRAM

