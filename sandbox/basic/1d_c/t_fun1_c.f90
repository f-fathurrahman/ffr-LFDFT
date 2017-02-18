!------------------------------------------------------------------------------
PROGRAM t_fun
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF1d
  IMPLICIT NONE
  !
  INTEGER :: N = 4
  REAL(8) :: L = 2.5d0
  !
  TYPE(LF1d_t) :: LF
  !
  INTEGER, PARAMETER :: NPTS_PLOT = 200
  INTEGER :: ibf, ip
  REAL(8) :: x, f(NPTS_PLOT)
  
  ! Initialize the basis functions
  CALL init_LF1d_c( LF, N, 0.d0, L )
  CALL info_LF1d( LF, .TRUE. )

  DO ibf = 1, N
    DO ip=1,NPTS_PLOT
      x = (ip-1)*L/(NPTS_PLOT-1)
      f(ip) = eval_LF1d_c( LF, ibf, x )
      WRITE(100+ibf,*) x, f(ip)
    ENDDO
    ! Test the normalization of the analytic expression of LF
    WRITE(*,'(1x,A,I8,F18.10)') 'ibf, integral = ', ibf, sum( f(:)*f(:) )*L/(NPTS_PLOT-1)
  ENDDO


END PROGRAM

