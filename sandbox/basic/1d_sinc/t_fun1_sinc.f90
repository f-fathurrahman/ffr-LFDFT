!------------------------------------------------------------------------------
PROGRAM t_fun
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF1d
  IMPLICIT NONE
  !
  INTEGER :: N = 4
  REAL(8) :: h = 1.d0
  !
  TYPE(LF1d_t) :: LF
  !
  REAL(8) :: Lplot
  INTEGER, PARAMETER :: NPTS_PLOT=100
  INTEGER :: ibf, ip
  REAL(8) :: x, f(NPTS_PLOT)
  
  Lplot = 4.d0

  ! Initialize the basis functions
  CALL init_LF1d_sinc( LF, N, h )

  CALL info_LF1d( LF, .TRUE. )


  DO ibf = 1, N
    DO ip=1,NPTS_PLOT
      x = -Lplot + (ip-1)*2.d0*Lplot/(NPTS_PLOT-1)
      f(ip) = eval_LF1d_sinc( LF, ibf, x )
      WRITE(100+ibf,*) x, f(ip)
    ENDDO
  ENDDO


END PROGRAM

