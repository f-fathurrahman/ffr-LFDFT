!------------------------------------------------------------------------------
PROGRAM t_fun
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF1d
  IMPLICIT NONE
  !
  INTEGER :: N = 5
  REAL(8) :: L = 2.d0
  !
  TYPE(LF1d_t) :: LF
  !
  INTEGER, PARAMETER :: NPTS_PLOT = 101
  INTEGER :: ibf, ip
  REAL(8) :: x, f(NPTS_PLOT)
  
  ! Initialize the basis functions
  CALL init_LF1d_p( LF, N, 0.d0, L )

  CALL info_LF1d( LF, .TRUE. )

  ! Prepare data for plotting each LF
  DO ibf = 1, N
    DO ip=1,NPTS_PLOT
      x = (ip-1)*L/(NPTS_PLOT-1)
      f(ip) = eval_LF1d_p( LF, ibf, x )
      !WRITE(100+ibf,*) x, f(ip)
    ENDDO
    ! Test the normalization
    WRITE(*,'(1x,A,I8,F18.10)') 'ibf, integral = ', ibf, sum( f(:)*f(:) )*L/(NPTS_PLOT-1)
  ENDDO


END PROGRAM

