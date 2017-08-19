SUBROUTINE test_sinc( N, h, NPTS_PLOT )
  USE m_LF3d, ONLY : grid_x => LF3d_grid_x
  IMPLICIT NONE
  INTEGER :: N
  REAL(8) :: h, L
  INTEGER :: i, ibf
  INTEGER :: NPTS_PLOT
  !
  REAL(8) :: x, dx, f
  REAL(8) :: eval_LF1d_sinc
  INTEGER :: ierr

  ! Initialize the basis functions
  ALLOCATE( grid_x(N) )
  CALL init_grid_1d_sinc( N, h, grid_x ) 
  L = grid_x(N) - grid_x(1)
  WRITE(*,*)
  WRITE(*,'(1x,A,I5)') 'N = ', N
  WRITE(*,'(1x,A,F18.10)') 'L = ', L
  WRITE(*,'(1x,A,F18.10)') 'Grid spacing = ', grid_x(2) - grid_x(1)
  WRITE(*,*)
  WRITE(*,*) 'Grid points for sinc LF'
  WRITE(*,*)
  DO i = 1,N
    WRITE(*,'(1x,I5,F18.10)') i, grid_x(i)
  ENDDO 

  ! Plot
  OPEN( unit=100, file='LF1d_sinc.dat', action='write', iostat=ierr )
  dx = L/(NPTS_PLOT-1.d0)
  WRITE(*,*) 'dx = ', dx
  ! loop over number of points to plot
  DO i = 1,NPTS_PLOT
    x = -0.5d0*L + (i-1)*dx
    WRITE(100,'(1x,F18.10)',advance='no') x
    ! loop over number of LFs
    DO ibf = 1,N
      f = eval_LF1d_sinc( N, grid_x, ibf, x ) 
      WRITE(100,'(1x,F18.10)',advance='no') f
    ENDDO 
    WRITE(100,*)
  ENDDO 
  CLOSE(100)
  DEALLOCATE( grid_x )
END SUBROUTINE 


SUBROUTINE test_c( N, L, NPTS_PLOT )
  USE m_LF3d, ONLY : grid_x => LF3d_grid_x
  IMPLICIT NONE 
  INTEGER :: N, NPTS_PLOT
  REAL(8) :: L
  INTEGER :: i, ibf
  !
  REAL(8) :: x, dx, f
  REAL(8) :: eval_LF1d_c
  INTEGER :: ierr

  ! Initialize the basis functions
  ALLOCATE( grid_x(N) )
  CALL init_grid_1d_c( N, -0.5d0*L, 0.5d0*L, grid_x ) 

  WRITE(*,*)
  WRITE(*,'(1x,A,I5)') 'N = ', N
  WRITE(*,'(1x,A,F18.10)') 'L = ', L
  WRITE(*,'(1x,A,F18.10)') 'Grid spacing = ', grid_x(2) - grid_x(1)
  WRITE(*,*)
  WRITE(*,*) 'Grid points for cluster/box LF'
  WRITE(*,*)
  DO i = 1,N
    WRITE(*,'(1x,I5,F18.10)') i, grid_x(i)
  ENDDO 

  ! Plot
  OPEN( unit=100, file='LF1d_c.dat', action='write', iostat=ierr )
  dx = L/(NPTS_PLOT-1.d0)
  WRITE(*,*) 'dx = ', dx
  ! loop over number of points to plot
  DO i = 1,NPTS_PLOT
    x = -0.5d0*L + (i-1)*dx
    WRITE(100,'(1x,F18.10)',advance='no') x
    ! loop over number of LFs
    DO ibf = 1,N
      f = eval_LF1d_c( N, L, -0.5d0*L, grid_x, ibf, x ) 
      WRITE(100,'(1x,F18.10)',advance='no') f
    ENDDO 
    WRITE(100,*)
  ENDDO 
  CLOSE(100)

  DEALLOCATE( grid_x )
END SUBROUTINE 



SUBROUTINE test_p( N, L, NPTS_PLOT )
  USE m_LF3d, ONLY : grid_x => LF3d_grid_x
  IMPLICIT NONE 
  INTEGER :: N, NPTS_PLOT
  REAL(8) :: L
  INTEGER :: i, ibf
  !
  REAL(8) :: x, dx, f
  REAL(8) :: eval_LF1d_p
  INTEGER :: ierr

  ! Initialize the basis functions
  ALLOCATE( grid_x(N) )
  CALL init_grid_1d_p( N, -0.5d0*L, 0.5d0*L, grid_x ) 

  WRITE(*,*)
  WRITE(*,'(1x,A,I5)') 'N = ', N
  WRITE(*,'(1x,A,F18.10)') 'L = ', L
  WRITE(*,'(1x,A,F18.10)') 'Grid spacing = ', grid_x(2) - grid_x(1)
  WRITE(*,*)
  WRITE(*,*) 'Grid points for periodic LF'
  WRITE(*,*)
  DO i = 1,N
    WRITE(*,'(1x,I5,F18.10)') i, grid_x(i)
  ENDDO 

  ! Plot
  OPEN( unit=100, file='LF1d_p.dat', action='write', iostat=ierr )
  dx = L/(NPTS_PLOT-1.d0)
  WRITE(*,*) 'dx = ', dx
  ! loop over number of points to plot
  DO i = 1,NPTS_PLOT
    x = -0.5d0*L + (i-1)*dx
    WRITE(100,'(1x,F18.10)',advance='no') x
    ! loop over number of LFs
    DO ibf = 1,N
      f = eval_LF1d_p( N, L, grid_x, ibf, x ) 
      WRITE(100,'(1x,F18.10)',advance='no') f
    ENDDO 
    WRITE(100,*)
  ENDDO 
  CLOSE(100)

  DEALLOCATE( grid_x )
END SUBROUTINE 


PROGRAM ex_init_1d

  CALL test_p( 5, 5.d0, 100 )
  CALL test_c( 5, 5.d0, 100 )
  CALL test_sinc( 5, 5.d0/4, 100 )

END PROGRAM 

