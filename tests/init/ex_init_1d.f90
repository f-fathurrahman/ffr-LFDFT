SUBROUTINE test_sinc( N, h )
  USE m_LF3d, ONLY : grid_x => LF3d_grid_x
  IMPLICIT NONE
  INTEGER :: N
  REAL(8) :: h
  INTEGER :: i

  ! Initialize the basis functions
  ALLOCATE( grid_x(N) )
  CALL init_grid_1d_sinc( N, h, grid_x ) 

  WRITE(*,*)
  WRITE(*,'(1x,A,I5)') 'N = ', N
  WRITE(*,'(1x,A,F18.10)') 'L = ', grid_x(N) - grid_x(1)
  WRITE(*,'(1x,A,F18.10)') 'Grid spacing = ', grid_x(2) - grid_x(1)
  WRITE(*,*)
  WRITE(*,*) 'Grid points for sinc LF'
  WRITE(*,*)
  DO i = 1,N
    WRITE(*,'(1x,I5,F18.10)') i, grid_x(i)
  ENDDO 

  DEALLOCATE( grid_x )
END SUBROUTINE 


SUBROUTINE test_c( N, L )
  USE m_LF3d, ONLY : grid_x => LF3d_grid_x
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: L
  INTEGER :: i

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

  DEALLOCATE( grid_x )
END SUBROUTINE 



SUBROUTINE test_p( N, L )
  USE m_LF3d, ONLY : grid_x => LF3d_grid_x
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: L
  INTEGER :: i

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

  DEALLOCATE( grid_x )
END SUBROUTINE 


PROGRAM ex_init_1d

  CALL test_p( 5, 5.d0 )
  CALL test_c( 5, 5.d0 )
  CALL test_sinc( 5, 5.d0/4 )

END PROGRAM 

