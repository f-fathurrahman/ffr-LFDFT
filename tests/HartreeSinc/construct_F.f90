
SUBROUTINE construct_F( axis, t_size, t_values, F_values )
  USE m_LF3d, ONLY : NN => LF3d_NN, &
                     grid_x => LF3d_grid_x, &
                     grid_y => LF3d_grid_y, &
                     grid_z => LF3d_grid_z, &
                     hh => LF3d_hh
  IMPLICIT NONE 
  INTEGER :: axis
  INTEGER :: t_size
  REAL(8) :: t_values(t_size)
  REAL(8) :: F_values(t_size,NN(axis),NN(axis))
  INTEGER :: i_t, i, j
  REAL(8), ALLOCATABLE :: grid(:)
  REAL(8) :: compute_F

  ! copy grid_* to grid
  ALLOCATE( grid(NN(axis)) )
  IF( axis == 1 ) THEN 
    grid(:) = grid_x(:)
  ELSEIF( axis == 2 ) THEN 
    grid(:) = grid_y(:)
  ELSEIF( axis == 3 ) THEN 
    grid(:) = grid_z(:)
  ELSE 
    WRITE(*,*) 'Invalid value to axis', axis
    STOP 
  ENDIF 

  DO i_t = 1,t_size
    DO i = 1,NN(axis)
      DO j = 1,NN(axis)
        F_values(i_t,i,j) = compute_F( t_values(i_t), abs( grid(i) - grid(j) ), hh(axis) )
      ENDDO 
    ENDDO 
  ENDDO 

  DEALLOCATE( grid )

END SUBROUTINE 


