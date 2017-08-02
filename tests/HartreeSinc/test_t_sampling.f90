PROGRAM test_t_sampling

  IMPLICIT NONE 
  INTEGER :: num_points1, num_points2
  INTEGER :: t_size, i
  REAL(8) :: t_i, t_l, t_f
  REAL(8), ALLOCATABLE :: t_values(:), w_t(:)

  num_points1 = 5
  num_points2 = 10

  t_size = num_points1 + num_points2

  t_i = 0.d0
  t_l = 2.d0
  t_f = 10000.d0

  ALLOCATE( t_values(t_size), w_t(t_size) )
  CALL HartreeSinc_t_sampling( num_points1, num_points2, t_i, t_l, t_f, &
                               t_values, w_t )
  DO i = 1, t_size
    WRITE(*,'(F18.10,1x,F18.10)') t_values(i), w_t(i)
  ENDDO 

  DEALLOCATE( t_values, w_t )
END PROGRAM 

