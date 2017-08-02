! New method by Sundholm (JCP 132, 024102, 2010).
! Within two-divided regions ([t_i,t_l], [t_l,t_f]),
! num_points1, num_points2 quadrature points are made, respectively.

SUBROUTINE HartreeSinc_t_sampling( &
    num_points1, num_points2, t_i, t_l, t_f, t_values, w_t )

  USE m_constants, ONLY : PI
  IMPLICIT NONE 
  INTEGER :: num_points1, num_points2
  REAL(8) :: t_i, t_l, t_f
  REAL(8) :: t_values(num_points1 + num_points2)
  REAL(8) :: w_t(num_points1 + num_points2)
  !
  INTEGER :: j
  REAL(8) :: s_p, w_p
  REAL(8), ALLOCATABLE :: x_leg(:), w_leg(:)

  ! Linear coord region:  [t_i, t_l]
  ALLOCATE( x_leg(num_points1), w_leg(num_points1) )
  CALL gauleg(t_i, t_l, num_points1, x_leg, w_leg)
  !
  DO j=1,num_points1
    t_values(j) = x_leg(j)
    w_t(j)      = w_leg(j)*2.d0/sqrt(PI)
  ENDDO 

  DEALLOCATE( x_leg, w_leg )

  ! Logarithmic coord region: [t_l, t_f]
  ALLOCATE( x_leg(num_points2), w_leg(num_points2) )
  CALL gauleg( log(t_l), log(t_f), num_points2, x_leg, w_leg )

  ! Return the log-coord-partitioned points back to linear t-space.
  s_p = 0.d0
  w_p = 0.d0
  DO j = 1,num_points2
    s_p = x_leg(j)
    w_p = w_leg(j)
    x_leg(j) = exp(s_p)
    w_leg(j) = w_p * exp(s_p)
  ENDDO 

  DO j = 1,num_points2
    t_values(num_points1+j) = x_leg(j)
    w_t(num_points1+j)      = w_leg(j)*2.d0/sqrt(PI)
  ENDDO 

  DEALLOCATE( x_leg, w_leg )

END SUBROUTINE 


