!---------------------------
SUBROUTINE r8_rand_vec(N, v)
!---------------------------
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: v(N)
  !
  INTEGER :: i

  DO i=1,N
    CALL random_number( v(i) )
  ENDDO

END SUBROUTINE


