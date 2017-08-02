PROGRAM test_compute_F
  IMPLICIT NONE 
  REAL(8) :: t, x_bar, h, F
  REAL(8) :: compute_F

  t = 0.1d0
  x_bar = 0.12d0
  h = 0.2d0
  F = compute_F(t, x_bar, h)

  WRITE(*,'(1x,A,F18.10)') 'F = ', F
END PROGRAM 
