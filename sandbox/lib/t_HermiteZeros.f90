PROGRAM t_HermiteZeros
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: N = 11
  REAL(8) :: zeros(n)
  INTEGER :: ip

  CALL h_polynomial_zeros( n, zeros )
  DO ip=1,N
    WRITE(*,'(I8,F18.10)') ip, zeros(ip)
  ENDDO

END PROGRAM

