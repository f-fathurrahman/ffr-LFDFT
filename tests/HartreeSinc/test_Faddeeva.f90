PROGRAM test_Faddeeva

  IMPLICIT NONE 
  REAL(8) :: f_re, f_im
  COMPLEX(8) :: z

  CALL Cwrap_faddeeva( 1.d0, 1.d0, f_re, f_im )
  WRITE(*,*)
  WRITE(*,'(1x,2F18.10)') f_re, f_im

END PROGRAM 

