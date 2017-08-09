PROGRAM test_Faddeeva
  use Faddeeva
  IMPLICIT NONE 
  COMPLEX(8) :: z
  real(8) :: relerr

  WRITE(*,*)
  WRITE(*,'(1x,2F18.10)') erfcx( cmplx(2.d0,1.1d0,kind=8) )
  !
  WRITE(*,*)
  WRITE(*,'(1x,2F18.10)') erfcx( cmplx(2.d0,-1.1d0,kind=8) )
  !
  WRITE(*,*)
  WRITE(*,'(1x,2F18.10)') erfcx( cmplx(2.d0,0.0d0,kind=8) )
END PROGRAM 

