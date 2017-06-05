#include "sinc.f90"

PROGRAM test_sinc

  IMPLICIT NONE 
  REAL(8) :: sinc

  WRITE(*,*) sinc(1.d-8)

END PROGRAM 
