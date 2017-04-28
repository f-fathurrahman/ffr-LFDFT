PROGRAM test_limit

  IMPLICIT NONE
  REAL(8) :: eps 

  eps = epsilon(1.d0)
  WRITE(*,*) 'eps = ', eps
  WRITE(*,*) 'limit = ', erf(eps)/eps

END PROGRAM 

