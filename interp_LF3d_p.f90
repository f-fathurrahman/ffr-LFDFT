FUNCTION interp_LF3d_p( NN, LL, lingrid, lin2xyz, coefs, point ) RESULT(f)

  IMPLICIT NONE 
  !
  INTEGER :: NN(3)
  REAL(8) :: LL(3)
  REAL(8) :: lingrid(3,*)
  INTEGER :: lin2xyz(3,*)
  REAL(8) :: coefs(*)
  REAL(8) :: point(3)
  REAL(8) :: f
  !
  INTEGER :: ip, i, j, k, Npoints
  REAL(8) :: eval_LF1d_p

  Npoints = NN(1)*NN(2)*NN(3)

  f = 0.d0
  DO ip = 1,Npoints
    i = lin2xyz(1,ip)
    j = lin2xyz(2,ip)
    k = lin2xyz(3,ip)
    f = f + eval_LF1d_p( NN(1), LL(1), lingrid(1,ip), i, point(1) ) * &
            eval_LF1d_p( NN(2), LL(2), lingrid(2,ip), j, point(2) ) * &
            eval_LF1d_p( NN(3), LL(3), lingrid(3,ip), k, point(3) ) * coefs(ip)
  ENDDO 

END FUNCTION 


