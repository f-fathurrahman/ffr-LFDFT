PROGRAM ex_init
  USE m_constants
  IMPLICIT NONE
  !
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  REAL(8) :: t1, t2
  REAL(8) :: hh(3)

  NN = (/ 63, 63, 63 /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)

  !CALL init_LF3d_c( NN, AA, BB )
  !CALL init_LF3d_p( NN, AA, BB )

  hh = (/ 0.3d0, 0.3d0, 0.3d0 /)
  CALL init_LF3d_sinc( NN, hh )

  CALL info_LF3d()

  CALL dealloc_LF3d()
  
END PROGRAM

