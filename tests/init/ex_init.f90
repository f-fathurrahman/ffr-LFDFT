PROGRAM ex_init
  USE m_constants
  IMPLICIT NONE
  !
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)

  NN = (/ 31, 31, 29 /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)

  CALL init_LF3d_p( NN, AA, BB )

  CALL info_LF3d()

  CALL dealloc_LF3d()
  
END PROGRAM

