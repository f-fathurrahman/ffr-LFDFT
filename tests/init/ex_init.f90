PROGRAM ex_init
  USE m_constants
  IMPLICIT NONE
  !
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  REAL(8) :: t1, t2

  NN = (/ 63, 63, 63 /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)

  CALL init_LF3d_p( NN, AA, BB )

  CALL info_LF3d()

  CALL cpu_time(t1)
  CALL init_nabla2_sparse( .TRUE. )
  CALL cpu_time(t2)
  WRITE(*,*) 'Time for init_nabla2_sparse: ', t2-t1, ' seconds'
  !CALL calc_nabla2_NNZ()

  CALL dealloc_nabla2_sparse()
  CALL dealloc_LF3d()
  
END PROGRAM
