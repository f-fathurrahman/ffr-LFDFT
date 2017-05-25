SUBROUTINE read_control(filein)
  USE m_input_vars
  IMPLICIT NONE 
  CHARACTER(*) :: filein

  ! default values
  pseudo_dir = './'
  etot_conv_thr = 1d-6  ! in Ha instead of Ry !!!

  OPEN(IU,file=filein, status='old')
  READ(IU,nml=CONTROL)
  CLOSE(IU)

  WRITE(*,*) 'pseudo_dir = ', trim(pseudo_dir)
  WRITE(*,'(1x,A,ES18.10)') 'etot_conv_thr = ', etot_conv_thr

END SUBROUTINE 

