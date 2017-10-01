! currently only for writing evecs to 
PROGRAM postproc
  IMPLICIT NONE 
  
  CALL write_KS_evecs_xsf()
END PROGRAM 


! wrapper to 
SUBROUTINE write_KS_evecs_xsf()
  USE m_checkpoint
  USE m_io_data
  USE m_states, ONLY : KS_evecs, Nstates
  IMPLICIT NONE 
  INTEGER :: ist
  CHARACTER(56) :: filename, str_ist

  CALL read_checkpoint()
  CALL read_KS_evecs('KS_evecs.dat')

  DO ist = 1, Nstates
    IF( ist < 10 ) THEN 
      WRITE(str_ist,'(I1)') ist
    ELSEIF( ist < 100 ) THEN 
      WRITE(str_ist,'(I2)') ist
    ELSEIF( ist < 1000 ) THEN 
      WRITE(str_ist,'(I3)') ist
    ENDIF 
    filename = 'evecs_'//trim(str_ist)//'.xsf'
    CALL write_data3d_xsf(KS_evecs(:,ist), filename)
  ENDDO
END SUBROUTINE 

