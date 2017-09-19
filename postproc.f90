! currently only for writing evecs to 
PROGRAM postproc
  IMPLICIT NONE 
  
  CALL write_KS_evecs_xsf()
END PROGRAM 


! wrapper to 
SUBROUTINE write_KS_evecs_xsf()
  USE m_checkpoint
  USE m_io_data
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_states, ONLY : KS_evecs, Nstates
  IMPLICIT NONE 

  CALL read_checkpoint()
  CALL read_KS_evecs('KS_evecs.dat')

  WRITE(*,*) 'Npoints = ', Npoints
  WRITE(*,*) 'Nstates = ', Nstates

  CALL write_data3d_xsf(KS_evecs(:,1), 'evecs.xsf')
END SUBROUTINE 

