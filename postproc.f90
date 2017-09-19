! currently only for writing evecs to 
PROGRAM postproc
  IMPLICIT NONE 
  
  CALL write_KS_evecs_xsf()
END PROGRAM 

SUBROUTINE write_KS_evecs_xsf()
  USE m_checkpoint
  USE m_io_data
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_states, ONLY : KS_evecs, Nstates
  IMPLICIT NONE 

  CALL
  CALL read_KS_evecs('KS_evecs.dat')

  WRITE(*,*) 'Npoints = ', Npoints
  WRITE(*,*) 'Nstates = ', Nstates
END SUBROUTINE 

