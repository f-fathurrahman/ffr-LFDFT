! various subroutines to write binary data

SUBROUTINE write_KS_evecs(filname)
  USE m_states, ONLY : KS_evecs
  IMPLICIT NONE 
  CHARACTER(*) :: filname
  INTEGER, PARAMETER :: IU = 55

  OPEN( unit=55, file=filname , action='write', form='unformatted' )
  WRITE(IU) KS_evecs
  CLOSE(IU)
END SUBROUTINE 

