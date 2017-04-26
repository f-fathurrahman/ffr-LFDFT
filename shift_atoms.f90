SUBROUTINE shift_atoms()
  
  USE m_atoms, ONLY : atpos => AtomicCoords, Natoms
  USE m_LF3d, ONLY : lingrid => LF3d_lingrid
  IMPLICIT NONE 
  INTEGER :: ia

  DO ia = 1,Natoms
    atpos(:,ia) = atpos(:,ia) - lingrid(:,1)
  ENDDO 

END SUBROUTINE 
