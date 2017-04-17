SUBROUTINE info_atoms()

  USE m_atoms
  IMPLICIT NONE 
  INTEGER :: ia, isp

  WRITE(*,*)
  WRITE(*,*) 'Atomic coordinates:'
  WRITE(*,*)
  DO ia = 1,Natoms
    isp = atm2species(ia)
    WRITE(*,'(1x,A,3F18.10)') trim(SpeciesSymbols(isp)), AtomicCoords(1:3,ia)
  ENDDO 

END SUBROUTINE 
