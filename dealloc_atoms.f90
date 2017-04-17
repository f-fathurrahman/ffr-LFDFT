SUBROUTINE dealloc_atoms()

  USE m_atoms
  IMPLICIT NONE 

  IF( allocated( AtomicCoords ) ) DEALLOCATE( AtomicCoords )
  IF( allocated( SpeciesSymbols) ) DEALLOCATE( SpeciesSymbols )
  IF( allocated( atm2species ) ) DEALLOCATE( atm2species )

END SUBROUTINE 
