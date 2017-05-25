SUBROUTINE setup_atoms()
  USE m_constants, ONLY : ANG2BOHR
  USE m_input_vars, ONLY : ntyp, nat, species, in_pos, masses, in_atmsymb
  USE m_atoms
  IMPLICIT NONE 
  INTEGER :: ia, isp

  Natoms = nat
  Nspecies = ntyp

  ALLOCATE( AtomicCoords(3,Natoms) )
  ALLOCATE( SpeciesSymbols(Nspecies) )
  ALLOCATE( AtomicMasses(Nspecies) )

  AtomicCoords(:,:) = in_pos(:,:)*ANG2BOHR
  SpeciesSymbols(:) = species(:)
  AtomicMasses(:)   = masses(:)

  ALLOCATE( atm2species(Natoms) )
  DO ia=1,Natoms
    DO isp=1,Nspecies
      IF( in_atmsymb(ia) == SpeciesSymbols(isp) ) THEN
        atm2Species(ia) = isp
      ENDIF
    ENDDO 
  ENDDO

  ALLOCATE( AtomicValences(Nspecies) ) ! will set when reading pseudopotentials
  AtomicValences(:) = 0.d0

  ! Free memory
  DEALLOCATE( species )
  DEALLOCATE( masses )
  DEALLOCATE( in_pos )
  DEALLOCATE( in_atmsymb )

END SUBROUTINE 

