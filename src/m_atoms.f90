MODULE m_atoms

  IMPLICIT NONE 

  INTEGER :: Natoms
  INTEGER :: Nspecies
  REAL(8), ALLOCATABLE :: AtomicCoords(:,:)
  INTEGER, ALLOCATABLE :: atm2species(:)
  CHARACTER(5), ALLOCATABLE :: SpeciesSymbols(:)  ! Nspecies
  REAL(8), ALLOCATABLE :: AtomicValences(:)  ! should be SpeciesValences
  REAL(8), ALLOCATABLE :: AtomicMasses(:)    ! should be SpeciesMasses

  COMPLEX(8), ALLOCATABLE :: StructureFactor(:,:)

END MODULE 

