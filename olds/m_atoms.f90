! eFeFeR, September 2014
! eFeFeR, June 2012

!--------------------!
MODULE m_atoms
!--------------------!

IMPLICIT NONE

TYPE atoms_t
  ! Lattice constant
  REAL(8) :: alat
  !
  ! Unit cell volume
  !
  REAL(8) :: cell_volume
  !
  ! Lattice vectors
  !
  REAL(8) :: A1(3), A2(3), A3(3)
  !
  ! Reciprocal lattice vectors
  !
  REAL(8) :: B1(3), B2(3), B3(3)
  !
  ! Number of species
  !
  INTEGER :: Nspecies
  !
  ! Number of total atoms
  !
  INTEGER :: Natoms
  !
  ! Valence charge of the species
  !
  REAL(8), ALLOCATABLE :: Zv(:) ! ZV(1:NSPECIES)
  !
  ! Atomic positions, in Bohr
  !
  REAL(8), ALLOCATABLE :: positions(:,:) ! positions(1:3,1:NATOMS)
  !
  ! These variables will be allocated after reading input files
  !
  INTEGER, ALLOCATABLE :: atmSpec(:) ! atmSpec(NATOMS)  FIXME: Not used?
  INTEGER, ALLOCATABLE :: atmToSpecies(:) ! atmToSpecies(NATOMS)
  CHARACTER(20), ALLOCATABLE :: atmSymb(:) ! atmSymb(NATOMS)
  CHARACTER(20), ALLOCATABLE :: speciesSymb(:) ! speciesSymb(NSPECIES)
  !
END TYPE


CONTAINS


! Initialize a PERIODIC structure.
! Atomic positions are read from fil_xyz
! Lattice vectors are given  by A1, A2, and A3
SUBROUTINE init_atoms_xyz_per(atoms, fil_xyz, A1, A2, A3)
  USE m_constants, ONLY: ANG2BOHR
  IMPLICIT NONE
  !
  ! Arguments
  !
  CHARACTER(*) :: fil_xyz
  REAL(8) :: A1(3), A2(3), A3(3)
  !
  ! Output
  !
  TYPE(atoms_t) :: atoms
  !
  ! Functions
  !
  REAL(8) :: det3
  !
  
  CALL init_atoms_xyz(atoms, fil_xyz)

  !
  ! Unit cell related setup
  !
  ! Set lattice vectors
  atoms%A1 = A1
  atoms%A2 = A2
  atoms%A3 = A3
  !
  atoms%cell_volume = det3(A1,A2,A3)
  !
  CALL calc_reciplat( atoms )

END SUBROUTINE


! Initialize a non-periodic structure.
! Atomic positions are read from fil_xyz
SUBROUTINE init_atoms_xyz(atoms, fil_xyz)
  USE m_constants, ONLY: ANG2BOHR
  IMPLICIT NONE
  !
  ! Arguments
  !
  CHARACTER(*) :: fil_xyz
  !
  ! Output
  !
  TYPE(atoms_t) :: atoms
  !
  ! Local
  !
  INTEGER :: unitxyz
  INTEGER :: ios, ia, k1, k2, isp, idx1
  INTEGER :: natoms, nspecies
  !
  OPEN( UNIT=unitxyz, FILE=fil_xyz, ACTION='read', STATUS='old', &
    FORM='formatted', IOSTAT=ios )
  !
  ! Read
  !
  READ( unitxyz, * ) natoms
  !WRITE(*,*) 'natoms = ', natoms
  atoms%natoms = natoms

  ALLOCATE( atoms%atmSymb( natoms ) )
  ALLOCATE( atoms%positions( 3, natoms ) )

  DO ia=1,natoms
    READ( unitxyz, * ) atoms%atmSymb(ia), atoms%positions(1,ia), &
      atoms%positions(2,ia), atoms%positions(3,ia)
    !WRITE(*,'(1x,A,3G18.10)') trim(atoms%atmSymb(ia)), atoms%positions(1,ia), &
    !  atoms%positions(2,ia), atoms%positions(3,ia)
  ENDDO
  !
  ! convert from angstrom to bohr
  atoms%positions(:,:) = atoms%positions(:,:)*ANG2BOHR

  ! Determine number of species
  Nspecies = 0
  DO ia=1,NATOMS
    k2 = 0
    DO k1=1,ia-1
      IF( atoms%atmSymb(k1) == atoms%atmSymb(ia) ) k2=1
    ENDDO
    ! find different
    IF(k2==0) THEN
      !WRITE(*,*) 'Find new species'
      NSPECIES = NSPECIES + 1
    ENDIF
  ENDDO
  atoms%nspecies = Nspecies
  !WRITE(*,*) 'nspecies = ', atoms%nspecies

  ALLOCATE( atoms%speciesSymb(Nspecies) )
  idx1 = 0
  DO ia=1,Natoms
    k2 = 0
    DO k1=1,ia-1
      IF( atoms%atmSymb(k1) == atoms%atmSymb(ia) ) k2=1
    ENDDO
    ! find different
    IF(k2==0) THEN
      idx1 = idx1 + 1
      atoms%speciesSymb(idx1) = atoms%atmSymb(ia)
    ENDIF
  ENDDO
  !
  !WRITE(*,*) 'Species:'
  !DO isp=1,Nspecies
  !  WRITE(*,*) trim(atoms%speciesSymb(isp))
  !ENDDO

  ! Mapping of atoms to species index
  ALLOCATE( atoms%atmToSpecies(Natoms) )
  DO ia=1,Natoms
    DO isp=1,Nspecies
      IF( atoms%atmSymb(ia) == atoms%speciesSymb(isp) ) THEN
        atoms%atmToSpecies(ia) = isp
      ENDIF
    ENDDO 
  ENDDO

  CLOSE( unitxyz )
END SUBROUTINE


SUBROUTINE calc_reciplat(at1)
  IMPLICIT NONE
  ! Argument
  TYPE(atoms_t) :: at1
  ! Local
  REAL(8) :: SCAL
  REAL(8) :: A1(3), A2(3), A3(3)
  REAL(8) :: B1(3), B2(3), B3(3)

  A1 = at1%A1
  A2 = at1%A2
  A3 = at1%A3
  
  SCAL = A1(1)/at1%cell_volume
  
  B1(1) = ( A2(2)*A3(3) - A2(3)*A3(2) )*SCAL
  B1(2) = ( A1(3)*A3(2) - A1(2)*A3(3) )*SCAL
  B1(3) = ( A1(2)*A2(3) - A1(3)*A2(2) )*SCAL
  
  B2(1) = ( A2(3)*A3(1) - A2(1)*A3(3) )*SCAL
  B2(2) = ( A1(1)*A3(3) - A1(3)*A3(1) )*SCAL
  B2(3) = ( A1(3)*A2(1) - A1(1)*A2(3) )*SCAL

  B3(1) = ( A2(1)*A3(2) - A2(2)*A3(1) )*SCAL
  B3(2) = ( A1(2)*A3(1) - A1(1)*A3(2) )*SCAL
  B3(3) = ( A1(1)*A2(2) - A1(2)*A2(1) )*SCAL

  at1%B1 = B1
  at1%B2 = B2
  at1%B3 = B3

END SUBROUTINE

!TODO
SUBROUTINE set_positions()
  IMPLICIT NONE
END SUBROUTINE

!TODO
SUBROUTINE set_pbc()
  IMPLICIT NONE
END SUBROUTINE

!
! Display information about atoms_t object
!
SUBROUTINE print_atoms_per(at1)
  IMPLICIT NONE
  !
  TYPE(atoms_t) :: at1
  INTEGER :: ia, ii

  WRITE(*,*)
  WRITE(*,*) 'Information about at1 object'
  !
  WRITE(*,*) 'NATOMS   = ', at1%natoms
  WRITE(*,*) 'NSPECIES = ', at1%nspecies
  !
  WRITE(*,*) 'Atomic positions:'
  DO ia=1,at1%natoms
    WRITE(*,'(1x,A,3F18.10)') trim(at1%atmSymb(ia)), at1%positions(1,ia), &
      at1%positions(2,ia), at1%positions(3,ia)
  ENDDO
  !
  WRITE(*,*) 'Unit cell vectors:'
  DO ii=1,3
    WRITE(*,'(3F18.10)') at1%A1(ii), at1%A2(ii), at1%A3(ii)
  ENDDO
  WRITE(*,*) 'Unit cell volume = ', at1%cell_volume
  !
  WRITE(*,*) 'Reciprocal unit cell vectors:'
  DO ii=1,3
    WRITE(*,'(3F18.10)') at1%B1(ii), at1%B2(ii), at1%B3(ii)
  ENDDO
 
END SUBROUTINE


!
! Display information about atoms_t object
!
SUBROUTINE print_atoms(at1)
  IMPLICIT NONE
  !
  TYPE(atoms_t) :: at1
  INTEGER :: ia

  WRITE(*,*)
  WRITE(*,*) 'Information about at1 object'
  !
  WRITE(*,*) 'NATOMS   = ', at1%Natoms
  WRITE(*,*) 'NSPECIES = ', at1%Nspecies
  !
  WRITE(*,*) 'Atomic positions: (in Bohr)'
  DO ia=1,at1%natoms
    WRITE(*,'(1x,A,3F18.10)') trim(at1%atmSymb(ia)), at1%positions(1,ia), &
      at1%positions(2,ia), at1%positions(3,ia)
  ENDDO 
END SUBROUTINE

END MODULE

