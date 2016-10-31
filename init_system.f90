! ffr

! Initialization of molecular structures, pseudopotentials, etc
SUBROUTINE init_system()
  USE m_atoms, ONLY : atoms_t, init_atoms_xyz, print_atoms
  USE m_globals, ONLY : ATOMS, PSPOTS, Nelec, Nstate, Focc
  USE m_ps_hgh, ONLY : hgh_init, hgh_info
  IMPLICIT NONE
  CHARACTER(len=256) :: filname
  INTEGER :: is, ia

  filname = 'Pd.xyz'

  WRITE(*,'(/,1x,2A)') 'Initializing atomic positions from file: ', trim(filname)
  CALL init_atoms_xyz( ATOMS, filname ) 
  CALL print_atoms( ATOMS )

  ALLOCATE( PSPOTS( ATOMS%Nspecies ) )
  ALLOCATE( ATOMS%Zv(ATOMS%Nspecies) )
  !
  WRITE(*,'(/,1x,A)') 'Reading pseudopotentials:'
  DO is = 1, ATOMS%Nspecies
    ! filename or path for pseudopotential file
    ! currently we only support HGH pseudopotential
    filname = 'HGH/'//trim(ATOMS%speciesSymb(is))//'.hgh'
    !WRITE(*,*) trim(filname)
    CALL hgh_init( PSPOTS(is), filname )
    !
    !CALL hgh_info( PSPOTS(is) )

    ATOMS%Zv(is) = PSPOTS(is)%z_val
    WRITE(*,'(1x,A,A3,F10.5)') 'Species, valency: ', trim(ATOMS%speciesSymb(is)), ATOMS%Zv(is)
  ENDDO

  Nelec = 0.d0
  DO ia = 1, ATOMS%Natoms
    Nelec = Nelec + ATOMS%Zv( ATOMS%atmToSpecies(ia) )
  ENDDO
  WRITE(*,'(1x,A,F15.5)') 'Number of electrons = ', Nelec
  IF( mod(int(Nelec),2) /= 0 ) THEN
    WRITE(*,*) 'WARNING: Odd number of electrons'
    Nstate = int(Nelec)/2 + 1
    ALLOCATE( Focc(Nstate) )
    Focc(:) = 2.d0
    Focc(Nstate) = 1.d0
  ELSE
    Nstate = int(Nelec)/2
    ALLOCATE( Focc(Nstate) )
    Focc(:) = 2.d0
  ENDIF
  WRITE(*,'(1x,A,I9)') 'Number of states    = ', Nstate
  WRITE(*,*) 'Occupation numbers'
  WRITE(*,*) 'focc = ', focc


END SUBROUTINE

