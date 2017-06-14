SUBROUTINE init_states()

  USE m_states
  USE m_atoms, ONLY : Zv => AtomicValences, Natoms, atm2species
  IMPLICIT NONE 
  INTEGER :: ist, isp, ia
 
  Nspin = 1  ! default
  Nstates_extra = 0 ! default

  Nelectrons = 0.d0
  DO ia = 1,Natoms
    isp = atm2species(ia)
    Nelectrons = Nelectrons + Zv(isp)
  ENDDO 
  WRITE(*,*)
  WRITE(*,'(1x,A,F18.5)') 'Number of electrons = ', Nelectrons

  Nstates_occ = ceiling( Nelectrons / 2.d0 )
  Nstates = Nstates_occ + Nstates_extra
  ALLOCATE( Focc(Nstates) )
  WRITE(*,'(1x,A,I8)') 'Nstates = ', Nstates

  Focc(:) = 0.d0
  IF( mod( int(Nelectrons), 2 ) == 0 ) THEN
    DO ist = 1,Nstates_occ
      Focc(ist) = 2.d0
    ENDDO 
  ELSE 
    WRITE(*,*) 'WARNING: Odd number of electrons = ', int(Nelectrons)
    DO ist = 1,Nstates_occ-1
      Focc(ist) = 2.d0
    ENDDO 
    Focc(Nstates_occ) = 1.d0
  ENDIF

  WRITE(*,*)
  WRITE(*,*) 'Initial occupation numbers'
  DO ist = 1,Nstates
    WRITE(*,'(1x,I8,F18.5)') ist, Focc(ist)
  ENDDO 

END SUBROUTINE 
