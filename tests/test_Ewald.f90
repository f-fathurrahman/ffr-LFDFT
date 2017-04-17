PROGRAM test_Ewald
  
  USE m_atoms, ONLY : Zv => AtomicValences, &
                      Nspecies, Natoms, atm2species, &
                      atpos => AtomicCoords, &
                      strf => StructureFactor

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     Gv => LF3d_Gv

  IMPLICIT NONE 
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)

  !CALL init_atoms_xyz('../structures/NH3.xyz')
  ! Manually set Zv => AtomicValences
  !Zv(1) = 5.d0
  !Zv(2) = 1.d0

  CALL init_atoms_xyz('../structures/H2.xyz')
  ! Manually set Zv => AtomicValences
  Zv(1) = 1.d0

  CALL info_atoms()

  NN = (/ 63, 63, 63 /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)
  CALL init_LF3d_p( NN, AA, BB )
  CALL info_LF3d()

  ALLOCATE( strf(Npoints, Nspecies) )

  CALL calc_strfact( Natoms, atpos, Nspecies, atm2species, Npoints, Gv, strf )

  CALL calc_Ewald()

  DEALLOCATE( strf )
  CALL dealloc_LF3d()
  CALL dealloc_atoms()

END PROGRAM 

