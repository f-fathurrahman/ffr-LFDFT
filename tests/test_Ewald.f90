PROGRAM test_Ewald
  
  USE m_PsPot, ONLY : PsPot_Dir
  USE m_atoms, ONLY : Zv => AtomicValences, &
                      Nspecies, Natoms, atm2species, &
                      atpos => AtomicCoords, &
                      strf => StructureFactor
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     Gv => LF3d_Gv

  IMPLICIT NONE 
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  INTEGER :: Narg, N_in
  INTEGER :: iargc
  CHARACTER(64) :: filexyz, arg_N
  !
  REAL(8) :: ewald

  Narg = iargc()
  IF( Narg /= 2 ) THEN 
    WRITE(*,*) 'ERROR: exactly two arguments must be given:'
    WRITE(*,*) '       N and path to structure file'
    STOP 
  ENDIF 

  CALL getarg( 1, arg_N )
  READ(arg_N, *) N_in

  CALL getarg( 2, filexyz )
  
  CALL init_atoms_xyz(filexyz)

  ! Override PsPot_Dir
  PsPot_Dir = '../HGH/'
  CALL init_PsPot()

  NN = (/ N_in, N_in, N_in /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)
  CALL init_LF3d_p( NN, AA, BB )

  CALL info_LF3d()
  CALL info_atoms()

  CALL init_strfact_shifted()

  CALL calc_Ewald()

  WRITE(*,*)
  WRITE(*,*) 'ewald = ', ewald()

  CALL dealloc_LF3d()
  CALL dealloc_atoms()
  CALL dealloc_LF3d()
  CALL dealloc_PsPot()

END PROGRAM 

