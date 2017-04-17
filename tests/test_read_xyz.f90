PROGRAM test_read_xyz

  IMPLICIT NONE 

  CALL init_atoms_xyz('../structures/NH3.xyz')
  CALL info_atoms()

  CALL dealloc_atoms()

END PROGRAM 

