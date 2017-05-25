SUBROUTINE read_input(filein)
  IMPLICIT NONE 
  CHARACTER(*) :: filein

  CALL read_control(filein)
  CALL read_system(filein)
  CALL read_electrons(filein)
  CALL read_atomic_species(filein)
  CALL read_atomic_positions(filein)

END SUBROUTINE 
