SUBROUTINE read_input(filein)
  IMPLICIT NONE 
  CHARACTER(*) :: filein

  WRITE(*,*)
  WRITE(*,'(1x,2A)') 'Reading input file: ', trim(filein)

  CALL read_control(filein)
  CALL read_system(filein)
  CALL read_electrons(filein)
  CALL read_atomic_species(filein)
  CALL read_atomic_positions(filein)

  WRITE(*,*)
  WRITE(*,*) 'Finished reading input'

END SUBROUTINE 
