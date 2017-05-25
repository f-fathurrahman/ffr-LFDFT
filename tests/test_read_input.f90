PROGRAM test_read_input

  IMPLICIT NONE 

  CALL read_control('INPUT')
  CALL read_system('INPUT')
  CALL read_electrons('INPUT')
  CALL read_atomic_species('INPUT')
  CALL read_atomic_positions('INPUT')

  WRITE(*,*) 'End of test_read_input'

END PROGRAM 

