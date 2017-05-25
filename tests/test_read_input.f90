PROGRAM test_read_input

  IMPLICIT NONE 

  CALL read_control('INPUT')
  CALL read_system('INPUT')
  CALL read_atomic_species('INPUT')
  CALL read_atomic_positions('INPUT')

  WRITE(*,*) 'Pass here ...'

END PROGRAM 

