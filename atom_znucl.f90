FUNCTION atom_znucl( symb ) RESULT( znucl )
  IMPLICIT NONE 
  CHARACTER(*) :: symb
  REAL(8) :: znucl

  znucl = 1.d0  !? to supress warning in g95

  SELECT CASE(trim(symb))
  CASE('H')
    znucl = 1.d0
  CASE('He')
    znucl = 2.d0
  CASE('Li')
    znucl = 3.d0
  CASE('Be')
    znucl = 4.d0
  CASE('B')
    znucl = 5.d0
  CASE('C')
    znucl = 6.d0
  CASE('N')
    znucl = 7.d0
  CASE DEFAULT
    znucl = 0.d0
  END SELECT 
END FUNCTION 

