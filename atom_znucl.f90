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
  CASE('O')
    znucl = 8.d0
  CASE('F')
    znucl = 9.d0
  CASE('Ne')
    znucl = 10.d0
  CASE('Na')
    znucl = 11.d0
  CASE('Mg')
    znucl = 12.d0
  CASE('Al')
    znucl = 13.d0
  CASE('Si')
    znucl = 14.d0
  CASE('P')
    znucl = 15.d0
  CASE('S')
    znucl = 16.d0
  CASE('Cl')
    znucl = 17.d0
  CASE('Ar')
    znucl = 18.d0
  CASE('K')
    znucl = 19.d0
  CASE('Ca')
    znucl = 20.d0
  CASE('Sc')
    znucl = 21.d0
  CASE('Ti')
    znucl = 22.d0
  CASE('V')
    znucl = 23.d0
  CASE('Cr')
    znucl = 24.d0
  CASE('Mn')
    znucl = 25.d0
  CASE('Fe')
    znucl = 26.d0
  CASE('Co')
    znucl = 27.d0
  CASE('Ni')
    znucl = 28.d0
  CASE('Cu')
    znucl = 29.d0
  CASE('Zn')
    znucl = 30.d0
  CASE('Ga')
    znucl = 31.d0
  CASE('Ge')
    znucl = 32.d0
  CASE('As')
    znucl = 33.d0
  CASE('Se')
    znucl = 34.d0
  CASE('Br')
    znucl = 35.d0
  CASE('Kr')
    znucl = 36.d0
  CASE('Rb')
    znucl = 37.d0
  CASE('Sr')
    znucl = 38.d0
  CASE('Y')
    znucl = 39.d0
  CASE('Zr')
    znucl = 40.d0
  CASE('Nb')
    znucl = 41.d0
  CASE('Mo')
    znucl = 42.d0
  CASE('Tc')
    znucl = 43.d0
  CASE('Ru')
    znucl = 44.d0
  CASE('Rh')
    znucl = 45.d0
  CASE('Pd')
    znucl = 46.d0
  CASE('Ag')
    znucl = 47.d0
  CASE('Cd')
    znucl = 48.d0
  CASE('In')
    znucl = 49.d0
  CASE('Sn')
    znucl = 50.d0
  CASE('Sb')
    znucl = 51.d0
  CASE('Te')
    znucl = 52.d0
  CASE('I')
    znucl = 53.d0
  CASE('Xe')
    znucl = 54.d0
  CASE('Cs')
    znucl = 55.d0
  CASE('Ba')
    znucl = 56.d0
  !
  ! Lanthanide skipped for the moment ...
  !
  CASE('Hf')
    znucl = 72.d0
  CASE('Ta')
    znucl = 73.d0
  CASE('W')
    znucl = 74.d0
  CASE('Re')
    znucl = 75.d0
  CASE('Os')
    znucl = 76.d0
  CASE('Ir')
    znucl = 77.d0
  CASE('Pt')
    znucl = 78.d0
  CASE('Au')
    znucl = 79.d0
  CASE('Hg')
    znucl = 80.d0
  CASE('Tl')
    znucl = 81.d0
  CASE('Pb')
    znucl = 82.d0
  CASE('Bi')
    znucl = 83.d0
  CASE('Po')
    znucl = 84.d0
  CASE('At')
    znucl = 85.d0
  CASE('Rn')
    znucl = 86.d0
  CASE DEFAULT
    WRITE(*,*) 'ERROR: Unknown/uncoded atomic symbols: ', trim(symb)
    STOP 
  END SELECT 
END FUNCTION 

