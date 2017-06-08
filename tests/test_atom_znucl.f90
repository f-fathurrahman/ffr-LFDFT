PROGRAM test_atom_znucl

  IMPLICIT NONE 
  REAL(8) :: atom_znucl

  WRITE(*,*) 'H  ', atom_znucl('H')
  WRITE(*,*) 'He ', atom_znucl('He')
  WRITE(*,*) 'Li ', atom_znucl('Li')
  WRITE(*,*) 'Be ', atom_znucl('Be')
  WRITE(*,*) 'B  ', atom_znucl('B')
  WRITE(*,*) 'C  ', atom_znucl('C')
  WRITE(*,*) 'N  ', atom_znucl('N')

END PROGRAM 

