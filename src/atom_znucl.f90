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
  CASE('Li', 'Li_sc')
    znucl = 3.d0
  CASE('Be', 'Be_sc')
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
  CASE('K', 'K_sc')
    znucl = 19.d0
  CASE('Ca', 'Ca_sc')
    znucl = 20.d0
  CASE('Sc', 'Sc_sc')
    znucl = 21.d0
  CASE('Ti', 'Ti_sc')
    znucl = 22.d0
  CASE('V', 'V_sc')
    znucl = 23.d0
  CASE('Cr', 'Cr_sc')
    znucl = 24.d0
  CASE('Mn', 'Mn_sc')
    znucl = 25.d0
  CASE('Fe', 'Fe_sc')
    znucl = 26.d0
  CASE('Co', 'Co_sc')
    znucl = 27.d0
  CASE('Ni', 'Ni_sc')
    znucl = 28.d0
  CASE('Cu', 'Cu_sc')
    znucl = 29.d0
  CASE('Zn', 'Zn_sc')
    znucl = 30.d0
  CASE('Ga', 'Ga_sc')
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
  CASE('Rb', 'Rb_sc')
    znucl = 37.d0
  CASE('Sr', 'Sr_sc')
    znucl = 38.d0
  CASE('Y', 'Y_sc')
    znucl = 39.d0
  CASE('Zr', 'Zr_sc')
    znucl = 40.d0
  CASE('Nb', 'Nb_sc')
    znucl = 41.d0
  CASE('Mo', 'Mo_sc')
    znucl = 42.d0
  CASE('Tc', 'Tc_sc')
    znucl = 43.d0
  CASE('Ru', 'Ru_sc')
    znucl = 44.d0
  CASE('Rh', 'Rh_sc')
    znucl = 45.d0
  CASE('Pd', 'Pd_sc')
    znucl = 46.d0
  CASE('Ag', 'Ag_sc')
    znucl = 47.d0
  CASE('Cd', 'Cd_sc')
    znucl = 48.d0
  CASE('In', 'In_sc')
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
  CASE('Cs', 'Cs_sc')
    znucl = 55.d0
  CASE('Ba', 'Ba_sc')
    znucl = 56.d0
  !
  ! Lanthanide and actinide are not complete
  !
  CASE('Ce', 'Ce_sc')
    znucl = 58.d0
  CASE('Hf', 'Hf_sc')
    znucl = 72.d0
  CASE('Ta', 'Ta_sc')
    znucl = 73.d0
  CASE('W', 'W_sc')
    znucl = 74.d0
  CASE('Re', 'Re_sc')
    znucl = 75.d0
  CASE('Os', 'Os_sc')
    znucl = 76.d0
  CASE('Ir', 'Ir_sc')
    znucl = 77.d0
  CASE('Pt', 'Pt_sc')
    znucl = 78.d0
  CASE('Au', 'Au_sc')
    znucl = 79.d0
  CASE('Hg', 'Hg_sc')
    znucl = 80.d0
  CASE('Tl', 'Tl_sc')
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

