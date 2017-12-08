SUBROUTINE read_system(filein)
  USE m_input_vars
  IMPLICIT NONE 
  CHARACTER(*) :: filein

  ! default values
  nat = 0
  ntyp = 0
  A = 1.d0
  B = 1.d0
  C = 1.d0
  nr1 = 1
  nr2 = 1
  nr3 = 1
  ibrav = 8
  _Nstates_extra = 0
  input_dft = ''  ! default value is set in m_xc  ???

  OPEN(unit=IU, file=filein, status='old')
  READ(IU, nml=SYSTEM)
  CLOSE(IU)

  IF( ibrav /= 8 .AND. ibrav /= 1 ) THEN 
    WRITE(*,*) 'ERROR: unsupported ibrav = ', ibrav
    STOP 
  ENDIF 

  !WRITE(*,*) 'ibrav = ', ibrav
  !WRITE(*,*) 'nat = ', nat
  !WRITE(*,*) 'ntyp = ', ntyp
  !WRITE(*,'(1x,A,F18.10)') 'A = ', A
  !WRITE(*,'(1x,A,F18.10)') 'B = ', B
  !WRITE(*,'(1x,A,F18.10)') 'C = ', C
  !WRITE(*,'(1x,A,3I8)') 'nr1, nr2, n3 = ', nr1, nr2, nr3
  
END SUBROUTINE 

