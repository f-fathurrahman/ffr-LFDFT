SUBROUTINE read_electrons(filein)
  USE m_input_vars
  IMPLICIT NONE 
  CHARACTER(*) :: filein

  OPEN(unit=IU,file=filein,status='old')
  READ(IU,nml=ELECTRONS)
  CLOSE(IU)

  !WRITE(*,*) KS_Solve_Method
  !WRITE(*,*) icg_beta
  !WRITE(*,*) electron_maxstep
  !WRITE(*,*) mixing_beta
  !WRITE(*,*) diagonalization

END SUBROUTINE 

