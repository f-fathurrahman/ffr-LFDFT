SUBROUTINE read_electrons(filein)
  USE m_input_vars
  IMPLICIT NONE 
  CHARACTER(*) :: filein
  
  ! Initial/default values. Most of them are not sensible inputs.
  ! This is used as mechanism to check whether user inputted the values or not.
  KS_Solve = 'undefined'
  cg_beta = 'undefined'
  electron_maxstep = -1
  mixing_beta = -1.d0
  mixing_mode = 'default'
  diagonalization = 'undefined'
  conv_thr = -1.d0
  startingwfc = 'atomic'

  OPEN(unit=IU,file=filein,status='old')
  READ(IU,nml=ELECTRONS)
  CLOSE(IU)

END SUBROUTINE 

