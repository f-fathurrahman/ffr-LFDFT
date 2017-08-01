SUBROUTINE gen_gaussian_evecs()
  IMPLICIT NONE 

  CALL update_potentials()

  WRITE(*,*)
  WRITE(*,*) 'Generating starting wavefunctions'
  WRITE(*,*) '---------------------------------'

  CALL Sch_solve_diag()

END SUBROUTINE 

