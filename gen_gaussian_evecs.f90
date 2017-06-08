SUBROUTINE gen_gaussian_evecs()
  IMPLICIT NONE 

  CALL update_potentials()
  CALL Sch_solve_diag()

END SUBROUTINE 

