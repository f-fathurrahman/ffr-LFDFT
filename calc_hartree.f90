! Calculate Hartree potential and energy
SUBROUTINE calc_hartree()
  USE m_globals, ONLY : Rho, Vhartree, Ehartree, PI, LF
  IMPLICIT NONE

  WRITE(*,'(/1x,A)') 'Calculating Hartree potential'
  WRITE(*,*)         '-----------------------------'

  CALL solve_poisson_cg( -4.d0*PI*Rho, Vhartree)

  Ehartree = 0.5d0*sum(Rho*Vhartree)*LF%LFx%h*LF%LFy%h*LF%LFz%h
  WRITE(*,*) 'Ehartree = ', Ehartree
END SUBROUTINE



