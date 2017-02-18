SUBROUTINE print_energies()
  USE m_globals, ONLY : Etot, Ekin, Epot, Ehartree, Exc
  IMPLICIT NONE
  WRITE(*,*)
  WRITE(*,'(1x,A,F18.10)') 'Etot     = ', Etot
  WRITE(*,'(1x,A,F18.10)') '-----------------------------'
  WRITE(*,'(1x,A,F18.10)') 'Ekin     = ', Ekin
  WRITE(*,'(1x,A,F18.10)') 'Epot     = ', Epot
  WRITE(*,'(1x,A,F18.10)') 'Ehartree = ', Ehartree
  WRITE(*,'(1x,A,F18.10)') 'Exc      = ', Exc
END SUBROUTINE


SUBROUTINE print_evals( Nstates, evals )
  IMPLICIT NONE
  !
  INTEGER :: Nstates
  REAL(8) :: evals(Nstates)
  !
  INTEGER :: is

  WRITE(*,*)
  WRITE(*,*) 'Eigenvalues:'
  WRITE(*,*) '-----------------------'
  DO is = 1, Nstates
    WRITE(*,'(1x,I5,F18.10)') is, evals(is)
  ENDDO
END SUBROUTINE
