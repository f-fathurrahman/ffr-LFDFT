!!>
!!> A driver routine for solving Kohn-Sham equations
!!>
SUBROUTINE do_KS_solve()

  USE m_constants, ONLY : Ry2eV
  USE m_input_vars, ONLY : startingwfc
  USE m_options, ONLY : I_KS_Solve
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  IMPLICIT NONE 
  INTEGER :: ist

  IF( I_KS_SOLVE == 1 ) THEN 

    CALL KS_solve_Emin_pcg( 3.d-5, .FALSE. )
    CALL info_energies()
    CALL calc_evals( Nstates, Focc, evecs, evals )

  ELSEIF( I_KS_SOLVE == 2 ) THEN

    IF( startingwfc == 'random' ) THEN
      ! Initial Rhoe and potentials
      CALL calc_rhoe( Focc, evecs )
      CALL update_potentials()
    ENDIF

    CALL KS_solve_SCF()
    CALL info_energies()

  ENDIF

  WRITE(*,*)
  WRITE(*,*) 'Final eigenvalues (Ha and eV)'
  WRITE(*,*)
  DO ist = 1,Nstates
    WRITE(*,'(1x,I8,2F18.10)') ist, evals(ist), evals(ist)*2.d0*Ry2eV
  ENDDO

  !FIXME Need tidy up
  CALL write_checkpoint()
  CALL write_KS_evecs('KS_evecs.dat')
END SUBROUTINE 
