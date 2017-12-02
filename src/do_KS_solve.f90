!!>
!!> \section{Subroutine \texttt{do\_KS\_solve}}
!!>
!!> A driver routine for solving Kohn-Sham equations.
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

!!>
!!> We are using direct minimization (this is the default):
!!>
  IF( I_KS_SOLVE == 1 ) THEN 
!!> Call the main computational routine for direct minimization.
!!> Note that for the moment the parameter $\alpha_t$ is hardcoded to $3\times10^{-5}$
    CALL KS_solve_Emin_pcg( 3.d-5, .FALSE. )
!!> Display total energy components
    CALL info_energies()
!!> We need to calculate eigenvalues explicitly
    CALL calc_evals( Nstates, Focc, evecs, evals )

!!>
!!> We are using the conventional SCF with electron density mixing
!!>
  ELSEIF( I_KS_SOLVE == 2 ) THEN
!!> We need to calculate electron density and update the local potentials
!!> if starting from random wavefunction
    IF( startingwfc == 'random' ) THEN
      ! Initial Rhoe and potentials
      CALL calc_rhoe( Focc, evecs )
      ! FIXME: need this ?
      CALL update_potentials()
    ENDIF
!!> Call the main computational routine for SCF
    CALL KS_solve_SCF()
    CALL info_energies()

  ENDIF

!!> Eigenvalues are displayed here
  WRITE(*,*)
  WRITE(*,*) 'Final eigenvalues (Ha and eV)'
  WRITE(*,*)
  DO ist = 1,Nstates
    WRITE(*,'(1x,I8,2F18.10)') ist, evals(ist), evals(ist)*2.d0*Ry2eV
  ENDDO

!!> Write restart data
  !FIXME Need tidy up
  CALL write_checkpoint()
  CALL write_KS_evecs('KS_evecs.dat')
END SUBROUTINE 
