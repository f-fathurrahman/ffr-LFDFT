!!>
!!> \section{Subroutine \texttt{guess\_KS\_solutions}}
!!>
!!> Generate guess solutions (density and orbitals) for \texttt{KS\_solve\_XXX}
!!> subroutines.
!!>
!!> Note that the global arrays \texttt{KS\_evals} and \texttt{KS\_evecs} are
!!> allocated in this subroutine.
!!>
SUBROUTINE guess_KS_solutions()

  USE m_input_vars, ONLY : startingwfc
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_states, ONLY : Nstates, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  !
  ! Guess density
  !
  ! FIXME: Need gen_guess_rho_gaussian for isolated system ?
  !        thus bypassing calculation of structure factors ?
  IF( startingwfc /= 'random' ) THEN
    CALL gen_guess_rho_gaussian()
  ENDIF

  ! Manually allocate KS eigenvectors and eigenvalues
  ALLOCATE( evecs(Npoints,Nstates), evals(Nstates) )

  CALL gen_random_evecs()  ! also needed for initial diagonalization routine

  IF( startingwfc /= 'random' ) THEN
    ! This will call diagonalization routine
    CALL gen_gaussian_evecs()
  ELSE
    WRITE(*,*)
    WRITE(*,*) 'Using random starting wavefunction'
  ENDIF
END SUBROUTINE
