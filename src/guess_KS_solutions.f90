!!>
!!> \section{Subroutine \texttt{guess\_KS\_solutions}}
!!>
!!> Generate guess solutions (density and orbitals) for \texttt{KS\_solve\_XXX}
!!> subroutines.
!!>
!!> Note that the global arrays \texttt{KS\_evals} and \texttt{KS\_evecs} are
!!> allocated in this subroutine.
!!>
!!> The algorithm used in this subroutine is not sophisticated. I have used ABINIT
!!> routine \texttt{atmlength} to generate gaussian charge density.
!!>
!!> For direct minimization (using \texttt{KS\_solve\_Emin\_XXX}),
!!> we need to generate initial wavefunction, so an
!!> iterative diagonalization step must be performed. This might cause significant
!!> additional time before the actual direct minimization step starts.
!!>
!!> For SCF, we can simply generate random initial wavefunction which will be used for
!!> iterative diagonalization step.
!!>
SUBROUTINE guess_KS_solutions()

!!> Here are the imported variables:
  USE m_input_vars, ONLY : startingwfc

!!> \begin{itemize}  
!!> \item
!!> Generate gaussian charge density.
!!>
  ! FIXME: Need gen_guess_rho_gaussian for isolated system ?
  !        thus bypassing calculation of structure factors ?
  IF( startingwfc /= 'random' ) THEN
    CALL gen_guess_rho_gaussian()
  ENDIF

!!> \item
!!> Generate random initial wavefunction.
!!> This step is needed for both SCF and direct minimization method.

  CALL gen_random_evecs()  ! also needed for initial diagonalization routine

!!> \item
!!> If the input variable \texttt{startingwfc} is not set to \texttt{'random'}
!!> then we need to do one iterative diagonalization step.
!!> If it is set to random then we don't need to do iterative diagonalization
!!> step as random wavefunction is suffice and subroutine will immediately return.
  IF( startingwfc /= 'random' ) THEN
    ! This will call diagonalization routine
    CALL gen_gaussian_evecs()
  ELSE
    WRITE(*,*)
    WRITE(*,*) 'Using random starting wavefunction'
  ENDIF
END SUBROUTINE
!!> \end{itemize}