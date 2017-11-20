!!>
!!> This is the main program for \texttt{ffr\_LFDFT}.
!!>
PROGRAM ffr_LFDFT

  IMPLICIT NONE

!!> \begin{itemize}
!!>
!!> \item
!!> Variables for timing stuffs
  INTEGER :: tstart, counts_per_second, tstop

!!> \item
!!> The following subroutine will print some information about the program.
  CALL welcome()

!!> \item
!!> We start timing from here.
  CALL system_clock( tstart, counts_per_second )

!!> \item
!!> This subroutine handles reading input file, setting-up
!!> basis set and grid points, pseudopotential-related variables, etc.
  CALL setup_ffr_LFDFT()

!!> \item Preparation for solving Kohn-Sham equation
  CALL guess_KS_solutions()

!!> \item This soubroutine is the driver for primary subroutines to solve Kohn-Sham
!!> equation.
  CALL do_KS_solve()

!!> \item Free allocated memory
  CALL cleanup_ffr_LFDFT()

!!> \item Stop timing and display elapsed time
  CALL system_clock( tstop )
  WRITE(*,*)
  WRITE(*,'(1x,A,ES18.10,A)') 'Total elapsed time: ', &
           dble(tstop - tstart)/counts_per_second, ' second.'

!!> \item Display goodbye message
  CALL goodbye()

END PROGRAM
!!>\end{itemize}
