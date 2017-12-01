!!>
!!> \section{Main program: \texttt{ffr\_LFDFT}}
!!>
!!> This is the main program for \texttt{ffr\_LFDFT}.
!!>
!!> Start from here if you want to learn the call sequence of the main program.
!!> Most of the subroutines called here are wrapper subroutines which call
!!> various core computational subroutines.
!!>
!!> If you want to quickly do a calculation for simple system or want to use directly
!!> the core computational subroutines, please see examples in directory \texttt{\tests}.
!!>
PROGRAM ffr_LFDFT

  IMPLICIT NONE

!!> Variables for timing stuffs
!!> Note that, for the moment the program only does timing for overall subroutines.
!!> It is very desirable, however, to measure timing of various steps in the program.
  INTEGER :: tstart, counts_per_second, tstop

!!>
!!> Here begins the calls of wrapper subroutines.
!!
!!> \begin{itemize}
!!>
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

!!> \item
!!> Preparation for solving Kohn-Sham equation
  CALL guess_KS_solutions()

!!> \item
!!> This soubroutine is the driver for primary subroutines to solve Kohn-Sham
!!> equation.
  CALL do_KS_solve()

!!> \item
!!> Free allocated memory
  CALL cleanup_ffr_LFDFT()

!!> \item
!!> Stop timing and display elapsed time
  CALL system_clock( tstop )
  WRITE(*,*)
  WRITE(*,'(1x,A,ES18.10,A)') 'Total elapsed time: ', &
           dble(tstop - tstart)/counts_per_second, ' second.'

!!> \item
!!> Display goodbye message
  CALL goodbye()

END PROGRAM
!!>\end{itemize}
