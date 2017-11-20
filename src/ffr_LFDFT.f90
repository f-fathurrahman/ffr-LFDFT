!!>
!!> This is the main program for \texttt{ffr\_LFDFT}.
!!>
!!
PROGRAM ffr_LFDFT

  IMPLICIT NONE

!!> Variables for timing stuffs
  INTEGER :: tstart, counts_per_second, tstop

!!> The following subroutine will print some information about the program.
  CALL welcome()

!!> We start timing from here.
  CALL system_clock( tstart, counts_per_second )

!!> Here is the first step. This subroutine handles reading input file, setting-up
!!> basis set and grid points, pseudopotential-related variables, etc.
  CALL setup_ffr_LFDFT()

!!> Preparation for solving Kohn-Sham equation
  CALL guess_KS_solutions()

!!> This soubroutine is the driver for primary subroutines to solve Kohn-Sham
!!> equation.
  CALL do_KS_solve()

!!> Free allocated memory
  CALL cleanup_ffr_LFDFT()

!!> Stop timing and display elapsed time
  CALL system_clock( tstop )
  WRITE(*,*)
  WRITE(*,'(1x,A,ES18.10,A)') 'Total elapsed time: ', &
           dble(tstop - tstart)/counts_per_second, ' second.'

!!> Display goodbye message
  CALL goodbye()

END PROGRAM
