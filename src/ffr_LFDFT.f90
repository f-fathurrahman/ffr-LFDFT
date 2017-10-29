!!>
!!> This is the main program for \texttt{ffr\_LFDFT}.
!!>
!!> TODO: The main program should only call several other drivers subroutines.
!!
PROGRAM ffr_LFDFT

  IMPLICIT NONE
  INTEGER :: tstart, counts_per_second, tstop

  CALL welcome()

  CALL system_clock( tstart, counts_per_second )

  CALL setup_ffr_LFDFT()

  CALL guess_KS_solutions()

  CALL do_KS_solve()

  CALL cleanup_ffr_LFDFT()

  CALL system_clock( tstop )

  WRITE(*,*)
  WRITE(*,'(1x,A,ES18.10,A)') 'Total elapsed time: ', &
           dble(tstop - tstart)/counts_per_second, ' second.'

  CALL goodbye()

END PROGRAM




