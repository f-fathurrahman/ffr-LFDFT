SUBROUTINE welcome()
  IMPLICIT NONE

  ! Generated using figlet
  WRITE(*,*)
  WRITE(*,*) "  __  __             _      ______ _____  ______ _______ "
#ifdef __PGI
  WRITE(*,*) " / _|/ _|           | |    |  ____|  __ \\|  ____|__   __|"
#else
  WRITE(*,*) " / _|/ _|           | |    |  ____|  __ \|  ____|__   __|"
#endif
  WRITE(*,*) "| |_| |_ _ __ ______| |    | |__  | |  | | |__     | |   "
  WRITE(*,*) "|  _|  _| '__|______| |    |  __| | |  | |  __|    | |   "
  WRITE(*,*) "| | | | | |         | |____| |    | |__| | |       | |   "
  WRITE(*,*) "|_| |_| |_|         |______|_|    |_____/|_|       |_|   "

  CALL timestamp('Program started at')

  CALL compile_info()

END SUBROUTINE
