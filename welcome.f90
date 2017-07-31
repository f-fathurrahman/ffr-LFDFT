SUBROUTINE welcome()
  IMPLICIT NONE

  WRITE(*,*)
  WRITE(*,*) "  __  __             _      ______ _____  ______ _______ "
  WRITE(*,*) " / _|/ _|           | |    |  ____|  __ \\|  ____|__   __|"
  WRITE(*,*) "| |_| |_ _ __ ______| |    | |__  | |  | | |__     | |   "
  WRITE(*,*) "|  _|  _| '__|______| |    |  __| | |  | |  __|    | |   "
  WRITE(*,*) "| | | | | |         | |____| |    | |__| | |       | |   "
  WRITE(*,*) "|_| |_| |_|         |______|_|    |_____/|_|       |_|   "

  CALL timestamp('Program started at')

END SUBROUTINE
