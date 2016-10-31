SUBROUTINE read_arguments()
  USE m_globals, ONLY : LF_type, N, A, B, h, Solution_Method
  IMPLICIT NONE
  CHARACTER(32) :: buffer

  CALL getarg(1,LF_type)
  WRITE(*,*) 'LF_type = ', LF_type

  CALL getarg(2,buffer)
  WRITE(*,*) 'buffer = ', buffer
  READ(buffer,*) N
  WRITE(*,*) 'N = ', N

  IF( LF_type == 'sinc' ) THEN
    CALL getarg(3,buffer)
    READ(buffer,*) h
    WRITE(*,*) 'h = ', h
  ELSEIF( LF_type == 'box' .OR. LF_type == 'per' ) THEN
    CALL getarg(3,buffer)
    READ(buffer,*) B
    A = -B
    !WRITE(*,*) 'A, B = ', A, B
  ELSE
    WRITE(*,*) 'Unrecognized LF_type = ', LF_type
    STOP
  ENDIF

  CALL getarg(4,Solution_Method)

END SUBROUTINE

