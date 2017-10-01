SUBROUTINE read_atomic_positions(filein)
  USE m_constants, ONLY : ANG2BOHR
  USE m_input_vars
  IMPLICIT NONE 
  CHARACTER(*) :: filein
  CHARACTER(128) :: line, line_unit
  INTEGER :: ia

  OPEN(unit=IU,file=filein,status='old')

  DO WHILE(.true.)
    READ(IU,*) line, line_unit
    IF( line(1:16) == 'ATOMIC_POSITIONS' ) GOTO 2909
  ENDDO 

  WRITE(*,*)
  WRITE(*,*) 'ERROR: no ATOMIC_POSITIONS is found in the input file:'
  STOP 

  2909 CONTINUE 
!  WRITE(*,*) 'ATOMIC_POSITIONS is read!'
!  WRITE(*,*) line_unit
  
  !
  ALLOCATE( in_atmsymb(nat) )
  ALLOCATE( in_pos(3,nat) )
  !
  DO ia = 1, nat
    READ(IU,*) in_atmsymb(ia), in_pos(1,ia), in_pos(2,ia), in_pos(3,ia)
    !WRITE(*,'(1x,A5,3F18.10)') trim(in_atmsymb(ia)), in_pos(1:3,ia)
  ENDDO 

  CLOSE(IU)

  IF( trim(line_unit) == 'angstrom' ) THEN 
    in_pos(:,:) = ANG2BOHR*in_pos(:,:)
    WRITE(*,*)
    WRITE(*,*) 'Atomic positions are given in angstrom.'
    WRITE(*,*) 'The program have converted them to bohr.'
  ELSE 
    WRITE(*,*)
    WRITE(*,*) 'Atomic positions are given in bohr.'
  ENDIF 

END SUBROUTINE 

