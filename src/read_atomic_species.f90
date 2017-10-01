SUBROUTINE read_atomic_species(filein)
  USE m_input_vars
  IMPLICIT NONE 
  !
  CHARACTER(*) :: filein
  CHARACTER(128) :: line
  INTEGER :: isp

  OPEN(unit=IU,file=filein,status='old')

  DO WHILE(.true.)
    READ(IU,*) line
    IF( trim(line) == 'ATOMIC_SPECIES') GOTO 2909
  ENDDO 

  WRITE(*,*)
  WRITE(*,*) 'ERROR: no ATOMIC_SPECIES is found in the input file:'
  STOP 

  2909 CONTINUE 
  !WRITE(*,*) 'ATOMIC_SPECIES is read!'
  !
  ALLOCATE( species(ntyp) )
  ALLOCATE( masses(ntyp) )
  ALLOCATE( pp_name(ntyp) )
  !
  DO isp = 1, ntyp
    READ(IU,*) species(isp), masses(isp), pp_name(isp)
    !WRITE(*,'(1x,A,F10.3,3x,A)') trim(species(isp)), masses(isp), trim(pp_name(isp))
  ENDDO 

  CLOSE(IU)

END SUBROUTINE 

