SUBROUTINE hgh_info( ps )
  USE ps_hgh_m
  IMPLICIT NONE
  !
  TYPE(hgh_t) :: ps
  INTEGER :: i, j, l

  WRITE(*,*) 'atom_name = ', trim(ps%atom_name)
  WRITE(*,'(1x,A,I5)') 'z_val = ', ps%z_val
  WRITE(*,'(1x,A,I5)') 'l_max = ', ps%l_max
  WRITE(*,'(1x,A,F18.10)') 'rlocal = ', ps%rlocal
  !
  WRITE(*,*)
  IF( ps%l_max > 0 ) THEN
    WRITE(*,*) 'rc (for non-local potential) = '
    DO i = 0, ps%l_max
      WRITE(*,'(1x,I1,F18.10)') i, ps%rc(i)
    ENDDO
  ENDIF
  !
  WRITE(*,*)
  IF( ps%l_max > 0 ) THEN
    WRITE(*,*) 'kbr (for non-local potential) = '
    DO l = 0, ps%l_max
      WRITE(*,'(1x,I1,F18.10)') l, ps%kbr(l)
    ENDDO
  ENDIF
  !
  WRITE(*,*)
  WRITE(*,*) 'c (for local potential ) = '
  DO i = 1, 4
    WRITE(*,'(1x,I1,F18.10)') i, ps%c(i)
  ENDDO
  !
  DO l = 0, ps%l_max
    WRITE(*,*)
    WRITE(*,*) 'Matrix h for l = ', l
    DO i = 1, 3
      WRITE(*,*)
      DO j = 1, 3
        WRITE(*,'(F18.10)',advance='no') ps%h(l,i,j)
      ENDDO
    ENDDO
    WRITE(*,*)
  ENDDO

  IF( allocated(ps%Vlocal) ) WRITE(*,*) 'shape(ps%Vlocal) = ', shape(ps%Vlocal)
  IF( allocated(ps%kb) ) WRITE(*,*) 'shape(ps%kb) = ', shape(ps%kb)

  WRITE(*,*) 'shape(ps%g%rofi) = ', shape(ps%g%rofi)

END SUBROUTINE

