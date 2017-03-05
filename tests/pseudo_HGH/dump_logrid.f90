
SUBROUTINE dump_logrid( ps )
  USE ps_hgh_m
  IMPLICIT NONE
  TYPE(hgh_t) :: ps
  !
  REAL(8) :: r, dr, f, r0
  INTEGER :: ir, iu, il
  CHARACTER(256) :: filname
  INTEGER :: NPTS

  NPTS = size(ps%g%rofi)

  ! Local potential
  filname = trim(ps%atom_name)//'_Vlocal.dat'
  iu = 331
  OPEN(unit=iu,file=filname,action='write')
  DO ir = 1, NPTS
    WRITE(iu,'(2F18.10)') ps%g%rofi(ir), ps%Vlocal(ir)
  ENDDO
  CLOSE(iu)
  
  ! KB projectors
  DO il = 0, ps%l_max
    WRITE(filname,'(A,I1,A)') trim(ps%atom_name)//'_proj_',il,'.dat'
    iu = 331
    OPEN(unit=iu,file=filname,action='write')
    DO ir = 1, NPTS
      WRITE(iu,'(4F18.10)') ps%g%rofi(ir), ps%kb(ir,il,1), ps%kb(ir,il,2), ps%kb(ir,il,3)
    ENDDO
    CLOSE(iu)
  ENDDO

END SUBROUTINE

