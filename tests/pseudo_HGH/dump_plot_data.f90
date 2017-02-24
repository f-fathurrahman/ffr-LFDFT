
SUBROUTINE dump_plot_data( ps )
  USE ps_hgh_m
  IMPLICIT NONE
  TYPE(hgh_t) :: ps
  !
  INTEGER, PARAMETER :: NPTS = 8000
  REAL(8) :: r, dr, f, r0, f1, f2, f3
  INTEGER :: ir, iu, il
  CHARACTER(256) :: filname

  dr = 0.01d0

  ! Local potential
  filname = trim(ps%atom_name)//'_Vlocal.dat'
  iu = 331
  OPEN(unit=iu,file=filname,action='write')
  r0 = 1.d-8
  DO ir = 1, NPTS
    r = r0 + dr*(ir-1)
    f = vlocalr_scalar( r, ps )
    WRITE(iu,'(2F18.10)') r, f
  ENDDO
  CLOSE(iu)
  
  ! A projector function
  !filname = trim(ps%atom_name)//'_proj.dat'
  DO il = 0, ps%l_max
    WRITE(filname,'(A,I1,A)') trim(ps%atom_name)//'_proj_',il,'.dat'
    iu = 331
    OPEN(unit=iu,file=filname,action='write')
    r0 = 1.d-8
    DO ir = 1, NPTS/20
      r = r0 + dr*(ir-1)
      f1 = projectorr_scalar( r, ps, 1, il )
      f2 = projectorr_scalar( r, ps, 2, il )
      f3 = projectorr_scalar( r, ps, 3, il )
      WRITE(iu,'(4F18.10)') r, f1, f2, f3
    ENDDO
    CLOSE(iu)
  ENDDO
END SUBROUTINE

