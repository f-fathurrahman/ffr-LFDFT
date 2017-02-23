PROGRAM test_hgh
  USE m_ps_hgh, ONLY : ps_hgh_t, &
                       hgh_init, &
                       hgh_info
  
  IMPLICIT NONE
  TYPE(ps_hgh_t) :: ps
  INTEGER :: i, j, l
  CHARACTER(64) :: filename

  CALL getarg(1,filename)

  CALL hgh_init( ps, filename )

  CALL hgh_info( ps )

  !WRITE(*,*) 'ps%atom_name = ', trim(ps%atom_name)
  !WRITE(*,'(1x,A,I5)') 'ps%z_val = ', ps%z_val
  !WRITE(*,'(1x,A,I5)') 'ps%l_max = ', ps%l_max
  !WRITE(*,'(1x,A,F18.10)') 'ps%rlocal = ', ps%rlocal
  !
  !WRITE(*,*) 'ps%rc = '
  !DO i = 0, 3
  !  WRITE(*,'(6x,I5,F18.10)') i, ps%rc(i)
  !ENDDO
  !
  !WRITE(*,*) 'ps%c = '
  !DO i = 1, 4
  !  WRITE(*,'(6x,I5,F18.10)') i, ps%c(i)
  !ENDDO
  !
  !DO l = 0, ps%l_max
  !  WRITE(*,*)
  !  WRITE(*,*) 'Matrix h for l = ', l
  !  DO i = 1, 3
  !    WRITE(*,*)
  !    DO j = 1, 3
  !      WRITE(*,'(F18.10)',advance='no') ps%h(l,i,j)
  !    ENDDO
  !  ENDDO
  !  WRITE(*,*)
  !ENDDO

  CALL dump_plot_data( ps )

END PROGRAM

SUBROUTINE dump_plot_data( ps )
  USE m_ps_hgh, ONLY : ps_hgh_t, hgh_eval_Vloc, hgh_eval_proj
  IMPLICIT NONE
  TYPE(ps_hgh_t) :: ps
  !
  INTEGER, PARAMETER :: NPTS = 8000
  REAL(8) :: r, dr, f, r0, f1, f2, f3
  INTEGER :: ir, iu
  CHARACTER(256) :: filname

  dr = 0.01d0

  ! Local potential
  filname = trim(ps%atom_name)//'_Vlocal.dat'
  iu = 331
  OPEN(unit=iu,file=filname,action='write')
  r0 = 1.d-8
  DO ir = 1, NPTS
    r = r0 + dr*(ir-1)
    f = hgh_eval_Vloc( ps, r )
    WRITE(iu,'(2F18.10)') r, f
  ENDDO
  CLOSE(iu)
  
  ! A projector function
  filname = trim(ps%atom_name)//'_proj.dat'
  iu = 331
  OPEN(unit=iu,file=filname,action='write')
  r0 = 1.d-8
  DO ir = 1, NPTS/20
    r = r0 + dr*(ir-1)
    f1 = hgh_eval_proj( ps, 0, 1, r )
    f2 = hgh_eval_proj( ps, 0, 2, r )
    f3 = hgh_eval_proj( ps, 0, 3, r )
    WRITE(iu,'(4F18.10)') r, f1, f2, f3
  ENDDO
  CLOSE(iu)
END SUBROUTINE

