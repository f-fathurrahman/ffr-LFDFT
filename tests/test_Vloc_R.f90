PROGRAM test_hgh

  USE m_Ps_HGH

  IMPLICIT NONE
  TYPE(Ps_HGH_Params_T) :: ps
  CHARACTER(64) :: filename
  REAL(8) :: dr, VlocR, r, r0
  INTEGER :: i
  INTEGER :: iargc

  IF( iargc() /= 1 ) THEN 
    WRITE(*,*) 'Exactly one argument is needed'
    STOP 
  ENDIF 

  CALL getarg(1,filename)

  CALL init_Ps_HGH_Params( ps, filename )

  CALL info_Ps_HGH_Params( ps )

  r0 = epsilon(1.d0)
  dr = 1d-8
  DO i = 1,10000
    r = r0 + (i-1)*dr
    VlocR = hgh_eval_Vloc_R( ps, r )
    WRITE(222,'(1x,2F18.10)') r, VlocR
  ENDDO 

END PROGRAM


