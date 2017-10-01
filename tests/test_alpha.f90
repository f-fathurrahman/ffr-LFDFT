PROGRAM test_alpha

  USE m_constants, ONLY : PI
  USE m_Ps_HGH

  IMPLICIT NONE
  TYPE(Ps_HGH_Params_T) :: ps
  CHARACTER(64) :: filename
  INTEGER, PARAMETER :: Nradial = 1000
  REAL(8), ALLOCATABLE :: ff_Vloc(:)
  REAL(8) :: r, dr
  INTEGER :: i
  INTEGER :: iargc

  IF( iargc() /= 1 ) THEN 
    WRITE(*,*) 'Exactly one argument is needed: path to pseudopotential file'
    STOP 
  ENDIF 

  CALL getarg(1,filename)

  CALL init_Ps_HGH_Params( ps, filename )

  CALL info_Ps_HGH_Params( ps )

  ALLOCATE( ff_Vloc(Nradial) )

  dr = 0.01
  DO i = 1, Nradial
    r = (i-1)*dr
    ff_Vloc(i) = hgh_eval_Vloc_4pi_r2( ps, r ) + 4.d0*PI*r*ps%zval
    WRITE(101,'(1x,F18.10,ES18.10)') r, ff_Vloc(i) 
  ENDDO 
  WRITE(*,*) 'integ = ', sum(ff_Vloc)*dr
  WRITE(*,*) 'alpha = ', ps%zval/16.d0**3 * sum(ff_Vloc)*dr

  DEALLOCATE( ff_Vloc )

END PROGRAM


