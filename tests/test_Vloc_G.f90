PROGRAM test_Vloc

  USE m_constants, ONLY : PI
  USE m_Ps_HGH

  IMPLICIT NONE
  TYPE(Ps_HGH_Params_T) :: ps
  CHARACTER(64) :: filename
  INTEGER, PARAMETER :: Nradial = 1000
  REAL(8), ALLOCATABLE :: ff_Vloc(:)
  REAL(8) :: r, dr, rstart
  INTEGER :: i

  IF( iargc() /= 1 ) THEN 
    WRITE(*,*) 'Exactly one argument is needed'
    STOP 
  ENDIF 

  CALL getarg(1,filename)

  CALL init_Ps_HGH_Params( ps, filename )

  CALL info_Ps_HGH_Params( ps )

  ALLOCATE( ff_Vloc(Nradial) )

  dr = 0.01
  rstart = 0.001
  DO i = 1, Nradial
    r = (i-1)*dr + rstart
    ff_Vloc(i) = hgh_eval_Vloc_G( ps, r )
    WRITE(222,'(1x,F18.10,ES18.10)') r, ff_Vloc(i) 
  ENDDO 

  DEALLOCATE( ff_Vloc )

END PROGRAM


