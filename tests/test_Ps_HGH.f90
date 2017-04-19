PROGRAM test_hgh

  USE m_Ps_HGH

  IMPLICIT NONE
  TYPE(Ps_HGH_Params_T) :: ps
  CHARACTER(64) :: filename

  IF( iargc() /= 1 ) THEN 
    WRITE(*,*) 'Exactly one argument is needed'
    STOP 
  ENDIF 

  CALL getarg(1,filename)

  CALL init_Ps_HGH_Params( ps, filename )

  CALL info_Ps_HGH_Params( ps )

  !CALL dump_plot_data( ps )
  !CALL dump_logrid( ps )

END PROGRAM


