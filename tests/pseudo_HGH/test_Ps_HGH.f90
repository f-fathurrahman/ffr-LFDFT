PROGRAM test_hgh

  USE m_Ps_HGH

  IMPLICIT NONE
  TYPE(Ps_HGH_Params) :: ps
  CHARACTER(64) :: filename

  IF( iargc() /= 1 ) THEN 
    WRITE(*,*) 'Exactly one argument is needed'
    STOP 
  ENDIF 

  CALL getarg(1,filename)

  CALL hgh_init( ps, filename )

  CALL hgh_info( ps )
  !CALL dump_plot_data( ps )
  !CALL dump_logrid( ps )

END PROGRAM


