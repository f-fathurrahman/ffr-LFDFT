PROGRAM test_hgh
  USE ps_hgh_m
  IMPLICIT NONE
  TYPE(hgh_t) :: ps
  CHARACTER(64) :: filename

  CALL getarg(1,filename)

  CALL hgh_init( ps, filename )

  CALL hgh_process(ps)

  CALL hgh_info( ps )

  CALL dump_plot_data( ps )

  CALL hgh_end( ps )
END PROGRAM


