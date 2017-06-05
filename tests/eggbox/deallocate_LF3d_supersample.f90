SUBROUTINE dealloc_LF3d_supersample()

  USE m_LF3d_supersample
  IMPLICIT NONE 

  IF( allocated(LF3d_grid_x_ss) ) DEALLOCATE(LF3d_grid_x_ss)
  IF( allocated(LF3d_grid_y_ss) ) DEALLOCATE(LF3d_grid_y_ss)
  IF( allocated(LF3d_grid_z_ss) ) DEALLOCATE(LF3d_grid_z_ss)

  IF( allocated(LF3d_lingrid_ss) ) DEALLOCATE(LF3d_lingrid_ss)
  IF( allocated(LF3d_xyz2lin_ss) ) DEALLOCATE(LF3d_xyz2lin_ss)
  IF( allocated(LF3d_lin2xyz_ss) ) DEALLOCATE(LF3d_lin2xyz_ss)

END SUBROUTINE 

