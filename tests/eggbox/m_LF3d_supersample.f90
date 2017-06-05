MODULE m_LF3d_supersample

  IMPLICIT NONE 

  INTEGER :: Nsupersample
  INTEGER :: LF3d_Npoints_ss
 
  REAL(8) :: LF3d_hh_ss(3)
  REAL(8), ALLOCATABLE :: LF3d_grid_x_ss(:)
  REAL(8), ALLOCATABLE :: LF3d_grid_y_ss(:)
  REAL(8), ALLOCATABLE :: LF3d_grid_z_ss(:)

  REAL(8), ALLOCATABLE :: LF3d_lingrid_ss(:,:)
  INTEGER, ALLOCATABLE :: LF3d_xyz2lin_ss(:,:,:)
  INTEGER, ALLOCATABLE :: LF3d_lin2xyz_ss(:,:)

END MODULE

