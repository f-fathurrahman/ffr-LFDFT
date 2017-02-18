MODULE m_LF3d

  INTEGER, DIMENSION(3) :: LF3d_NN
  REAL(8), DIMENSION(3) :: LF3d_LL
  REAL(8), DIMENSION(3) :: LF3d_AA, LF3d_BB
  REAL(8), DIMENSION(3) :: LF3d_hh

  INTEGER :: LF3d_Npoints
  REAL(8) :: LF3d_dVol
 
  REAL(8), ALLOCATABLE :: LF3d_grid_x(:)
  REAL(8), ALLOCATABLE :: LF3d_grid_y(:)
  REAL(8), ALLOCATABLE :: LF3d_grid_z(:)

  REAL(8), ALLOCATABLE :: LF3d_D1jl_x(:,:)
  REAL(8), ALLOCATABLE :: LF3d_D1jl_y(:,:)
  REAL(8), ALLOCATABLE :: LF3d_D1jl_z(:,:)

  REAL(8), ALLOCATABLE :: LF3d_D2jl_x(:,:)
  REAL(8), ALLOCATABLE :: LF3d_D2jl_y(:,:)
  REAL(8), ALLOCATABLE :: LF3d_D2jl_z(:,:)

  REAL(8), ALLOCATABLE :: LF3d_lingrid(:,:)
  INTEGER, ALLOCATABLE :: LF3d_xyz2lin(:,:,:)
  INTEGER, ALLOCATABLE :: LF3d_lin2xyz(:,:)

  REAL(8), ALLOCATABLE :: LF3d_G2(:)

END MODULE

