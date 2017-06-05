SUBROUTINE init_LF3d_supersample( Nsupersample_ )

  USE m_LF3d, ONLY : AA => LF3d_AA, BB => LF3d_BB, LF3d_NN, LF3d_Npoints, LF3d_hh
  USE m_LF3d_supersample
  !
  IMPLICIT NONE
  !
  INTEGER :: Nsupersample_
  !! Local
  INTEGER :: Nx, Ny, Nz
  INTEGER :: i, j, k, ip

  ! Shortcuts
  Nx = LF3d_NN(1)
  Ny = LF3d_NN(2)
  Nz = LF3d_NN(3)

  Nsupersample = Nsupersample_
  LF3d_Npoints_ss = LF3d_Npoints*Nsupersample**3

  LF3d_hh_ss(1) = LF3d_hh(1)/Nsupersample
  LF3d_hh_ss(2) = LF3d_hh(2)/Nsupersample
  LF3d_hh_ss(3) = LF3d_hh(3)/Nsupersample

  ! Initialize grid points
  ALLOCATE( LF3d_grid_x_ss( Nsupersample*Nx ) )
  ALLOCATE( LF3d_grid_y_ss( Nsupersample*Ny ) )
  ALLOCATE( LF3d_grid_z_ss( Nsupersample*Nz ) )
  !
  CALL init_grid_1d_p( Nsupersample*Nx, AA(1), BB(1), LF3d_grid_x_ss )
  CALL init_grid_1d_p( Nsupersample*Ny, AA(2), BB(2), LF3d_grid_y_ss )
  CALL init_grid_1d_p( Nsupersample*Nz, AA(3), BB(3), LF3d_grid_z_ss )

  !
  ! 3D mapping stuffs
  !
  ALLOCATE( LF3d_lingrid_ss( 3, LF3d_Npoints_ss ) )
  ALLOCATE( LF3d_xyz2lin_ss( Nx*Nsupersample, Ny*Nsupersample, Nz*Nsupersample ) )
  ALLOCATE( LF3d_lin2xyz_ss( 3, LF3d_Npoints_ss ) )
  ip = 0
  DO k = 1, Nz*Nsupersample
    DO j = 1, Ny*Nsupersample
      DO i = 1, Nx*Nsupersample
        ip = ip + 1
        LF3d_lingrid_ss( 1, ip ) = LF3d_grid_x_ss(i)
        LF3d_lingrid_ss( 2, ip ) = LF3d_grid_y_ss(j)
        LF3d_lingrid_ss( 3, ip ) = LF3d_grid_z_ss(k)
        !
        LF3d_xyz2lin_ss( i, j, k ) = ip
        LF3d_lin2xyz_ss( 1:3, ip ) = (/ i, j, k /)
      ENDDO
    ENDDO
  ENDDO

  WRITE(*,*) 'Supersampling info:'
  WRITE(*,*) 'Nsupersample = ', Nsupersample
  WRITE(*,*) 'Npoints_ss = ', LF3d_Npoints_ss

END SUBROUTINE
