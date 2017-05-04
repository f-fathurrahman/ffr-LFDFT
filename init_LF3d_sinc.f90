!! PURPOSE:
!!
!!   This subroutine initialize sinc Lagrange basis functions
!!   in 3D system
!!
!! AUTHOR:
!!
!!   Fadjar Fathurrahman

SUBROUTINE init_LF3d_sinc( NN, hh )

  USE m_LF3d

  IMPLICIT NONE
  
  !! Number of sampling points in x, y, and z directions
  INTEGER :: NN(3)

  REAL(8) :: hh(3)

  !! Local
  INTEGER :: Nx, Ny, Nz
  INTEGER :: i, j, k, ip

  ! Set the type of LF3d
  LF3d_TYPE = LF3d_SINC

  LF3d_NN(:) = NN(:)

  ! Is there any point of having different NN for x, y, and z ?
  LF3d_AA(:) = -( NN(:) - 1 )/2.d0 * hh(:)
  LF3d_BB(:) =  ( NN(:) - 1 )/2.d0 * hh(:)
  LF3d_hh(:) = hh(:)

  LF3d_LL(:) = LF3d_BB(:) - LF3d_AA(:)

  LF3d_dVol = LF3d_hh(1) * LF3d_hh(2) * LF3d_hh(3)
  LF3d_Npoints = LF3d_NN(1) * LF3d_NN(2) * LF3d_NN(3)

  ! Shortcuts
  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  ! Initialize grid points
  ALLOCATE( LF3d_grid_x( Nx ) )
  ALLOCATE( LF3d_grid_y( Ny ) )
  ALLOCATE( LF3d_grid_z( Nz ) )

  !
  CALL init_grid_1d_sinc( Nx, hh(1), LF3d_grid_x )
  CALL init_grid_1d_sinc( Ny, hh(2), LF3d_grid_y )
  CALL init_grid_1d_sinc( Nz, hh(3), LF3d_grid_z )
  
  ! Initialize matrices D2jl
  !
  ! NOTE: We don't need D1jl matrices for the moment
  ALLOCATE( LF3d_D2jl_x( Nx, Nx ) )
  ALLOCATE( LF3d_D2jl_y( Ny, Ny ) )
  ALLOCATE( LF3d_D2jl_z( Nz, Nz ) )

  CALL init_deriv_matrix_sinc( Nx, hh(1), LF3d_grid_x, LF3d_D2jl_x )
  CALL init_deriv_matrix_sinc( Ny, hh(2), LF3d_grid_y, LF3d_D2jl_y )
  CALL init_deriv_matrix_sinc( Nz, hh(3), LF3d_grid_z, LF3d_D2jl_z )

  !
  ! 3D mapping stuffs
  !
  ALLOCATE( LF3d_lingrid( 3, Nx*Ny*Nz ) )
  ALLOCATE( LF3d_xyz2lin( Nx, Ny, Nz ) )
  ALLOCATE( LF3d_lin2xyz( 3, Nx*Ny*Nz ) )
  ip = 0
  DO k = 1, Nz
    DO j = 1, Ny
      DO i = 1, Nx
        ip = ip + 1
        LF3d_lingrid( 1, ip ) = LF3d_grid_x(i)
        LF3d_lingrid( 2, ip ) = LF3d_grid_y(j)
        LF3d_lingrid( 3, ip ) = LF3d_grid_z(k)
        !
        LF3d_xyz2lin( i, j, k ) = ip
        LF3d_lin2xyz( 1:3, ip ) = (/ i, j, k /)
      ENDDO
    ENDDO
  ENDDO

  ! G-vectors
  CALL init_gvec()

END SUBROUTINE

