!! PURPOSE:
!!
!!   This subroutine initialize cluster Lagrange basis functions
!!   in 3D system
!!
!! AUTHOR:
!!
!!   Fadjar Fathurrahman

SUBROUTINE init_LF3d_c( NN, AA, BB )

  USE m_LF3d

  IMPLICIT NONE
  
  !! Number of sampling points in x, y, and z directions
  INTEGER :: NN(3)

  !! End points of the simulation box (or unit cell in the periodic case)
  REAL(8) :: AA(3), BB(3)
  
  !! Local
  INTEGER :: Nx, Ny, Nz
  REAL(8) :: Lx, Ly, Lz
  INTEGER :: i, j, k, ip

  ! Set the type of LF3d
  LF3d_TYPE = LF3d_CLUSTER

  LF3d_NN(:) = NN(:)
  LF3d_AA(:) = AA(:)
  LF3d_BB(:) = BB(:)

  ! Check if BB(i) > AA(i)
  DO i = 1,3
    IF( BB(i) <= AA(i) ) THEN
      WRITE(*,*) 'ERROR in init_LF3d_c:'
      WRITE(*,'(1x,A,I5,2F18.10)') 'AA(i) <= BB(i): ', i, AA(i), BB(i)
      STOP
    ENDIF 
  ENDDO

  LF3d_LL(:) = BB(:) - AA(:)
  LF3d_hh(:) = LF3d_LL(:)/( NN(:) + 1 )

  LF3d_dVol = LF3d_hh(1) * LF3d_hh(2) * LF3d_hh(3)
  LF3d_Npoints = LF3d_NN(1) * LF3d_NN(2) * LF3d_NN(3)

  ! Shortcuts
  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)
  !
  Lx = LF3d_LL(1)
  Ly = LF3d_LL(2)
  Lz = LF3d_LL(3)

  ! Initialize grid points
  ALLOCATE( LF3d_grid_x( Nx ) )
  ALLOCATE( LF3d_grid_y( Ny ) )
  ALLOCATE( LF3d_grid_z( Nz ) )

  !
  CALL init_grid_1d_c( Nx, AA(1), BB(1), LF3d_grid_x )
  CALL init_grid_1d_c( Ny, AA(2), BB(2), LF3d_grid_y )
  CALL init_grid_1d_c( Nz, AA(3), BB(3), LF3d_grid_z )
  
  ! Initialize matrices D2jl
  !
  ! NOTE: We don't need D1jl matrices for the moment
  ALLOCATE( LF3d_D2jl_x( Nx, Nx ) )
  ALLOCATE( LF3d_D2jl_y( Ny, Ny ) )
  ALLOCATE( LF3d_D2jl_z( Nz, Nz ) )

  CALL init_deriv_matrix_c( Nx, Lx, LF3d_D2jl_x )
  CALL init_deriv_matrix_c( Ny, Ly, LF3d_D2jl_y )
  CALL init_deriv_matrix_c( Nz, Lz, LF3d_D2jl_z )

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

