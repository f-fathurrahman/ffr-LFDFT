!--------------------------------
SUBROUTINE init_LF3d_p( NN, AA, BB )
!--------------------------------
  USE m_LF3d
  IMPLICIT NONE
  !
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  !
  INTEGER :: Nx, Ny, Nz
  REAL(8) :: Lx, Ly, Lz
  !
  INTEGER :: i, j, k, ip


  LF3d_NN(:) = NN(:)
  LF3d_AA(:) = AA(:)
  LF3d_BB(:) = BB(:)
  LF3d_LL(:) = BB(:) - AA(:)  ! TODO: Check if BB > AA

  LF3d_hh(:) = LF3d_LL(:)/N

  LF3d_dVol = LF3d_hh(1) * LF3d_hh(2) * LF3d_hh(3)
  LF3d_Npoints = LF3d_NN(1) * LF3d_NN(2) * LF3d_NN(3)

  ! Shortcuts
  Nx = NN(1); Ny = NN(2); Nz = NN(3)
  Lx = BB(1) - AA(1);  Ly = BB(2) - AA(2);  Lz = BB(3) - AA(3)
  
  ALLOCATE( LF3d_grid_x( Nx ) )
  ALLOCATE( LF3d_grid_y( Ny ) )
  ALLOCATE( LF3d_grid_z( Nz ) )
  !
  CALL init_grid_1d_p( Nx, AA(1), BB(1), LF3d_grid_x )
  CALL init_grid_1d_p( Ny, AA(2), BB(2), LF3d_grid_y )
  CALL init_grid_1d_p( Nz, AA(3), BB(3), LF3d_grid_z )
  
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

END SUBROUTINE

