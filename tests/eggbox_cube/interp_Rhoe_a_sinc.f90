SUBROUTINE interp_Rhoe_a_sinc( Rhoe, Rhoe_a)
  USE m_LF3d, ONLY : LF3d_TYPE, LF3d_SINC
  USE m_grid_atom_cube
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     NN => LF3d_NN, &
                     xyz2lin => LF3d_xyz2lin, &
                     x => LF3d_grid_x, &
                     y => LF3d_grid_y, &
                     z => LF3d_grid_z

  USE bspline
  IMPLICIT NONE 
  REAL(8) :: Rhoe(Npoints)
  REAL(8) :: Rhoe_a(Npoints_a)
  !
  INTEGER :: Nx, Ny, Nz, kx, ky, kz
  INTEGER :: iknot
  REAL(8), ALLOCATABLE :: tx(:), ty(:), tz(:)
  REAL(8), ALLOCATABLE :: interp_Rhoe(:,:,:)
  REAL(8), ALLOCATABLE :: bcoef(:,:,:)
  INTEGER :: iflag
  INTEGER :: ip
  !
  INTEGER :: i,j,k, ii,jj,kk, ip_a
  REAL(8) :: shiftx, shifty, shiftz
  INTEGER :: idx, idy, idz, iloy, iloz, inbvx, inbvy, inbvz
  REAL(8) :: val, dx, dy, dz

  IF( LF3d_TYPE /= LF3d_SINC ) THEN 
    WRITE(*,*) 'ERROR in interp_Rhoe_a_sinc, need sinc LF'
    STOP 
  ENDIF 

  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  iknot = 0

  kx = 4
  ky = 4
  kz = 4
  ALLOCATE( tx(Nx+kz), ty(Ny+ky), tz(Nz+kz) )

  Rhoe_a(:) = 0.d0

  ALLOCATE( interp_Rhoe(Nx,Ny,Nz) )
  ALLOCATE( bcoef(Nx,Ny,Nz) )

  iloy = 1
  iloz = 1
  inbvx = 1
  inbvy = 1
  inbvz = 1
  idx = 0
  idy = 0
  idz = 0

  ! copy to interp_Rhoe
  DO kk = 1,Nz
  DO jj = 1,Ny
  DO ii = 1,Nx
    i = ii
    j = jj
    k = kk
    ip = xyz2lin(i,j,k)
    interp_Rhoe(ii,jj,kk) = Rhoe(ip)
  ENDDO 
  ENDDO 
  ENDDO 

  CALL db3ink( x, Nx, y, Ny, z, Nz, interp_Rhoe, kx,ky,kz, iknot, tx,ty,tz, bcoef, iflag )
  
  DO ip_a = 1, Npoints_a
    dx = grid_a(1,ip_a)
    dy = grid_a(2,ip_a)
    dz = grid_a(3,ip_a)
    CALL db3val( dx, dy, dz, idx,idy,idz, tx,ty,tz, Nx,Ny,Nz,kx,ky,kz, bcoef,&
                 val, iflag, inbvx, inbvy, inbvz, iloy, iloz )
    Rhoe_a(ip_a) = val
  ENDDO 

  WRITE(*,*) 'integ(Rhoe_a) = ', sum(Rhoe_a)*dVol_a

  DEALLOCATE( tx, ty, tz )
  DEALLOCATE( interp_Rhoe )
  DEALLOCATE( bcoef )

END SUBROUTINE 
