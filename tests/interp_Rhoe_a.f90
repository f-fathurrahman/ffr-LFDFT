SUBROUTINE interp_Rhoe_a( Rhoe, Rhoe_a)
  USE m_grid_atom_cube
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     NN => LF3d_NN, &
                     xyz2lin => LF3d_xyz2lin, &
                     grid_x => LF3d_grid_x, &
                     grid_y => LF3d_grid_y, &
                     grid_z => LF3d_grid_z

  USE bspline
  IMPLICIT NONE 
  REAL(8) :: Rhoe(Npoints)
  REAL(8) :: Rhoe_a(Npoints_a)
  !
  INTEGER :: Nx, Ny, Nz, kx, ky, kz
  INTEGER :: iknot
  REAL(8), ALLOCATABLE :: x(:), y(:), z(:)
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

  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  iknot = 0

  kx = 4
  ky = 4
  kz = 4
  ALLOCATE( tx(Nx+1+kz), ty(Ny+1+ky), tz(Nz+1+kz) )
  ALLOCATE( x(Nx+1), y(Ny+1), z(Nz+1) )

  shiftx = 0.5d0*( grid_x(2) - grid_x(1) )
  !shiftx = 0.d0
  DO i = 1,Nx
    x(i) = grid_x(i) - shiftx
  ENDDO 
  x(Nx+1) = x(Nx) + 2*shiftx

  shifty = 0.5d0*( grid_y(2) - grid_y(1) )
  !shifty = 0.d0
  DO j = 1,Ny
    y(j) = grid_y(j) - shifty
  ENDDO 
  y(Ny+1) = y(Ny) + 2*shifty

  shiftz = 0.5d0*( grid_z(2) - grid_z(1) )
  !shiftz = 0.d0
  DO k = 1,Nz
    z(k) = grid_z(k) - shiftz
  ENDDO 
  z(Nz+1) = z(Nz) + 2*shiftz

  Rhoe_a(:) = 0.d0

  ALLOCATE( interp_Rhoe(Nx+1,Ny+1,Nz+1) )
  ALLOCATE( bcoef(Nx+1,Ny+1,Nz+1) )

  iloy = 1
  iloz = 1
  inbvx = 1
  inbvy = 1
  inbvz = 1
  idx = 0
  idy = 0
  idz = 0

  ! copy to interp_Rhoe
  DO kk = 1,Nz+1
  DO jj = 1,Ny+1
  DO ii = 1,Nx+1
    i = ii
    j = jj
    k = kk
    IF( kk == Nz+1 ) k = 1
    IF( jj == Ny+1 ) j = 1
    IF( ii == Nx+1 ) i = 1
    ip = xyz2lin(i,j,k)
    !
    interp_Rhoe(ii,jj,kk) = Rhoe(ip)
  ENDDO 
  ENDDO 
  ENDDO 

  CALL db3ink( x, Nx+1, y, Ny+1, z, Nz+1, interp_Rhoe, kx,ky,kz, iknot, tx,ty,tz, bcoef, iflag )
!  WRITE(*,*) 'iflag = ', iflag
  
  DO ip_a = 1, Npoints_a
    dx = grid_a(1,ip_a) - shiftx
    dy = grid_a(2,ip_a) - shifty
    dz = grid_a(3,ip_a) - shiftz
    CALL db3val( dx, dy, dz, idx,idy,idz, tx,ty,tz, Nx+1,Ny+1,Nz+1,kx,ky,kz, bcoef,&
                 val, iflag, inbvx, inbvy, inbvz, iloy, iloz )
    IF( iflag /= 0 ) THEN 
      WRITE(*,*) 'ERROR in calling db3val: iflag = ', iflag
      STOP 
    ENDIF 
    Rhoe_a(ip_a) = val
  ENDDO 

  WRITE(*,*) 'integ(Rhoe_a) = ', sum(Rhoe_a)*dVol_a

  DEALLOCATE( x, y, z )
  DEALLOCATE( tx, ty, tz )
  DEALLOCATE( interp_Rhoe )
  DEALLOCATE( bcoef )

END SUBROUTINE 
