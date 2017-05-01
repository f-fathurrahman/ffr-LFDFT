SUBROUTINE add_V_ps_loc_G( isp )
  
  USE m_PsPot, ONLY : Ps => Ps_HGH_Params

  USE m_Ps_HGH, ONLY : hgh_eval_Vloc_G

  USE m_hamiltonian, ONLY : V_ps_loc

  USE m_atoms, ONLY : Nspecies, strf => StructureFactor

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     G2 => LF3d_G2, &
                     NN => LF3d_NN, &
                     LL => LF3d_LL, &
                     xyz2lin => LF3d_xyz2lin, &
                     grid_x => LF3d_grid_x, &
                     grid_y => LF3d_grid_y, &
                     grid_z => LF3d_grid_z

  USE m_constants, ONLY : PI
  USE bspline_module

  IMPLICIT NONE 
  INTEGER :: ip, isp, Nx, Ny, Nz
  REAL(8) :: Gm, Omega
  COMPLEX(8), ALLOCATABLE :: ctmp(:)
  !
  REAL(8), ALLOCATABLE :: interp_V(:,:,:)
  INTEGER :: i,j,k, ii,jj,kk, iflag
  REAL(8), ALLOCATABLE :: tx(:), ty(:), tz(:)
  INTEGER :: kx, ky, kz
  REAL(8), ALLOCATABLE :: x(:), y(:), z(:)
  REAL(8) :: shiftx, shifty, shiftz
  INTEGER :: iknot
  INTEGER :: idx, idy, idz, iloy, iloz, inbvx, inbvy, inbvz
  REAL(8) :: val

  WRITE(*,*)
  WRITE(*,*) 'Adding V_ps_loc via G-space'

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
  DO i = 1,Nx
    x(i) = grid_x(i) - shiftx
  ENDDO 
  x(Nx+1) = x(Nx) + 2*shiftx

  shifty = 0.5d0*( grid_y(2) - grid_y(1) )
  DO j = 1,Ny
    y(j) = grid_y(j) - shifty
  ENDDO 
  y(Ny+1) = y(Ny) + 2*shifty

  shiftz = 0.5d0*( grid_z(2) - grid_z(1) )
  DO k = 1,Nz
    z(k) = grid_z(k) - shiftz
  ENDDO 
  z(Nz+1) = z(Nz) + 2*shiftz

  WRITE(*,'(1x,A,3F18.10)') 'shift:', shiftx, shifty, shiftz

  ! Cell volume
  Omega = LL(1) * LL(2) * LL(3)

  ALLOCATE( ctmp(Npoints) )
  ctmp(:) = cmplx(0.d0,0.d0,kind=8)

  WRITE(*,*) 'isp = ', isp
  WRITE(*,*) 'zval = ', Ps(isp)%zval
  WRITE(*,*) 'sum(strf) = ', sum(strf)
  DO ip = 1,Npoints
    Gm = sqrt(G2(ip))
    ctmp(ip) = hgh_eval_Vloc_G( Ps(isp), Gm ) / Omega
    !ctmp(ip) = -4.d0*PI*Ps(isp)%zval/G2(ip)/Omega
  ENDDO

  ! inverse FFT: G -> R
  CALL fft_fftw3( ctmp, Nx, Ny, Nz, .true. )

  ALLOCATE( interp_V(Nx+1,Ny+1,Nz+1) )
  ! copy to interp_V
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
    interp_V(ii,jj,kk) = real( ctmp(ip), kind=8 )
  ENDDO 
  ENDDO 
  ENDDO 

  WRITE(*,*) 'sum(ctmp) = ', sum(real(ctmp))

  j = 1
  k = 1
  DO i = 1,Nx+1
    WRITE(101,*) x(i), interp_V(i,j,k)
  ENDDO 

  i = 1
  k = 1
  DO j = 1,Ny+1
    WRITE(102,*) y(j), interp_V(i,j,k)
  ENDDO 

  i = 1
  j = 1
  DO k = 1,Nz+1
    WRITE(103,*) z(k), interp_V(i,j,k)
  ENDDO 

  CALL db3ink( x, Nx+1, y, Ny+1, z, Nz+1, interp_V, kx,ky,kz, iknot, tx,ty,tz, interp_V, iflag )
  WRITE(*,*) 'iflag = ', iflag

  iloy = 1
  iloz = 1
  inbvx = 1
  inbvy = 1
  inbvz = 1
  idx = 0
  idy = 0
  idz = 0

  DO i = 1,Nx
    CALL db3val( abs(grid_x(i)-8.d0), grid_y(1), grid_z(1), &
                 idx,idy,idz, tx,ty,tz, Nx+1,Ny+1,Nz+1,kx,ky,kz, interp_V,&
                 val, iflag, inbvx, inbvy, inbvz, iloy, iloz )
    !WRITE(*,*) 'iflag = ', iflag
    WRITE(104,*) grid_x(i), val
  ENDDO 


  DEALLOCATE( interp_V, tx, ty, tz )
  DEALLOCATE( ctmp )

END SUBROUTINE 

