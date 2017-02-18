! apply Laplacian to input vector v
SUBROUTINE apply_laplacian( v, Lv )

  USE m_LF3d, ONLY : &
                    Npoints => LF3d_Npoints, &
                    NN => LF3d_NN, &
                    lin2xyz => LF3d_lin2xyz, &
                    xyz2lin => LF3d_xyz2lin, &
                    D2jl_x => LF3d_D2jl_x, &
                    D2jl_y => LF3d_D2jl_y, &
                    D2jl_z => LF3d_D2jl_z
  IMPLICIT NONE
  ! arguments
  REAL(8) :: v(Npoints)
  REAL(8) :: Lv(Npoints)
  !
  INTEGER :: i, j, k, ip, ii, jj, kk
  INTEGER :: Nx, Ny, Nz
  !

  ! Shortcuts
  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  DO ip = 1, Npoints
    i = lin2xyz(1,ip)
    j = lin2xyz(2,ip)
    k = lin2xyz(3,ip)
    !
    Lv(ip) = 0.d0
    !
    DO ii = 1, Nx
      Lv(ip) = Lv(ip) + D2jl_x(ii,i) * v( xyz2lin(ii,j,k) )
    ENDDO
    !
    DO jj = 1, Ny
      Lv(ip) = Lv(ip) + D2jl_y(jj,j) * v( xyz2lin(i,jj,k) )
    ENDDO
    !
    DO kk = 1, Nz
      Lv(ip) = Lv(ip) + D2jl_z(kk,k) * v( xyz2lin(i,j,kk) )
    ENDDO
  ENDDO

END SUBROUTINE

