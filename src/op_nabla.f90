! apply gradient operator to input vector v
SUBROUTINE op_nabla( v, Lv )

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     NN => LF3d_NN, &
                     lin2xyz => LF3d_lin2xyz, &
                     xyz2lin => LF3d_xyz2lin, &
                     D1jl_x => LF3d_D1jl_x, &
                     D1jl_y => LF3d_D1jl_y, &
                     D1jl_z => LF3d_D1jl_z
  IMPLICIT NONE
  ! arguments
  REAL(8) :: v(Npoints)
  REAL(8) :: Lv(3,Npoints)
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
    Lv(1,ip) = 0.d0
    DO ii = 1, Nx
      Lv(1,ip) = Lv(1,ip) + D1jl_x(ii,i) * v( xyz2lin(ii,j,k) )
    ENDDO
    !
    Lv(2,ip) = 0.d0
    DO jj = 1, Ny
      Lv(2,ip) = Lv(2,ip) + D1jl_y(jj,j) * v( xyz2lin(i,jj,k) )
    ENDDO
    !
    Lv(3,ip) = 0.d0
    DO kk = 1, Nz
      Lv(3,ip) = Lv(3,ip) + D1jl_z(kk,k) * v( xyz2lin(i,j,kk) )
    ENDDO

  ENDDO

END SUBROUTINE



! apply divergence (nabla dot) operator to input vector v
SUBROUTINE op_nabla_dot( v, Lv )

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     NN => LF3d_NN, &
                     lin2xyz => LF3d_lin2xyz, &
                     xyz2lin => LF3d_xyz2lin, &
                     D1jl_x => LF3d_D1jl_x, &
                     D1jl_y => LF3d_D1jl_y, &
                     D1jl_z => LF3d_D1jl_z
  IMPLICIT NONE
  ! arguments
  REAL(8) :: v(3,Npoints)
  REAL(8) :: Lv(Npoints)
  !
  INTEGER :: i, j, k, ip, ii, jj, kk
  INTEGER :: Nx, Ny, Nz
  REAL(8) :: Lx, Ly, Lz
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
    Lx = 0.d0
    DO ii = 1, Nx
      Lx = Lx + D1jl_x(ii,i) * v( 1, xyz2lin(ii,j,k) )
    ENDDO
    !
    Ly = 0.d0
    DO jj = 1, Ny
      Ly = Ly + D1jl_y(jj,j) * v( 2, xyz2lin(i,jj,k) )
    ENDDO
    !
    Lz = 0.d0
    DO kk = 1, Nz
      Lz = Lz + D1jl_z(kk,k) * v( 3, xyz2lin(i,j,kk) )
    ENDDO

    Lv(ip) = Lx + Ly + Lz

  ENDDO

END SUBROUTINE

