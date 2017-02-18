! efefer 28 December 2015


! Real part of a planewave exp(G*x)
FUNCTION pw1d( G, x )
  IMPLICIT NONE
  !
  REAL(8) :: pw1d, G, x
  !
  pw1d = cos(G*x)
END FUNCTION


! only implements the cos functions, no cross terms
FUNCTION pw3d( Gx, Gy, Gz, x, y, z )
  IMPLICIT NONE
  !
  REAL(8) :: pw3d, Gx,Gy,Gz, x,y,z
  !
  pw3d = cos(Gx*x) * cos(Gy*y) * cos(Gz*z)
END FUNCTION


FUNCTION dx_pw3d( Gx, Gy, Gz, x, y, z )
  IMPLICIT NONE
  !
  REAL(8) :: dx_pw3d, Gx,Gy,Gz, x,y,z
  !
  dx_pw3d = -Gx*sin(Gx*x) * cos(Gy*y) * cos(Gz*z)
END FUNCTION


FUNCTION d2x_pw3d( Gx, Gy, Gz, x, y, z )
  IMPLICIT NONE
  !
  REAL(8) :: d2x_pw3d, Gx,Gy,Gz, x,y,z
  !
  d2x_pw3d = -Gx**2 * cos(Gx*x) * cos(Gy*y) * cos(Gz*z)
END FUNCTION



PROGRAM t_LF3d
  USE m_LF3d
  USE m_LF1d
  IMPLICIT NONE
  !
  INTEGER :: Nx, Ny, Nz
  REAL(8) :: Lx, Ly, Lz
  TYPE(LF3d_t) :: LF3
  !
  REAL(8), ALLOCATABLE :: coef(:), cx(:), cy(:), cz(:)
  REAL(8), ALLOCATABLE :: dcoef(:), dcx(:)
  !
  REAL(8) :: Gx, Gy, Gz, x, y, z
  INTEGER :: ip, i, j, k
  !
  REAL(8) :: pw1d, dx_pw3d

  Nx = 5
  Ny = 7
  Nz = 3

  Lx = 5.d0
  Ly = 5.d0
  Lz = 3.d0

  CALL init_LF3d_p( LF3, (/Nx, Ny, Nz/), (/0.d0,0.d0,0.d0/), (/Lx, Ly, Lz/) )

  CALL info_LF1d( LF3%LFx, .TRUE. )
  CALL info_LF1d( LF3%LFy, .TRUE. )
  CALL info_LF1d( LF3%LFz, .TRUE. )

  ALLOCATE( coef(Nx*Ny*Nz), dcoef(Nx*Ny*Nz) )
  ALLOCATE( cx(Nx), cy(Ny), cz(Nz) )
  ALLOCATE( dcx(Nx) )

  Gx = 1/Lx*2.d0*PI
  Gy = 2/Ly*2.d0*PI
  Gz = 3/Lz*2.d0*PI

  DO i = 1,Nx
    x = LF3%LFx%grid(i)
    cx(i) = pw1d( Gx, x )
    WRITE(*,*) i, x, cx(i)
  ENDDO

  DO j = 1,Ny
    y = LF3%LFy%grid(j)
    cy(j) = pw1d( Gy, LF3 % LFy % grid(j) )
    WRITE(*,*) j, y, cy(j)
  ENDDO

  DO k = 1,Nz
    z = LF3%LFy%grid(k)
    cz(k) = pw1d( Gz, LF3%LFz%grid(k) )
    WRITE(*,*) k, z, cz(k)
  ENDDO

  WRITE(*,*) 'Linear grid'
  DO ip = 1,Nx*Ny*Nz
    i = LF3%lin2xyz( 1, ip )
    j = LF3%lin2xyz( 2, ip )
    k = LF3%lin2xyz( 3, ip )
    coef(ip) = cx(i) * cy(j) * cz(k)
    WRITE(*,'(1x,I8,F18.10)') ip, coef(ip)
  ENDDO

  ! operate d/dx to coef
  dcx = matmul( LF3%LFx%D1jl, cx )
  WRITE(*,*) dcx

  ! FIXME access x, y, z through lingrid
  DO ip=1,Nx*Ny*Nz
    i = LF3%lin2xyz( 1, ip )
    j = LF3%lin2xyz( 2, ip )
    k = LF3%lin2xyz( 3, ip )
    dcoef(ip) = dcx(i) * cy(j) * cz(k)
    !
    x = LF3%LFx%grid(i)
    y = LF3%LFy%grid(j)
    z = LF3%LFz%grid(k)
    WRITE(*,'(1x,2F18.10)') dcoef(ip), dx_pw3d(Gx,Gy,Gz,x,y,z)
  ENDDO

  DEALLOCATE( coef, dcoef )
  DEALLOCATE( cx, cy, cz )
  DEALLOCATE( dcx )

END PROGRAM

