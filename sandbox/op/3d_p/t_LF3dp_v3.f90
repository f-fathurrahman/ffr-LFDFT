! efefer 30 December 2015


! 1D planewave exp(i*Gx*x)
FUNCTION pw1d( Gx, x )
  IMPLICIT NONE
  !
  REAL(8) :: Gx, x
  COMPLEX(8) :: pw1d
  !
  pw1d =exp( (0.d0,1.d0) * Gx*x )
END FUNCTION


! 3D planewave
FUNCTION pw3d( G, r )
  IMPLICIT NONE
  !
  REAL(8) :: G(3), r(3)
  COMPLEX(8) :: pw3d
  !
  pw3d = exp( (0.d0,1.d0)*dot_product(G,r) )
END FUNCTION


! 1st partial derivative
FUNCTION do_pw3d( G, r, dir )
  IMPLICIT NONE
  !
  REAL(8) :: G(3), r(3)
  COMPLEX(8) :: do_pw3d
  INTEGER :: dir
  !
  do_pw3d = (0.d0,1.d0)*G(dir) * exp( (0.d0,1.d0)*dot_product(G,r) )
END FUNCTION

! 2nd partial derivative
FUNCTION do2_pw3d( G, r, dir )
  IMPLICIT NONE
  !
  REAL(8) :: G(3), r(3)
  COMPLEX(8) :: do2_pw3d
  INTEGER :: dir
  !
  do2_pw3d = -G(dir)**2 * exp( (0.d0,1.d0)*dot_product(G,r) )
  !do2_pw3d = -G(dir)**2 * exp( (0.d0,1.d0)*( G(1)*r(1) + G(2)*r(2) + G(3)*r(3) ) )
END FUNCTION


FUNCTION nabla2_pw3d( G, r )
  IMPLICIT NONE
  !
  REAL(8) :: G(3), r(3)
  COMPLEX(8) :: nabla2_pw3d
  INTEGER :: dir
  COMPLEX(8) :: do2_pw3d
  !
  nabla2_pw3d = -norm2(G)**2 * exp( (0.d0,1.d0)*dot_product(G,r) )
  !nabla2_pw3d = do2_pw3d(G,r,1) + do2_pw3d(G,r,2) + do2_pw3d(G,r,3)
END FUNCTION


!------------------------------------------------------------------------------
PROGRAM t_LF3d
!------------------------------------------------------------------------------
  USE m_LF3d
  USE m_LF1d
  IMPLICIT NONE
  !
  INTEGER :: Nx, Ny, Nz
  REAL(8) :: Lx, Ly, Lz
  TYPE(LF3d_t) :: LF3
  !
  COMPLEX(8), ALLOCATABLE :: coef(:), cx(:), cy(:), cz(:)
  COMPLEX(8), ALLOCATABLE :: nabla2coef(:), d2cx(:), d2cy(:), d2cz(:)
  !
  REAL(8) :: Gx, Gy, Gz, x, y, z
  INTEGER :: ip, i, j, k, ii, jj, kk
  !
  COMPLEX(8) :: pw1d, pw3d, nabla2_pw3d, do_pw3d, do2_pw3d

  !
  ! Initialization on LF
  !
  Nx = 3
  Ny = 3
  Nz = 3
  !
  Lx = 5.d0
  Ly = 5.d0
  Lz = 5.d0
  !
  CALL init_LF3d_p( LF3, (/Nx,Ny,Nz/), (/0.d0,0.d0,0.d0/), (/Lx,Ly,Lz/) )
  !
  CALL info_LF1d( LF3%LFx, .TRUE. )
  CALL info_LF1d( LF3%LFy, .TRUE. )
  CALL info_LF1d( LF3%LFz, .TRUE. )

  ! Initialization of expansion coefficients
  ALLOCATE( coef(Nx*Ny*Nz), nabla2coef(Nx*Ny*Nz) )
  ALLOCATE( cx(Nx), cy(Ny), cz(Nz) )
  ALLOCATE( d2cx(Nx), d2cy(Ny), d2cz(Nz) )
  !
  ! Fill up the coefficients for each direction
  ! FIXME: will not work for G's that (are very close to) integer multiple of 2*PI
  !
  Gx = 1/Lx*2.d0*PI
  Gy = 1/Ly*2.d0*PI
  Gz = 1/Lz*2.d0*PI
  ! FIXME: we omit the factor sqrt(L/N) for the coefficients
  WRITE(*,'(/,1x,A)') 'x direction:'
  DO i = 1,Nx
    x = LF3%LFx%grid(i)
    cx(i) = pw1d( Gx, x )
    WRITE(*,'(1x,I5,3F18.10)') i, x, cx(i)
  ENDDO
  !
  WRITE(*,'(/,1x,A)') 'y direction:'
  DO j = 1,Ny
    y = LF3%LFy%grid(j)
    cy(j) = pw1d( Gy, y )
    WRITE(*,'(1x,I5,3F18.10)') j, y, cy(j)
  ENDDO
  !
  WRITE(*,'(/,1x,A)') 'z direction:'
  DO k = 1,Nz
    z = LF3%LFz%grid(k)
    cz(k) = pw1d( Gz, z )
    WRITE(*,'(1x,I5,3F18.10)') k, z, cz(k)
  ENDDO

  ! Initialize the coefficients for linear grid via lingrid points
  WRITE(*,*) 'Linear grid'
  DO ip = 1,LF3%N
    coef(ip) = pw3d( (/Gx,Gy,Gz/), LF3%lingrid(:,ip) )
    ! For debugging
    i = LF3%lin2xyz( 1, ip )
    j = LF3%lin2xyz( 2, ip )
    k = LF3%lin2xyz( 3, ip )
    WRITE(*,'(1x,I8,4F18.10)') ip, coef(ip), cx(i)*cy(j)*cz(k)
  ENDDO

  ! Calculate Laplacian
  !d2cx = matmul( LF3%LFx%D2jl, cx )
  !d2cy = matmul( LF3%LFy%D2jl, cy )
  !d2cz = matmul( LF3%LFz%D2jl, cz )
  WRITE(*,'(/,1x,A)') 'Calculating Laplacian:'
  DO ip = 1, LF3%N
    i = LF3%lin2xyz( 1, ip )
    j = LF3%lin2xyz( 2, ip )
    k = LF3%lin2xyz( 3, ip )
    !nabla2coef(ip) = d2cx(i)*cy(j)*cz(k) + cx(i)*d2cy(j)*cz(k) + cx(i)*cy(j)*d2cz(k)
    nabla2coef(ip) = (0.d0,0.d0)
    !
    DO ii=1,Nx
      nabla2coef(ip) = nabla2coef(ip) + LF3%LFx%D2jl(ii,i)*coef(LF3%xyz2lin(ii,j,k))
    ENDDO
    !
    DO jj=1,Ny
      nabla2coef(ip) = nabla2coef(ip) + LF3%LFy%D2jl(jj,j)*coef(LF3%xyz2lin(i,jj,k))
    ENDDO
    !
    DO kk=1,Nz
      nabla2coef(ip) = nabla2coef(ip) + LF3%LFz%D2jl(kk,k)*coef(LF3%xyz2lin(i,j,kk))
    ENDDO
    ! compare
    WRITE(*,'(1x,I5,4F18.10)') ip, nabla2coef(ip), nabla2_pw3d( (/Gx,Gy,Gz/), LF3%lingrid(:,ip), 3 )
  ENDDO

  DEALLOCATE( coef, nabla2coef )
  DEALLOCATE( cx, cy, cz )
  DEALLOCATE( d2cx, d2cy, d2cz )

END PROGRAM

