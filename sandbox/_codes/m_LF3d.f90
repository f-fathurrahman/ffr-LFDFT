! efefer, 28-Dec-2015

MODULE m_LF3d

  USE m_LF1d
  IMPLICIT NONE

  TYPE LF3d_t
    TYPE(LF1d_t) :: LFx, LFy, LFz
    INTEGER :: Nx, Ny, Nz, N
    REAL(8) :: Lx, Ly, Lz
    REAL(8), ALLOCATABLE :: lingrid(:,:)
    INTEGER, ALLOCATABLE :: xyz2lin(:,:,:)
    INTEGER, ALLOCATABLE :: lin2xyz(:,:)
  END TYPE

  CONTAINS

  !--------------------------------------
  SUBROUTINE init_LF3d_p( LF3d, N, A, B )
  !--------------------------------------
    IMPLICIT NONE
    !
    TYPE(LF3d_t) :: LF3d
    INTEGER :: N(3)
    REAL(8) :: A(3), B(3)
    !
    INTEGER :: Nx, Ny, Nz
    REAL(8) :: Lx, Ly, Lz
    !
    INTEGER :: i, j, k, ip
    !
    Nx = N(1)
    Ny = N(2)
    Nz = N(3)
    !
    Lx = B(1) - A(1)
    Ly = B(2) - A(2)
    Lz = B(3) - A(3)
    !
    CALL init_LF1d_p( LF3d%LFx, Nx, A(1), B(1) )
    CALL init_LF1d_p( LF3d%LFy, Ny, A(2), B(2) )
    CALL init_LF1d_p( LF3d%LFz, Nz, A(3), B(3) )
    !
    LF3d%Nx = Nx
    LF3d%Ny = Ny
    LF3d%Nz = Nz
    LF3d%N  = Nx*Ny*Nz
    !
    LF3d%Lx = Lx
    LF3d%Ly = Ly
    LF3d%Lz = Lz
    !
    ALLOCATE( LF3d%lingrid( 3, Nx*Ny*Nz ) )
    ALLOCATE( LF3d%xyz2lin( Nx, Ny, Nz ) )
    ALLOCATE( LF3d%lin2xyz( 3, Nx*Ny*Nz ) )
    ip = 0
    DO k = 1, Nz
      DO j = 1, Ny
        DO i = 1, Nx
          ip = ip + 1
          LF3d%lingrid( 1, ip ) = LF3d%LFx%grid(i)
          LF3d%lingrid( 2, ip ) = LF3d%LFy%grid(j)
          LF3d%lingrid( 3, ip ) = LF3d%LFz%grid(k)
          !
          LF3d%xyz2lin( i, j, k ) = ip
          LF3d%lin2xyz( 1:3, ip ) = (/ i, j, k /)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE


  SUBROUTINE init_LF3d_c( LF3d, N, A, B )
    IMPLICIT NONE
    !
    TYPE(LF3d_t) :: LF3d
    INTEGER :: N(3)
    REAL(8) :: A(3), B(3)
    !
    INTEGER :: Nx, Ny, Nz
    REAL(8) :: Lx, Ly, Lz
    INTEGER :: i, j, k, ip
    !
    Nx = N(1)
    Ny = N(2)
    Nz = N(3)
    !
    Lx = B(1) - A(1)
    Ly = B(2) - A(2)
    Lz = B(3) - A(3)
    !
    CALL init_LF1d_c( LF3d%LFx, Nx, A(1), B(1) )
    CALL init_LF1d_c( LF3d%LFy, Ny, A(2), B(2) )
    CALL init_LF1d_c( LF3d%LFz, Nz, A(3), B(3) )
    !
    LF3d%Nx = Nx
    LF3d%Ny = Ny
    LF3d%Nz = Nz
    LF3d%N  = Nx*Ny*Nz
    !
    LF3d%Lx = Lx
    LF3d%Ly = Ly
    LF3d%Lz = Lz
    !
    ALLOCATE( LF3d%lingrid( 3, Nx*Ny*Nz ) )
    ALLOCATE( LF3d%xyz2lin( Nx, Ny, Nz ) )
    ALLOCATE( LF3d%lin2xyz( 3, Nx*Ny*Nz ) )
    ip = 0
    DO k = 1, Nz
      DO j = 1, Ny
        DO i = 1, Nx
          ip = ip + 1
          LF3d%lingrid( 1, ip ) = LF3d%LFx%grid(i)
          LF3d%lingrid( 2, ip ) = LF3d%LFy%grid(j)
          LF3d%lingrid( 3, ip ) = LF3d%LFz%grid(k)
          !
          LF3d%xyz2lin( i, j, k ) = ip
          LF3d%lin2xyz( 1:3, ip ) = (/ i, j, k /)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE init_LF3d_c



  SUBROUTINE init_LF3d_sinc( LF3d, N, h )
    IMPLICIT NONE
    !
    TYPE(LF3d_t) :: LF3d
    INTEGER :: N(3)
    REAL(8) :: h(3)
    !
    INTEGER :: Nx, Ny, Nz
    REAL(8) :: Lx, Ly, Lz
    INTEGER :: i, j, k, ip
    !
    Nx = N(1)
    Ny = N(2)
    Nz = N(3)
    !
    CALL init_LF1d_sinc( LF3d%LFx, Nx, h(1) )
    CALL init_LF1d_sinc( LF3d%LFy, Ny, h(2) )
    CALL init_LF1d_sinc( LF3d%LFz, Nz, h(3) )
    !
    Lx = LF3d%LFx%L
    Ly = LF3d%LFy%L
    Lz = LF3d%LFz%L
    !
    LF3d%Nx = Nx
    LF3d%Ny = Ny
    LF3d%Nz = Nz
    LF3d%N  = Nx*Ny*Nz
    !
    LF3d%Lx = Lx
    LF3d%Ly = Ly
    LF3d%Lz = Lz
    !
    ALLOCATE( LF3d%lingrid( 3, Nx*Ny*Nz ) )
    ALLOCATE( LF3d%xyz2lin( Nx, Ny, Nz ) )
    ALLOCATE( LF3d%lin2xyz( 3, Nx*Ny*Nz ) )
    ip = 0
    DO k = 1, Nz
      DO j = 1, Ny
        DO i = 1, Nx
          ip = ip + 1
          LF3d%lingrid( 1, ip ) = LF3d%LFx%grid(i)
          LF3d%lingrid( 2, ip ) = LF3d%LFy%grid(j)
          LF3d%lingrid( 3, ip ) = LF3d%LFz%grid(k)
          !
          LF3d%xyz2lin( i, j, k ) = ip
          LF3d%lin2xyz( 1:3, ip ) = (/ i, j, k /)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE init_LF3d_sinc


END MODULE

