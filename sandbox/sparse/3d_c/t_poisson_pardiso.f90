! efefer 1 January 2016

MODULE gbl_poisson
  USE m_LF3d
  IMPLICIT NONE
  INTEGER :: Nx, Ny, Nz
  REAL(8) :: Lx, Ly, Lz
  TYPE(LF3d_t) :: LF
END MODULE


MODULE gbl_laplacian
  TYPE SparseMat
    REAL(8), ALLOCATABLE :: values(:)
    INTEGER, ALLOCATABLE :: column(:)
    INTEGER, ALLOCATABLE :: rowIdx(:)
  END TYPE

  TYPE(SparseMat) :: Laplacian
END MODULE


SUBROUTINE init_Laplacian()
  USE gbl_poisson, ONLY : LF, Nx
  USE gbl_laplacian, ONLY : Laplacian
  IMPLICIT NONE
  INTEGER :: istep, ip1, ip2, i1,i2, j1,j2, k1,k2
  INTEGER :: N, NNZ
  REAL(8) :: nabla2
  LOGICAL :: Tnz(3)

  ! SPECIAL CASE !!! Nx = Ny = Nz
  N = Nx

  ! Nabla2 matrix
  DO istep = 1, 2
    IF( istep == 2 ) THEN
      ALLOCATE( Laplacian%values(NNZ) )
      ALLOCATE( Laplacian%column(NNZ) )
      ALLOCATE( Laplacian%rowIdx(N**3+1) )
      Laplacian%rowIdx(1) = 1 ! heuristic?
    ENDIF
    NNZ = 0
    DO ip2 = 1, N**3
      DO ip1 = ip2, N**3
        Tnz(:) = .FALSE.
        i1 = LF%lin2xyz(1,ip1)
        i2 = LF%lin2xyz(1,ip2)
        !
        j1 = LF%lin2xyz(2,ip1)
        j2 = LF%lin2xyz(2,ip2)
        !
        k1 = LF%lin2xyz(3,ip1)
        k2 = LF%lin2xyz(3,ip2)
        !
        nabla2 = 0.d0
        !
        IF( j1 == j2 .AND. k1 == k2 ) THEN
          nabla2 = nabla2 + LF%LFx%D2jl(i1,i2)
          Tnz(1) = .TRUE.
        ENDIF
        !
        IF( i1 == i2 .AND. k1 == k2 ) THEN
          nabla2 = nabla2 + LF%LFy%D2jl(j1,j2)
          Tnz(2) = .TRUE.
        ENDIF
        !
        IF( i1 == i2 .AND. j1 == j2 ) THEN
          nabla2 = nabla2 + LF%LFz%D2jl(k1,k2)
          Tnz(3) = .TRUE.
        ENDIF
        !WRITE(*,*) ip1,ip2, any(Tnz)
        IF( ANY(Tnz) ) THEN
          NNZ = NNZ + 1
          IF( istep == 2 ) THEN
            Laplacian%values(NNZ) = nabla2
            Laplacian%column(NNZ) = ip1
          ENDIF
        ENDIF
        ! 
      ENDDO
      IF( istep == 2 ) Laplacian%rowIdx(ip2+1) = NNZ + 1 
    ENDDO
  ENDDO ! istep

  WRITE(*,*) 'NNZ, N**3 = ', NNZ, N**3, NNZ - N**3, dble(N**6)
  WRITE(*,*) 'Percentage = ', ( (NNZ - N**3)*2.d0 + N**3 ) / dble(N**6) * 100d0
END SUBROUTINE


PROGRAM t_poisson
  USE m_constants
  USE m_LF3d
  USE gbl_poisson
  IMPLICIT NONE
  !

  REAL(8), ALLOCATABLE :: rho(:) ! density
  REAL(8), ALLOCATABLE :: phi(:) ! potential
  !
  REAL(8) :: sigma1, sigma2, r, x0, y0, z0, deltaV
  INTEGER :: ip
  REAL(8) :: Uana, Unum, t1, t2

  Nx = 35
  Ny = 35
  Nz = 35
  !
  Lx = 16.d0
  Ly = 16.d0
  Lz = 16.d0
  !
  CALL init_LF3d_c( LF, (/Nx,Ny,Nz/), (/0.d0,0.d0,0.d0/), (/Lx,Ly,Lz/) )

  CALL cpu_time(t1)
  CALL init_Laplacian()
  CALL cpu_time(t2)
  WRITE(*,'(1x,A,F18.10)') 'Allocation of Laplacian = ', t2-t1
  
  ALLOCATE( rho(Nx*Ny*Nz) )
  ALLOCATE( phi(Nx*Ny*Nz) )

  ! center of the box
  x0 = Lx/2.d0
  y0 = Ly/2.d0
  z0 = Lz/2.d0
  ! Initialize
  sigma1 = 0.75d0
  sigma2 = 0.50d0
  DO ip = 1, LF%N
    r = norm2( LF%lingrid(:,ip) - (/x0,y0,z0/) )
    rho(ip) = exp(-r**2/(2*sigma2**2))/(2*pi*sigma2**2)**1.5d0 - &
              exp(-r**2/(2*sigma1**2))/(2*pi*sigma1**2)**1.5d0
  ENDDO

  ! For cluster LF
  deltaV = LF%LFx%h * LF%LFy%h * LF%LFz%h
  WRITE(*,*) 'deltaV = ', deltaV
  WRITE(*,*) sum( rho(:) )*deltaV

  ! Solve Poisson equation
  CALL solve_poisson_pardiso( -4.d0*PI*rho, phi )

  !
  Unum = 0.5d0*sum( rho(:)*phi(:) )*deltaV
  Uana = ( (1.d0/sigma1 + 1.d0/sigma2)/2.d0 - sqrt(2.d0)/sqrt(sigma1**2 + sigma2**2) )/sqrt(PI)
  WRITE(*,*) 'Unum = ', Unum
  WRITE(*,*) 'Uana = ', Uana

  DEALLOCATE( rho, phi )

END PROGRAM

