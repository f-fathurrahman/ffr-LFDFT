!------------------------------------------------------------------------------
PROGRAM nabla2_sparse
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF3d
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: N = 8
  REAL(8), PARAMETER :: L = 2.d0
  !
  TYPE(LF3d_t) :: LF
  REAL(8) :: nabla2
  INTEGER :: NNZ
  LOGICAL :: Tnz(3)
  INTEGER :: ip1, ip2, i1,i2,j1,j2,k1,k2, istep
  REAL(8), ALLOCATABLE :: values(:)
  INTEGER, ALLOCATABLE :: column(:)
  INTEGER, ALLOCATABLE :: rowIdx(:)
  
  ! Initialize the basis functions
  CALL init_LF3d_c( LF, (/N,N,N/), (/0.d0,0.d0,0.d0/), (/L,L,L/) )

  ! Nabla2 matrix
  DO istep = 1, 2
    IF( istep == 2 ) THEN
      ALLOCATE( values(NNZ) )
      ALLOCATE( column(NNZ) )
      ALLOCATE( rowIdx(N**3+1) )
      rowIdx(1) = 1 ! heuristic?
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
          NNz = NNz + 1
          IF( istep == 2 ) THEN
            values(NNz) = nabla2
            column(NNz) = ip1
          ENDIF
        ENDIF
        ! 
      ENDDO
      IF( istep == 2 ) rowIdx(ip2+1) = NNz + 1 
    ENDDO
  ENDDO ! istep

  WRITE(*,*) 'NNz = ', NNz
  WRITE(*,*) 'Percentage = ', ( (NNz - N**3)*2.d0 + N**3 ) / N**6 * 100d0

  DO ip1 = 1, NNz
    WRITE(100,*) values(ip1)
    WRITE(200,*) column(ip1)
  ENDDO
  !
  DO ip1 = 1, N**3+1
    WRITE(300,*) rowIdx(ip1)
  ENDDO

END PROGRAM


