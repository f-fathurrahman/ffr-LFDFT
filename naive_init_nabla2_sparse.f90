
! A naive subroutine to initialize Laplacian matrix
! TODO Speed up the calculation
! 
! Alternative: use kronecker product
SUBROUTINE init_nabla2_sparse()
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : lin2xyz => LF3d_lin2xyz, &
                     D2jl_x => LF3d_D2jl_x, &
                     D2jl_y => LF3d_D2jl_y, &
                     D2jl_z => LF3d_D2jl_z, &
                     Npoints => LF3d_Npoints
  USE m_nabla2_sparse, ONLY : NNZ => nabla2_NNZ, &
                              values => nabla2_values, &
                              column => nabla2_column, &
                              rowIdx => nabla2_rowIdx
  IMPLICIT NONE
  !
  REAL(8) :: nabla2
  LOGICAL :: Tnz(3)
  INTEGER :: ip1, ip2, i1,i2,j1,j2,k1,k2, istep
  
  WRITE(*,'(/,1x,A)') 'Initializing nabla2_sparse'

  ! Nabla2 matrix
  DO istep = 1, 2
    IF( istep == 2 ) THEN
      ALLOCATE( values( NNZ ) )
      ALLOCATE( column( NNZ ) )
      ALLOCATE( rowIdx( Npoints+1 ) )
      rowIdx(1) = 1 !
    ENDIF
    NNZ = 0
    DO ip2 = 1, Npoints
      DO ip1 = ip2, Npoints
        Tnz(:) = .FALSE.
        i1 = lin2xyz(1,ip1)
        i2 = lin2xyz(1,ip2)
        !
        j1 = lin2xyz(2,ip1)
        j2 = lin2xyz(2,ip2)
        !
        k1 = lin2xyz(3,ip1)
        k2 = lin2xyz(3,ip2)
        !
        nabla2 = 0.d0
        !
        IF( j1 == j2 .AND. k1 == k2 ) THEN
          nabla2 = nabla2 + D2jl_x(i1,i2)
          Tnz(1) = .TRUE.
        ENDIF
        !
        IF( i1 == i2 .AND. k1 == k2 ) THEN
          nabla2 = nabla2 + D2jl_y(j1,j2)
          Tnz(2) = .TRUE.
        ENDIF
        !
        IF( i1 == i2 .AND. j1 == j2 ) THEN
          nabla2 = nabla2 + D2jl_z(k1,k2)
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

  WRITE(*,*) 'NNZ = ', NNZ
  WRITE(*,*) 'Percentage = ', ( (NNz - Npoints)*2.d0 + Npoints ) / Npoints**2 * 100d0
  WRITE(*,*) 'Done initializing nabla2_sparse'

END SUBROUTINE


