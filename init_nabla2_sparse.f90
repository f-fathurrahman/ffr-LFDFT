! A naive subroutine to initialize Laplacian matrix
! TODO Speed up the calculation
! 
! Alternative: use kronecker product
SUBROUTINE init_nabla2_sparse( t_read )
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
  ! argument
  LOGICAL :: t_read
  !
  REAL(8) :: nabla2
  LOGICAL :: Tnz(3)
  INTEGER :: ip1, ip2, i1,i2,j1,j2,k1,k2, ii, ip
  !
  INTEGER, PARAMETER :: IU_values=111
  INTEGER, PARAMETER :: IU_column=112
  INTEGER, PARAMETER :: IU_rowIdx=113
  INTEGER, PARAMETER :: IU_nnz=110
  
  WRITE(*,'(/,1x,A)') 'Initializing nabla2_sparse'
 
  IF( .NOT. t_read ) THEN

    OPEN( IU_values, file='values.dat', form='unformatted', status='replace', access='stream' )
    OPEN( IU_column, file='column.dat', form='unformatted', status='replace', access='stream' )
    OPEN( IU_rowIdx, file='rowIdx.dat', form='unformatted', status='replace', access='stream' )

    OPEN( IU_NNZ, file='NNZ.dat', form='formatted', status='replace' )
  
    WRITE( IU_rowIdx ) 1
  
    ! Nabla2 matrix
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
          WRITE(IU_values) nabla2
          WRITE(IU_column) ip1
        ENDIF
        ! 
      ENDDO
      WRITE(IU_rowIdx) NNz + 1
    ENDDO
  
    WRITE( IU_NNZ, * ) NNZ
    
    CLOSE( IU_values )
    CLOSE( IU_column )
    CLOSE( IU_rowIdx )
    CLOSE( IU_NNZ )

  ENDIF ! t_read = .FALSE.

  IF( t_read ) THEN
    OPEN( IU_NNZ, file='NNZ.dat', form='formatted', status='old' )
    READ( IU_NNZ, * ) NNZ
    CLOSE( IU_NNZ )
  ENDIF

  WRITE(*,*) 'NNZ = ', NNZ
  WRITE(*,*) 'Npoints = ', Npoints
  WRITE(*,*) 'Npoints**2 = ', dble(Npoints)**2
  WRITE(*,*) 'Percentage = ', ( (NNZ - Npoints)*2.d0 + Npoints ) / dble(Npoints)**2 * 100d0

  ALLOCATE( values( NNZ ) )
  ALLOCATE( column( NNZ ) )
  ALLOCATE( rowIdx( Npoints+1 ) )

  ! Open again
  OPEN( IU_values, file='values.dat', form='unformatted', status='old', access='stream' )
  OPEN( IU_column, file='column.dat', form='unformatted', status='old', access='stream' )
  OPEN( IU_rowIdx, file='rowIdx.dat', form='unformatted', status='old', access='stream' )

  DO ii = 1, NNZ
    READ( IU_values ) values(ii)
    READ( IU_column ) column(ii)
  ENDDO

  DO ip = 1, Npoints+1
    READ( IU_rowIdx ) rowIdx(ip)
  ENDDO

  CLOSE( IU_values )
  CLOSE( IU_column )
  CLOSE( IU_rowIdx )

  WRITE(*,*) 'Done initializing nabla2_sparse'

END SUBROUTINE

