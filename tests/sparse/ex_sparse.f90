PROGRAM ex_sparse

  USE m_constants
  USE m_nabla2_sparse, ONLY : nzval => nabla2_nzval, &
                              colptr => nabla2_colptr, &
                              rowval => nabla2_rowval

  IMPLICIT NONE
  !
  INTEGER :: NN(3), Npoints
  REAL(8) :: AA(3), BB(3)
  REAL(8), ALLOCATABLE :: v(:), Av(:)
  INTEGER :: ip
  INTEGER, ALLOCATABLE :: iwork(:)
  INTEGER :: nwork
  REAL(8) :: t1, t2
  REAL(8), ALLOCATABLE :: alu_ilu0(:)
  INTEGER, ALLOCATABLE :: jlu_ilu0(:), ju_ilu0(:), iw_ilu0(:)
  INTEGER :: ierr
  INTEGER :: Nx, Ny, Nz

  !NN = (/ 3, 2, 4 /)
  NN = (/ 64, 64, 64 /)
  !NN = (/ 80, 80, 80 /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)

  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)
  Npoints = Nx * Ny * Nz

  CALL init_LF3d_c( NN, AA, BB )
  CALL info_LF3d()

  CALL cpu_time(t1)
  CALL init_nabla2_sparse()
  CALL cpu_time(t2)
  WRITE(*,*) 'init_nabla2_sparse = ', t2-t1

  ! sort
  CALL cpu_time(t1)
  nwork = max( Npoints+1, 2*(colptr(Npoints+1)-colptr(1)) )
  WRITE(*,*) 'Npoints+1 = ', Npoints+1
  WRITE(*,*) '2*nnz = ', 2*(colptr(Npoints+1)-colptr(1)) 
  WRITE(*,*) 'nwork = ', nwork
  ALLOCATE( iwork( nwork ) )
  CALL csort( Npoints, nzval, rowval, colptr, iwork, .TRUE. )
  DEALLOCATE( iwork )
  CALL cpu_time(t2)
  WRITE(*,*) 'csort: ', t2-t1

  !DO ip = 1, size(nzval)
  !  WRITE(*,'(1x,2I5,F18.10)') ip, rowval(ip), nzval(ip)
  !ENDDO 
  
  ALLOCATE( v(Npoints) )
  ALLOCATE( Av(Npoints) )

  v(:) = 1.1d0

  ! use SPARSKIT amux for matrix multiplication
  CALL cpu_time(t1)
  CALL amux( Npoints, v, Av, nzval, rowval, colptr )
  CALL cpu_time(t2)
  WRITE(*,*) 'amux:', t2-t1

  !DO ip = 1, Npoints
  !  WRITE(*,'(1x,I5,2F18.10)') ip, v(ip), Av(ip)
  !ENDDO 

  ! Build ILU0 preconditioner
  CALL cpu_time(t1)
  ALLOCATE( alu_ilu0(Npoints*(Nx+Ny+Nz-2)) )
  ALLOCATE( jlu_ilu0(Npoints*(Nx+Ny+Nz-2)) )
  ALLOCATE( ju_ilu0(Npoints) )
  ALLOCATE( iw_ilu0(Npoints) )
  !
  CALL ilu0( Npoints, nzval, rowval, colptr, alu_ilu0, jlu_ilu0, ju_ilu0, iw_ilu0, ierr )
  !WRITE(*,*) 'ierr = ', ierr
  CALL cpu_time(t2)
  WRITE(*,*) 'ilu0 = ', t2-t1

  ! Call preconditioner
  WRITE(*,*) 'Before prec: sum(v) = ', sum(v)
  CALL cpu_time(t1)
  CALL lusol( Npoints, v, v, alu_ilu0, jlu_ilu0, ju_ilu0 )
  CALL cpu_time(t2)
  WRITE(*,*) 'lusol = ', t2-t1
  WRITE(*,*) 'After prec: sum(v) = ', sum(v)

  DEALLOCATE( v, Av )
  DEALLOCATE( alu_ilu0, jlu_ilu0, ju_ilu0, iw_ilu0 )

  CALL dealloc_nabla2_sparse()
  CALL dealloc_LF3d()
  
  WRITE(*,*) 'Program ended normally'
END PROGRAM

