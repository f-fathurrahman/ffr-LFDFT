PROGRAM ex_PETSC
  USE m_constants
  USE m_LF3d, ONLY : D2jl_x => LF3d_D2jl_x, &
                     D2jl_y => LF3d_D2jl_y, &
                     D2jl_z => LF3d_D2jl_z, &
                     Npoints => LF3d_Npoints
  IMPLICIT NONE

#include <petsc/finclude/petsc.h>
  !
  INTEGER :: NN(3), Nx, Ny, Nz
  REAL(8) :: AA(3), BB(3)
  INTEGER :: nnzc
  INTEGER, ALLOCATABLE :: nnz(:)
  !
  PetscErrorCode :: ierr
  PetscMPIInt :: rank
  Mat :: nabla2

  CALL PetscInitialize( PETSC_NULL_CHARACTER, ierr )

  CALL MPI_Comm_rank( PETSC_COMM_WORLD, rank, ierr )

  !NN = (/ 101, 101, 101 /)
  NN = (/ 63, 63, 63 /)
  !NN = (/ 3, 3, 5 /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)

  ! Shortcuts
  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  ! Initialize LF3d
  CALL init_LF3d_p( NN, AA, BB )
  IF ( rank == 0 ) THEN 
    CALL info_LF3d()
  ENDIF 

  ! Number of nonzeros per each row
  nnzc = Nx + Ny + Nz - 2
  ALLOCATE( nnz(Npoints) )
  nnz(:) = nnzc

  WRITE(*,*) 'Number of nonzeros = ', nnzc*Npoints
  WRITE(*,*) 'Memory required to those nonzeros is ', nnzc*Npoints*8.d0/1024.d0/1024.d0/1024.d0, ' GB'

  !
  CALL MatCreateSeqAIJ( PETSC_COMM_SELF, Npoints, Npoints, nnzc, PETSC_NULL_INTEGER, nabla2, ierr )
  !CALL MatCreateSeqAIJ( PETSC_COMM_SELF, Npoints, Npoints, nnzc, nnz, nabla2, ierr )
  !CALL MatSetOption( nabla2, MAT_SYMMETRIC, PETSC_TRUE, ierr )
  CALL MatSetUp( nabla2, ierr )

  IF ( rank == 0 ) THEN

    CALL add_nabla2_z( Nx, Ny, Nz, D2jl_z, nabla2 )
    CALL add_nabla2_y( Nx, Ny, Nz, D2jl_y, nabla2 )
    CALL add_nabla2_x( Nx, Ny, Nz, D2jl_x, nabla2 )

  ENDIF ! rank == 0

  CALL MatAssemblyBegin( nabla2, MAT_FINAL_ASSEMBLY, ierr )
  CALL MatAssemblyEnd( nabla2, MAT_FINAL_ASSEMBLY, ierr )

  !CALL MatView( nabla2, PETSC_VIEWER_STDOUT_WORLD, ierr )

  ! Free memory
  DEALLOCATE( nnz )
  CALL dealloc_LF3d()
  ! Shutdown PETSC
  CALL PetscFinalize( ierr )
  
END PROGRAM

