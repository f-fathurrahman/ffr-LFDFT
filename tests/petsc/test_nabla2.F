#include "add_nabla2_x.F"
#include "add_nabla2_y.F"

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
  INTEGER :: rowGbl
  INTEGER :: ix, irow
  INTEGER :: nnzc
  INTEGER, ALLOCATABLE :: nnz(:)
  !
  PetscErrorCode :: ierr
  PetscMPIInt :: rank
  Mat :: nabla2

  CALL PetscInitialize( PETSC_NULL_CHARACTER, ierr )

  CALL MPI_Comm_rank( PETSC_COMM_WORLD, rank, ierr )

  !NN = (/ 101, 101, 101 /)
  NN = (/ 81, 81, 81 /)
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

  !
  CALL MatCreateSeqAIJ( PETSC_COMM_SELF, Npoints, Npoints, nnzc, PETSC_NULL_INTEGER, nabla2, ierr )
  !CALL MatCreateSeqAIJ( PETSC_COMM_SELF, Npoints, Npoints, nnzc, nnz, nabla2, ierr )
  CALL MatSetOption( nabla2, MAT_SYMMETRIC, PETSC_TRUE, ierr )
  CALL MatSetUp( nabla2, ierr )


  IF ( rank == 0 ) THEN

    CALL add_nabla2_x( Nx, Ny, Nz, D2jl_x, nabla2 )

    CALL add_nabla2_y( Nx, Ny, Nz, D2jl_y, nabla2 )

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

