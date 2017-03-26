include "build_nabla2_x.F"

PROGRAM ex_PETSC
  USE m_constants
  USE m_LF3d, ONLY : D2jl_x => LF3d_D2jl_x, &
                     Npoints => LF3d_Npoints
  IMPLICIT NONE

#include <petsc/finclude/petsc.h>
  !
  INTEGER :: NN(3), Nx, Ny, Nz
  REAL(8) :: AA(3), BB(3)
  !
  PetscErrorCode :: ierr
  PetscMPIInt :: rank
  Mat :: nabla2_x

  CALL PetscInitialize( PETSC_NULL_CHARACTER, ierr )

  CALL MPI_Comm_rank( PETSC_COMM_WORLD, rank, ierr )

  NN = (/ 5, 5, 5 /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)

  ! Initialize LF3d
  CALL init_LF3d_p( NN, AA, BB )
  IF ( rank == 0 ) THEN 
    CALL info_LF3d()
  ENDIF 

  CALL build_nabla2_x( NN, D2jl_x, nabla2_x ) 

  CALL MatView( nabla2_x, PETSC_VIEWER_STDOUT_WORLD, ierr )

  ! Free memory
  CALL dealloc_LF3d()
  ! Shutdown PETSC
  CALL PetscFinalize( ierr )
  
END PROGRAM

