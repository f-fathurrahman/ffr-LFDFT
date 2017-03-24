PROGRAM ex_PETSC
  USE m_constants
  USE m_LF3d, ONLY : D2jl_x => LF3d_D2jl_x, &
                     Npoints => LF3d_Npoints
  IMPLICIT NONE

#include <petsc/finclude/petsc.h>
  !
  INTEGER :: NN(3), Nx, Ny, Nz
  REAL(8) :: AA(3), BB(3)
  INTEGER :: irow, icol
  INTEGER :: i, j, ix, jx, dd
  !
  PetscErrorCode :: ierr
  PetscMPIInt :: rank
  Mat :: nabla2_x

  CALL PetscInitialize( PETSC_NULL_CHARACTER, ierr )

  CALL MPI_Comm_rank( PETSC_COMM_WORLD, rank, ierr )

  NN = (/ 51, 51, 51 /)
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

  !
  CALL MatCreate( PETSC_COMM_WORLD, nabla2_x, ierr )
  CALL MatSetSizes( nabla2_x, PETSC_DECIDE, PETSC_DECIDE, Npoints, Npoints, ierr )
  CALL MatSetUp( nabla2_x, ierr )

  IF ( rank == 0 ) THEN 
    ! Petsc Matrix is row-major !
    DO jx = 1,Nx
      j = (jx-1)*Ny*Nz + 1
      DO ix = 1,Nx
        i = (ix-1)*Ny*Nz + 1
        DO dd = 1,Ny*Nz
          irow = i + dd - 1 - 1 ! convert to 0-based index
          icol = j + dd - 1 - 1
          WRITE(*,'(5I8)') ix, jx, dd, irow, icol
          !CALL MatSetValues( nabla2_x, 1, irow, 1, icol, D2jl_x(ix,jx), INSERT_VALUES, ierr )
          CALL MatSetValue( nabla2_x, irow, icol, D2jl_x(ix,jx), INSERT_VALUES, ierr )
        ENDDO 
      ENDDO
    ENDDO 
  ENDIF 

  CALL MatAssemblyBegin( nabla2_x, MAT_FINAL_ASSEMBLY, ierr )
  CALL MatAssemblyEnd( nabla2_x, MAT_FINAL_ASSEMBLY, ierr )

  !CALL MatView( nabla2_x, PETSC_VIEWER_STDOUT_WORLD, ierr )

  ! Free memory
  CALL dealloc_LF3d()
  ! Shutdown PETSC
  CALL PetscFinalize( ierr )
  
END PROGRAM

