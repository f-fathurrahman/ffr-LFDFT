PROGRAM ex_PETSC
  USE m_constants
  USE m_LF3d, ONLY : D2jl_x => LF3d_D2jl_x, &
                     Npoints => LF3d_Npoints
  IMPLICIT NONE

#include <petsc/finclude/petsc.h>
  !
  INTEGER :: NN(3), Nx, Ny, Nz
  REAL(8) :: AA(3), BB(3)
  INTEGER :: rowGbl
  INTEGER :: ix, irow
  INTEGER, ALLOCATABLE :: colGbl(:), colGbl_orig(:)
  !
  PetscErrorCode :: ierr
  PetscMPIInt :: rank
  Mat :: nabla2_x

  CALL PetscInitialize( PETSC_NULL_CHARACTER, ierr )

  CALL MPI_Comm_rank( PETSC_COMM_WORLD, rank, ierr )

  NN = (/ 101, 101, 101 /)
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
  !CALL MatCreate( PETSC_COMM_WORLD, nabla2_x, ierr )
  !CALL MatSetSizes( nabla2_x, PETSC_DECIDE, PETSC_DECIDE, Npoints, Npoints, ierr )
  CALL MatCreateSeqAIJ( PETSC_COMM_SELF, Npoints, Npoints, Nx, PETSC_NULL_INTEGER, nabla2_x, ierr )
  CALL MatSetUp( nabla2_x, ierr )


  IF ( rank == 0 ) THEN
    ALLOCATE( colGbl(Nx) )
    ALLOCATE( colGbl_orig(Nx) )
    !
    colGbl(1) = 1
    DO ix = 2,Nx
      colGbl(ix) = colGbl(1) + Ny*Nz*(ix-1)
    ENDDO
    ! convert to zero-based indexing
    !colGbl(:) = colGbl(:) - 1
    !WRITE(*,*) 'initial colGbl: ', colGbl(:)
    !WRITE(*,*)
    colGbl_orig(:) = colGbl(:)
    !
    DO ix = 1,Nx
      !WRITE(*,*)
      DO irow = 1,Ny*Nz
        colGbl(:) = colGbl_orig(:) + irow - 1
        rowGbl = irow + (ix-1)*Ny*Nz
        !WRITE(*,*) 'rowGbl, colGbl ', rowGbl, ', ', colGbl(:)
        CALL MatSetValues( nabla2_x, 1, rowGbl-1, Nx, colGbl-1, D2jl_x(1:Nx,ix), INSERT_VALUES, ierr )
        ! I transposed idx for D2jl_x
      ENDDO
    ENDDO
    DEALLOCATE( colGbl )
    DEALLOCATE( colGbl_orig )
  ENDIF ! rank == 0


  CALL MatAssemblyBegin( nabla2_x, MAT_FINAL_ASSEMBLY, ierr )
  CALL MatAssemblyEnd( nabla2_x, MAT_FINAL_ASSEMBLY, ierr )

  !CALL MatView( nabla2_x, PETSC_VIEWER_STDOUT_WORLD, ierr )

  ! Free memory
  CALL dealloc_LF3d()
  ! Shutdown PETSC
  CALL PetscFinalize( ierr )
  
END PROGRAM

