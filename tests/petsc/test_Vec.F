PROGRAM test_Vec

  IMPLICIT NONE 

#include <petsc/finclude/petsc.h>

  INTEGER :: Nstates
  Vec, ALLOCATABLE :: vv(:)
  Vec :: x, y, w, z(5)
  PetscMPIInt :: ierr
  PetscInt :: NN(3), Nx, Ny, Nz
  PetscInt :: Npoints
  PetscScalar :: ONE, TWO, THREE
  CHARACTER(40) :: name_
  PetscMPIInt :: rank
  INTEGER :: ii

  Nstates = 5

  ALLOCATE( vv(Nstates) )

  NN = (/ 5, 3, 3 /)
  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)
  Npoints = Nx * Ny * Nz

  ONE = 1.0
  TWO = 2.0
  THREE = 3.0

  CALL PetscInitialize( PETSC_NULL_CHARACTER, ierr )

  CALL MPI_Comm_rank( PETSC_COMM_WORLD, rank, ierr )

  ! Initialize vectors
  CALL VecCreate( PETSC_COMM_WORLD, x, ierr )
  CALL VecSetSizes( x, PETSC_DECIDE, Npoints, ierr )
  CALL VecSetFromOptions( x, ierr )
  CALL VecGetType( x, name_, ierr )

  CALL VecSet( x, ONE, ierr )

  CALL VecDuplicate( x, y, ierr )
  CALL VecDuplicate( x, w, ierr )
  CALL VecDuplicateVecs( x, 5, z, ierr )
  CALL VecDuplicateVecs( x, Nstates, vv, ierr )

  CALL VecSet( vv(1), TWO, ierr )
  CALL VecSet( vv(2), 1.10125d0, ierr )

  DO ii = 1, Npoints
    CALL VecSetValue( y, ii-1, 0.5d0*ii, INSERT_VALUES, ierr )
    CALL VecSetValue( w, ii-1, 0.4d0*(ii-1), INSERT_VALUES, ierr )
  ENDDO

  CALL VecAssemblyBegin( y, ierr )
  CALL VecAssemblyEnd( y, ierr )

  CALL VecAssemblyBegin( w, ierr )
  CALL VecAssemblyEnd( w, ierr )

  WRITE(*,*) 'name_ = ', name_

  CALL PetscViewerPushFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INDEX, ierr )
  !CALL PetscViewerPushFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO_DETAIL, ierr )

  !CALL VecView( z(1), PETSC_VIEWER_STDOUT_WORLD, ierr )
  !CALL VecView( y, PETSC_VIEWER_STDOUT_WORLD, ierr )
  !CALL VecView( w, PETSC_VIEWER_STDOUT_WORLD, ierr )
  !CALL VecView( vv(1), PETSC_VIEWER_STDOUT_WORLD, ierr )
  CALL VecView( vv(2), PETSC_VIEWER_STDOUT_WORLD, ierr )

  CALL PetscFinalize( ierr )

END PROGRAM 

