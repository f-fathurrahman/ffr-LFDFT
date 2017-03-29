! Should be called only by master process
SUBROUTINE add_nabla2_z( Nx, Ny, Nz, D2jl_z, nabla2 )

  IMPLICIT NONE  

#include <petsc/finclude/petsc.h>
  
  INTEGER :: Nx, Ny, Nz
  REAL(8) :: D2jl_z(Nz,Nz)
  Mat :: nabla2
  !
  INTEGER :: irow, ii
  PetscInt, ALLOCATABLE :: colGbl(:), colGbl_orig(:)
  PetscInt :: rowGbl, rowLoc

  PetscErrorCode :: ierr

  WRITE(*,*) 'Adding nabla2_z'

  ALLOCATE( colGbl(Nz) )
  ALLOCATE( colGbl_orig(Nz) )
  
  ! Initialize pattern for column indices
  DO ii = 1,Nz
    colGbl_orig(ii) = ii
  ENDDO
  colGbl(:) = colGbl_orig(:)

  !WRITE(*,*)
  !WRITE(*,*) 'colGbl_orig = ', colGbl_orig

  rowGbl = 0
  DO irow = 1,Nx*Ny
    DO rowLoc = 1,Nz
      rowGbl = rowGbl + 1
      colGbl(:) = colGbl_orig(:) + (irow-1)*Nz
      !WRITE(*,*) rowGbl, ':', colGbl(:)
      CALL MatSetValues( nabla2, 1, rowGbl-1, Nz, colGbl-1, D2jl_z(1:Nz,rowLoc), ADD_VALUES, ierr )
    ENDDO 
  ENDDO 

  DEALLOCATE( colGbl )
  DEALLOCATE( colGbl_orig )

END SUBROUTINE 

