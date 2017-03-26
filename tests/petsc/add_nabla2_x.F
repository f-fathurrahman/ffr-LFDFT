
! Should be called only by master process
SUBROUTINE add_nabla2_x( Nx, Ny, Nz, D2jl_x, nabla2 )

  IMPLICIT NONE  

#include <petsc/finclude/petsc.h>
  
  INTEGER :: Nx, Ny, Nz
  REAL(8) :: D2jl_x(Nx,Nx)
  Mat :: nabla2
  !
  INTEGER :: ix, irow
  INTEGER, ALLOCATABLE :: colGbl(:), colGbl_orig(:)
  INTEGER :: rowGbl

  PetscErrorCode :: ierr

  WRITE(*,*) 'Adding nabla2_x'

  ALLOCATE( colGbl(Nx) )
  ALLOCATE( colGbl_orig(Nx) )
  
  ! Initialize pattern for column indices
  colGbl_orig(1) = 1
  DO ix = 2,Nx
    colGbl_orig(ix) = colGbl_orig(1) + Ny*Nz*(ix-1)
  ENDDO

  colGbl(:) = colGbl_orig(:)
  !
  DO ix = 1,Nx
    DO irow = 1,Ny*Nz
      colGbl(:) = colGbl_orig(:) + irow - 1
      rowGbl = irow + (ix-1)*Ny*Nz
      ! Because D2jl_x is symmetric, we access them by column instead of by row
      CALL MatSetValues( nabla2, 1, rowGbl-1, Nx, colGbl-1, D2jl_x(1:Nx,ix), ADD_VALUES, ierr )
    ENDDO
  ENDDO

  DEALLOCATE( colGbl )
  DEALLOCATE( colGbl_orig )

END SUBROUTINE 

