
! Should be called only by master process
SUBROUTINE add_nabla2_y( Nx, Ny, Nz, D2jl_y, nabla2 )

  IMPLICIT NONE  

#include <petsc/finclude/petsc.h>
  
  INTEGER :: Nx, Ny, Nz
  REAL(8) :: D2jl_y(Ny,Ny)
  Mat :: nabla2
  !
  INTEGER :: irow, iy, ii, ix, iz
  PetscInt, ALLOCATABLE :: colGbl(:), colGbl_orig(:)
  PetscInt :: rowGbl, rowLoc, rowGbl_start, rowGbl_stop
  PetscInt :: Npoints

  PetscErrorCode :: ierr

  WRITE(*,*) 'Adding nabla2_y'

  ALLOCATE( colGbl(Ny) )
  ALLOCATE( colGbl_orig(Ny) )
  
    
  ! Initialize pattern for column indices
  colGbl_orig(1) = 1
  DO ii = 2,Ny
    colGbl_orig(ii) = colGbl_orig(ii-1) + Nz
  ENDDO
  colGbl(:) = colGbl_orig(:)
  !WRITE(*,*) 'colGbl_orig = ', colGbl_orig

  rowGbl = 0

!  rowLoc = 1
!  DO iy = 1, Ny
!    rowGbl = rowGbl + 1
!    colGbl(:) = colGbl_orig(:) + iy - 1
!    WRITE(*,*) rowGbl, ':', colGbl
!    CALL MatSetValues( nabla2, 1, rowGbl-1, Ny, colGbl-1, D2jl_y(1:Ny,rowLoc), ADD_VALUES, ierr )
!  ENDDO

!  rowLoc = 2
!  DO iy = 1, Ny
!    rowGbl = rowGbl + 1
!    colGbl(:) = colGbl_orig(:) + iy - 1
!    WRITE(*,*) rowGbl, ':', colGbl
!    CALL MatSetValues( nabla2, 1, rowGbl-1, Ny, colGbl-1, D2jl_y(1:Ny,rowLoc), ADD_VALUES, ierr )
!  ENDDO

  !
!  ix = 2
!  rowLoc = 1
!  DO iy = 1, Ny
!    rowGbl = rowGbl + 1
!    colGbl(:) = colGbl_orig(:) + iy - 1 + (ix-1)*Ny*Nz
!    WRITE(*,*) rowGbl, ':', colGbl
!    CALL MatSetValues( nabla2, 1, rowGbl-1, Ny, colGbl-1, D2jl_y(1:Ny,rowLoc), ADD_VALUES, ierr )
!  ENDDO

  DO ix = 1,Nx
    !WRITE(*,*)
    DO rowLoc = 1,Ny
      DO iz = 1, Nz
        rowGbl = rowGbl + 1
        colGbl(:) = colGbl_orig(:) + iz - 1 + (ix-1)*Ny*Nz
        !WRITE(*,*) rowGbl, ':', colGbl
        CALL MatSetValues( nabla2, 1, rowGbl-1, Ny, colGbl-1, D2jl_y(1:Ny,rowLoc), ADD_VALUES, ierr )
      ENDDO
    ENDDO
  ENDDO

  DEALLOCATE( colGbl )
  DEALLOCATE( colGbl_orig )

END SUBROUTINE 

