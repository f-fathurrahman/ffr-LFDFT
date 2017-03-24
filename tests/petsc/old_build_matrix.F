! Very slow :-(
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


  ! Also slow :-(
  IF ( rank == 0 ) THEN 
  DO colLoc = 1,Nx
    colGbl_start = (colLoc-1)*(Ny*Nz)+1
    colGbl_stop = colLoc*Ny*Nz
    DO rowLoc = 1,Nx
      DO colGbl = colGbl_start, colGbl_stop
        rowGbl = (rowLoc-1)*Ny*Nz + colGbl - (colLoc-1)*Ny*Nz
        CALL MatSetValue( nabla2_x, rowGbl-1, colGbl-1, D2jl_x(rowLoc,colLoc), INSERT_VALUES, ierr )
      ENDDO
    ENDDO
  ENDDO ! colLoc
