SUBROUTINE info_LF3d()
  USE m_LF3d
  IMPLICIT NONE
  INTEGER :: N, Nx, Ny, Nz

  WRITE(*,*) 'Grid related variables'
  WRITE(*,*) '^^^^^^^^^^^^^^^^^^^^^^'
  WRITE(*,*)
  WRITE(*,'(1x,A,3F10.5)') 'Box size            = ', LF3d_LL(:)
  WRITE(*,'(1x,A,3F10.5)') 'Grid spacing        = ', LF3d_hh
  WRITE(*,'(1x,A,3I5)')    'Sampling Nx, Ny, Nz = ', LF3d_NN(:)
  WRITE(*,'(1x,A,I10)')    'Number of points    = ', LF3d_Npoints
  WRITE(*,'(1x,A,F10.5)')  'dVol                = ', LF3d_dVol

  ! x-direction --------------------------------------
  Nx = LF3d_NN(1)
  Ny = LF3d_NN(2)
  Nz = LF3d_NN(3)
  WRITE(*,'(/,1x,A,/)') 'Some grid points in x, y, and z directions:'
  WRITE(*,fmt=999) 1, LF3d_grid_x(1), 1, LF3d_grid_y(1), 1, LF3d_grid_z(1)
  WRITE(*,fmt=999) 2, LF3d_grid_x(2), 2, LF3d_grid_y(2), 2, LF3d_grid_z(2)
  WRITE(*,*) '      ..   .......      ..   .......      ..   .......'
  WRITE(*,fmt=999) Nx-1, LF3d_grid_x(Nx-1), Ny-1, LF3d_grid_y(Ny-1), Nz-1, LF3d_grid_z(Nz-1)
  WRITE(*,fmt=999) Nx,   LF3d_grid_x(Nx), Ny, LF3d_grid_y(Ny), Nz, LF3d_grid_z(Nz)

  999 FORMAT(1x,I8,F10.5,I8,F10.5,I8,F10.5)

  ! This should be called only for periodic case
  IF( LF3d_TYPE == LF3d_PERIODIC ) THEN
    N = LF3d_Npoints
    WRITE(*,'(/,1x,A,/)') 'Some G2 values:'
    WRITE(*,'(1x,I10,F10.5)') 1, LF3d_G2(1)
    WRITE(*,'(1x,I10,F10.5)') 2, LF3d_G2(2)
    WRITE(*,*) '       ...   .......'
    WRITE(*,'(1x,I10,F10.5)') N-1, LF3d_G2(N-1)
    WRITE(*,'(1x,I10,F10.5)') N,   LF3d_G2(N)
  ENDIF 

END SUBROUTINE

