SUBROUTINE info_LF3d()
  USE m_LF3d
  IMPLICIT NONE
  INTEGER :: N

  WRITE(*,*) '--------------------------'
  WRITE(*,*) 'Basis function information'
  WRITE(*,*) '--------------------------'
  WRITE(*,*)
  WRITE(*,'(1x,A,3F18.10)') 'Box size = ', LF3d_LL(:)
  WRITE(*,'(1x,A,3I10)') 'Sampling Nx, Ny, Nz = ', LF3d_NN(:)
  WRITE(*,'(1x,A,I10)') 'Number of points = ', LF3d_Npoints
  WRITE(*,'(1x,A,3F18.10)') 'Grid spacing = ', LF3d_hh
  WRITE(*,'(1x,A,F18.10)') 'dVol = ', LF3d_dVol

  ! x-direction --------------------------------------
  N = LF3d_NN(1)
  WRITE(*,'(/,1x,A)') 'Some grid points in x-direction:'
  WRITE(*,*) 1, LF3d_grid_x(1)
  WRITE(*,*) 2, LF3d_grid_x(2)
  WRITE(*,*) '....'
  WRITE(*,*) N-1, LF3d_grid_x(N-1)
  WRITE(*,*) N,   LF3d_grid_x(N)

  ! y-direction --------------------------------------
  N = LF3d_NN(2)
  WRITE(*,'(/,1x,A)') 'Some grid points in y-direction:'
  WRITE(*,*) 1, LF3d_grid_y(1)
  WRITE(*,*) 2, LF3d_grid_y(2)
  WRITE(*,*) '....'
  WRITE(*,*) N-1, LF3d_grid_y(N-1)
  WRITE(*,*) N,   LF3d_grid_y(N)

  ! z-direction --------------------------------------
  N = LF3d_NN(3)
  WRITE(*,'(/,1x,A)') 'Some grid points in z-direction:'
  WRITE(*,*) 1, LF3d_grid_z(1)
  WRITE(*,*) 2, LF3d_grid_z(2)
  WRITE(*,*) '....'
  WRITE(*,*) N-1, LF3d_grid_z(N-1)
  WRITE(*,*) N,   LF3d_grid_z(N)

  ! FIXME This should be called only for periodic case
  ! G2 ------------------------------------------
  N = LF3d_Npoints
  WRITE(*,'(/,1x,A)') 'Some G2 values:'
  WRITE(*,*) 1, LF3d_G2(1)
  WRITE(*,*) 2, LF3d_G2(2)
  WRITE(*,*) '....'
  WRITE(*,*) N-1, LF3d_G2(N-1)
  WRITE(*,*) N,   LF3d_G2(N)

END SUBROUTINE

