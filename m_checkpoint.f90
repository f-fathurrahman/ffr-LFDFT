MODULE m_checkpoint

  IMPLICIT NONE 

  INTEGER :: CHK_Npoints
  INTEGER :: CHK_NN(3)
  INTEGER :: CHK_Nstates

  INTEGER :: CHK_LF3d_TYPE

  REAL(8) :: CHK_LL(3)
  REAL(8) :: CHK_AA(3), CHK_BB(3)
  REAL(8) :: CHK_hh(3)

  REAL(8) :: CHK_dVol
  REAL(8) :: CHK_GRID_SHIFT(3)

  NAMELIST /CHK_DATA/ CHK_Npoints, CHK_NN, CHK_Nstates, CHK_LF3d_TYPE, &
                      CHK_LL, CHK_AA, CHK_BB, CHK_hh, CHK_dVol, CHK_GRID_SHIFT

  INTEGER, PARAMETER :: IU_CHK = 100
  CHARACTER(56) :: CHK_FILENAME = 'CHK.dat'

  INTEGER, PARAMETER :: IU_GRID = 99
  CHARACTER(56) :: CHK_GRID_FILENAME = 'CHK_GRID.dat'

END MODULE 


SUBROUTINE write_checkpoint()
  USE m_checkpoint
  USE m_LF3d
  USE m_states
  IMPLICIT NONE 

  CHK_Npoints = LF3d_Npoints
  CHK_NN(:) = LF3d_NN(:)
  CHK_Nstates = Nstates

  CHK_LF3d_TYPE = LF3d_TYPE

  CHK_LL(:) = LF3d_LL(:)
  CHK_AA(:) = LF3d_AA(:)
  CHK_BB(:) = LF3d_BB(:)
  CHK_hh(:) = LF3d_hh(:)

  CHK_dVol = LF3d_dVol
  CHK_GRID_SHIFT(:) = LF3d_GRID_SHIFT(:)

  OPEN(unit=IU_CHK, file=CHK_FILENAME, form='formatted', action='write')
  WRITE(IU_CHK,nml=CHK_DATA) 
  CLOSE(IU_CHK)

  CALL write_grid_1d()

END SUBROUTINE 

SUBROUTINE write_grid_1d()
!  USE m_LF3d, ONLY : LF3d_grid_x, LF3d_grid_y, LF3d_grid_z
  IMPLICIT NONE 
!  NAMELIST /CHK_GRID_POINTS_1d/ LF3d_grid_x, LF3d_grid_y, LF3d_grid_z
  
!  OPEN(unit=IU_CHK, file=CHK_GRID_FILENAME, form='formatted', action='write')
!  WRITE(IU_CHK,nml=CHK_GRID_POINTS_1d) 
!  CLOSE(IU_CHK)
END SUBROUTINE 

SUBROUTINE read_checkpoint()
  USE m_checkpoint
  IMPLICIT NONE 
  !
  OPEN(unit=IU_CHK, file=CHK_FILENAME, form='formatted', status='old', action='read')
  READ(IU_CHK,nml=CHK_DATA)
  CLOSE(IU_CHK)
END SUBROUTINE 


