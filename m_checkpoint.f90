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

  NAMELIST /CHK_DATA/ CHK_Npoints, CHK_NN, CHK_Nstates, CHK_LF3d_TYPE, &
                      CHK_LL, CHK_AA, CHK_BB, CHK_hh, CHK_dVol

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

  WRITE(111,nml=CHK_DATA) 

END SUBROUTINE 



