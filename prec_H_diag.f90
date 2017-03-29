SUBROUTINE prec_H_diag( v, pv )

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_hamiltonian, ONLY : K_diag, V_ps_loc, V_Hartree, V_xc

  IMPLICIT NONE
  ! Arguments
  REAL(8) :: v(Npoints)
  REAL(8) :: pv(Npoints)
  ! Local
  INTEGER :: ip

  DO ip = 1, Npoints
    ! use absolute value ?
    pv(ip) = v(ip) / abs( K_diag(ip) + V_ps_loc(ip) + V_Hartree(ip) + V_xc(ip) )
  ENDDO 

END SUBROUTINE 
