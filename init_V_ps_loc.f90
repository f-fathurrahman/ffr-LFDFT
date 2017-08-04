SUBROUTINE init_V_ps_loc()
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     lingrid => LF3d_lingrid
  USE m_hamiltonian, ONLY : V_ps_loc
  USE m_PsPot, ONLY : Ps_HGH_Params
  USE m_atoms, ONLY : atpos => AtomicCoords, atm2species, Natoms
  USE m_Ps_HGH, ONLY : hgh_eval_Vloc_R
  IMPLICIT NONE 
  INTEGER :: ip, ia, isp
  REAL(8) :: dr

  V_ps_loc(:) = 0.d0

  DO ia = 1,Natoms
    isp = atm2species(ia)
    DO ip = 1, Npoints
      CALL calc_dr_1pnt( atpos(:,ia), lingrid(:,ip), dr )
      V_ps_loc(ip) = V_ps_loc(ip) + hgh_eval_Vloc_R( Ps_HGH_Params(isp), dr )
    ENDDO 
  ENDDO 

  WRITE(*,*) 'sum(V_ps_loc) = ', sum(V_ps_loc)

END SUBROUTINE 

