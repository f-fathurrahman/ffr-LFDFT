SUBROUTINE init_betaNL()
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     lingrid => LF3d_lingrid, &
                     LL => LF3d_LL
  USE m_PsPot, ONLY : betaNL, &
                      Ps => Ps_HGH_Params
  USE m_atoms, ONLY : atpos => AtomicCoords
  USE m_Ps_HGH, ONLY : hgh_eval_proj_R
  IMPLICIT NONE 
  INTEGER :: ia, isp, l, m, iprj
  INTEGER :: Np_beta, ip, ibeta
  REAL(8) :: dr_vec(3)
  REAL(8) :: dr
  REAL(8) :: Ylm_real
  
  !!! SPECIAL CASE !!!!
  ia = 1
  isp = 1
  l = 0
  m = 0
  iprj = 1
  betaNL(:,:) = 0.d0
  Np_beta = 0
  ibeta = 1

  DO ip = 1,Npoints
    CALL calc_dr_periodic_1pnt( LL, atpos(:,ia), lingrid(:,ip), dr_vec )
    dr = sqrt( dr_vec(1)**2 + dr_vec(2)**2 + dr_vec(3)**2 )
    IF( dr <= Ps(isp)%rcut_NL(l) ) THEN 
      Np_beta = Np_beta + 1
      betaNL(ip,ibeta) = hgh_eval_proj_R( Ps(isp), l, iprj, dr ) * Ylm_real( l, m, dr_vec )
    ENDIF 
  ENDDO 
  WRITE(*,*) 'Np_beta = ', Np_beta
  WRITE(*,*) 'sum(betaNL) = ', sum(betaNL)
END SUBROUTINE 

