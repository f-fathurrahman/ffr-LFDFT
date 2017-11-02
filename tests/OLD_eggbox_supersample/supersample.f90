SUBROUTINE supersample( fin, fout )
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     lingrid => LF3d_lingrid, &
                     LL => LF3d_LL, &
                     hh => LF3d_hh
  USE m_LF3d_supersample, ONLY : lingrid_ss => LF3d_lingrid_ss, &
                                 Npoints_ss => LF3d_Npoints_ss, &
                                 Nsupersample
  USE m_Ps_HGH, ONLY : hgh_eval_Vloc_R_short
  USE m_PsPot, ONLY : Ps => Ps_HGH_Params
  USE m_atoms, ONLY : atpos => AtomicCoords
  USE m_grid_atom, ONLY : Ngrid_atom, idxa => idx_grid_atom
  USE m_grid_ss_atom, ONLY : Ngrid_ss_atom, idx_ss => idx_grid_ss_atom
  IMPLICIT NONE 
  REAL(8) :: fin(Npoints)
  REAL(8) :: fout(Npoints)
  !
  INTEGER :: ip, ip_ss, ix, iy, iz
  INTEGER :: isp, ip_ss_a, ip_a
  !
  REAL(8) :: x, y, z, x_ss, y_ss, z_ss
  REAL(8) :: center(3), dr_vec(3), dr
  REAL(8) :: ff
  REAL(8) :: sinc

  isp = 1 

  fout(:) = 0.d0
  DO ip_a = 1,Ngrid_atom
    ff = 0.d0
    ip = idxa(ip_a)
    x = lingrid(1,ip)
    y = lingrid(2,ip)
    z = lingrid(3,ip)
    DO ip_ss_a = 1,Ngrid_ss_atom
      ip_ss = idx_ss(ip_ss_a)
      x_ss = lingrid_ss(1,ip_ss)
      y_ss = lingrid_ss(2,ip_ss)
      z_ss = lingrid_ss(3,ip_ss)
      !
      CALL calc_dr_periodic_1pnt( LL, atpos(:,1), (/x_ss,y_ss,z_ss/), dr_vec )
      dr = sqrt( dr_vec(1)**2 + dr_vec(2)**2 + dr_vec(3)**2 )
      !
      ff = ff + sinc((x-x_ss)/hh(1)) * sinc((y-y_ss)/hh(2)) * sinc((z-z_ss)/hh(3)) * &
                hgh_eval_Vloc_R_short(Ps(isp),dr) / Nsupersample**3 
    ENDDO 
    fout(ip) = ff
  ENDDO 

END SUBROUTINE 


