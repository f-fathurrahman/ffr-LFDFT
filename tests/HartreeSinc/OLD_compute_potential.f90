SUBROUTINE compute_potential( t_size, t_values, w_t, F_xs, F_ys, F_zs, &
                              density, potential )
  !
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : NN => LF3d_NN, &
                     Npoints => LF3d_Npoints, &
                     lin2xyz => LF3d_lin2xyz, &
                     xyz2lin => LF3d_xyz2lin
  IMPLICIT NONE 
  INTEGER :: t_size
  REAL(8) :: t_values(t_size)
  REAL(8) :: w_t(t_size)
  REAL(8) :: F_xs(t_size,NN(1),NN(1))
  REAL(8) :: F_ys(t_size,NN(2),NN(2))
  REAL(8) :: F_zs(t_size,NN(3),NN(3))
  REAL(8) :: density(Npoints)
  REAL(8) :: potential(Npoints)
  !
  INTEGER :: i_t, a, b, g, aa, bb, gg, ip, ipp
  REAL(8) :: s_F_xs, s_F_ys, s_F_zs, ss

  DO ip = 1, Npoints
    a = lin2xyz(1,ip)
    b = lin2xyz(2,ip)
    g = lin2xyz(3,ip)
    ss = 0.d0
    !
    DO i_t = 1, t_size
      !
      s_F_xs = 0.d0
      DO aa = 1,NN(1)
        ipp = xyz2lin(aa,b,g)
        s_F_xs = s_F_xs + F_xs(i_t,a,aa)*density(ipp)
      ENDDO 
      !
      s_F_ys = 0.d0
      DO bb = 1,NN(2)
        ipp = xyz2lin(a,bb,g)
        s_F_ys = s_F_ys + F_ys(i_t,b,bb)*density(ipp)
      ENDDO 
      !
      s_F_zs = 0.d0
      DO gg = 1,NN(3)
        ipp = xyz2lin(a,b,gg)
        s_F_zs = s_F_zs + F_zs(i_t,g,gg)*density(ipp)
      ENDDO
      !
      ss = ss + 2.d0*w_t(i_t)/sqrt(PI) * s_F_xs * s_F_ys * s_F_zs
    ENDDO 
    potential(ip) = ss
  ENDDO 

END SUBROUTINE 

