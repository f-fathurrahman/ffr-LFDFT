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
  REAL(8) :: s_F

  DO ip = 1, Npoints
    a = lin2xyz(1,ip)
    b = lin2xyz(2,ip)
    g = lin2xyz(3,ip)
    DO i_t = 1, t_size
      !
      s_F = 0.d0
      DO gg = 1,NN(3)
      DO bb = 1,NN(2)
      DO aa = 1,NN(1)
        ipp = xyz2lin(aa,bb,gg)
        s_F = s_F + density(ipp) * F_xs(i_t,a,aa) * F_ys(i_t,b,bb) * F_zs(i_t,g,gg)
      ENDDO 
      ENDDO 
      ENDDO 
      s_F = s_F*2.d0/sqrt(PI)*w_t(t_size)
      !
    ENDDO 
  ENDDO 

END SUBROUTINE 

