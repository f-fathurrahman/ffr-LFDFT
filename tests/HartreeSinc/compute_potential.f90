SUBROUTINE compute_potential( t_size, w_t, F_xs, F_ys, F_zs, &
                              density, potential )
  !
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : NN => LF3d_NN, &
                     Npoints => LF3d_Npoints, &
                     lin2xyz => LF3d_lin2xyz
  IMPLICIT NONE 
  INTEGER :: t_size
  REAL(8) :: w_t(t_size)
  REAL(8) :: F_xs(NN(1),NN(1),t_size)
  REAL(8) :: F_ys(NN(2),NN(2),t_size)
  REAL(8) :: F_zs(NN(3),NN(3),t_size)
  REAL(8) :: density(NN(1),NN(2),NN(3))
  REAL(8) :: potential(Npoints)
  !
  INTEGER :: i_t, a, b, g, bb, gg, ip
  REAL(8), ALLOCATABLE :: T_g(:,:,:), T_g2(:,:,:), T_b(:,:,:), T_b2(:,:,:)
  
  ALLOCATE( T_g(NN(1),NN(2),NN(3)) )
  ALLOCATE( T_g2(NN(1),NN(2),NN(3)) )
  ALLOCATE( T_b(NN(1),NN(3),NN(2)) )
  ALLOCATE( T_b2(NN(1),NN(3),NN(2)) )

  T_g(:,:,:)  = 0.d0
  T_g2(:,:,:) = 0.d0
  T_b(:,:,:)  = 0.d0
  T_b2(:,:,:) = 0.d0

  potential(:) = 0.d0

  WRITE(*,*) 'Using matmul:'

  DO i_t = 1,t_size

    DO gg = 1,NN(3)
      T_g(:,:,gg) = matmul( F_xs(:,:,i_t), density(:,:,gg) )
      T_g2(:,:,gg) = matmul( T_g(:,:,gg), F_ys(:,:,i_t) )
    ENDDO 

    ! reorder ?
    DO bb = 1,NN(2)
    DO gg = 1,NN(3)
      T_b(:,gg,bb) = T_g2(:,bb,gg)
    ENDDO 
    ENDDO 

    DO bb = 1,NN(2)
      T_b2(:,:,bb) = matmul( T_b(:,:,bb), F_zs(:,:,i_t) )
    ENDDO 


    DO ip = 1,Npoints
      a = lin2xyz(1,ip)
      b = lin2xyz(2,ip)
      g = lin2xyz(3,ip)
      potential(ip) = potential(ip) + w_t(i_t)*T_b2(a,b,g)
    ENDDO 
    ! no need to multiply by 2.0/sqrt(pi), it is included already in w_t

  ENDDO 

END SUBROUTINE 


SUBROUTINE transpose_yz( NN, T )
  IMPLICIT NONE 
  INTEGER :: NN(3)
  REAL(8) :: T(NN(1),NN(2),NN(3))
  INTEGER :: aa, bb, gg
  REAL(8) :: t1, t2

  DO gg = 1,NN(3)
  DO bb = 1,NN(2)
  DO aa = 1,NN(1)
    t1 = T(aa,bb,gg)
    t2 = T(aa,gg,bb)
    T(aa,bb,gg) = t1
    T(aa,gg,bb) = t2
  ENDDO 
  ENDDO 
  ENDDO 

END SUBROUTINE 
