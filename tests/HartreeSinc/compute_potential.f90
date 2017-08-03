SUBROUTINE compute_potential( t_size, w_t, F_xs, F_ys, F_zs, &
                              density, potential )
  !
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : NN => LF3d_NN, &
                     Npoints => LF3d_Npoints, &
                     lin2xyz => LF3d_lin2xyz, &
                     xyz2lin => LF3d_xyz2lin
  IMPLICIT NONE 
  INTEGER :: t_size
  REAL(8) :: w_t(t_size)
  REAL(8) :: F_xs(t_size,NN(1),NN(1))
  REAL(8) :: F_ys(t_size,NN(2),NN(2))
  REAL(8) :: F_zs(t_size,NN(3),NN(3))
  REAL(8) :: density(Npoints)
  REAL(8) :: potential(Npoints)
  !
  INTEGER :: i_t, a, b, g, aa, bb, gg, ip, ipp
  REAL(8), ALLOCATABLE :: Fd(:,:)
  REAL(8), ALLOCATABLE :: FFd(:)
  REAL(8) :: FFFd, wFFFd

!  WRITE(*,*)
!  WRITE(*,*) 'In compute_potential:'
!  WRITE(*,*) 'sum(F_xs) = ', sum(F_xs)
!  WRITE(*,*) 'sum(F_ys) = ', sum(F_ys)
!  WRITE(*,*) 'sum(F_zs) = ', sum(F_zs)
!  WRITE(*,*) 'maxval(F_xs) = ', maxval(F_xs)
!  WRITE(*,*) 'minval(F_xs) = ', minval(F_xs)

!  WRITE(*,*) 'NN(:) = ', NN(:)

!  STOP 

  ALLOCATE( Fd(NN(1),NN(2)) )
  ALLOCATE( FFd(NN(1)) )

  DO ip = 1, Npoints
    a = lin2xyz(1,ip)
    b = lin2xyz(2,ip)
    g = lin2xyz(3,ip)
    wFFFd = 0.d0
    !
    DO i_t = 1, t_size
      !
      Fd(:,:) = 0.d0
      DO gg = 1,NN(3)
        DO aa = 1,NN(1)
        DO bb = 1,NN(2)
          ipp = xyz2lin(aa,bb,gg)
          Fd(aa,bb) = Fd(aa,bb) + F_zs(i_t,g,gg)*density(ipp)
        ENDDO 
        ENDDO 
      ENDDO 
      !
      FFd(:) = 0.d0
      DO bb = 1,NN(2)
        DO aa = 1,NN(1)
          FFd(aa) = FFd(aa) + F_ys(i_t,b,bb)*Fd(aa,bb)
        ENDDO
      ENDDO 
      !
      FFFd = 0.d0
      DO aa = 1,NN(1)
        FFFd = FFFd + F_xs(i_t,a,aa)*FFd(aa)
        !WRITE(*,'(1x,I4,2ES18.10)') aa, FFFd, F_xs(i_t,a,aa)
      ENDDO 
      wFFFd = wFFFd + w_t(i_t)*FFFd
      !WRITE(*,'(I8,5ES18.10)') i_t, w_t(i_t), sum(Fd), sum(FFd), FFFd, wFFFd
    ENDDO 
    potential(ip) = 2.d0*wFFFd/sqrt(PI)
    !WRITE(*,'(1x,A,I8)') 'done ip = ', ip
  ENDDO 

  DEALLOCATE( Fd )
  DEALLOCATE( FFd )

END SUBROUTINE 

