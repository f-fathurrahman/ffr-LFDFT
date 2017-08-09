MODULE m_Poisson_solve_DAGE

  IMPLICIT NONE 

  INTEGER :: N_t
  REAL(8) :: t_i, t_l, t_f
  REAL(8), ALLOCATABLE :: w_t(:), x_t(:)
  REAL(8), ALLOCATABLE :: F_xs(:,:,:), F_ys(:,:,:), F_zs(:,:,:)

END MODULE 


FUNCTION compute_F( t, x_bar, h ) RESULT( f )
!
  USE m_Faddeeva, ONLY : erfcx
  USE m_constants, ONLY : PI, EPS_SMALL
  IMPLICIT NONE 
  REAL(8) :: t, x_bar, h
  REAL(8) :: f
  !
  COMPLEX(8) :: z, w_iz

  f = 0.d0

  IF(x_bar < EPS_SMALL) THEN 
    f = sqrt(h) * erf(PI/(2.d0*h*t))
  ELSE 
    z = cmplx( PI/(2.d0*h*t), t*x_bar, kind=8 )
    w_iz = erfcx( z )
    f = exp( -t*t*x_bar*x_bar )
    f = f - REAL( exp( -t*t*x_bar*x_bar - z*z )*w_iz, kind=8) 
    f = sqrt(h)*f
  ENDIF 
END FUNCTION 


SUBROUTINE construct_F( axis, t_size, t_values, F_values )
  USE m_LF3d, ONLY : NN => LF3d_NN, &
                     grid_x => LF3d_grid_x, &
                     grid_y => LF3d_grid_y, &
                     grid_z => LF3d_grid_z, &
                     hh => LF3d_hh
  IMPLICIT NONE 
  INTEGER :: axis
  INTEGER :: t_size
  REAL(8) :: t_values(t_size)
  REAL(8) :: F_values(NN(axis),NN(axis),t_size)
  INTEGER :: i_t, i, j
  REAL(8), ALLOCATABLE :: grid(:)
  REAL(8) :: compute_F

  ! copy grid_* to grid
  ALLOCATE( grid(NN(axis)) )
  IF( axis == 1 ) THEN 
    grid(:) = grid_x(:)
  ELSEIF( axis == 2 ) THEN 
    grid(:) = grid_y(:)
  ELSEIF( axis == 3 ) THEN 
    grid(:) = grid_z(:)
  ELSE 
    WRITE(*,*) 'Invalid value to axis', axis
    STOP 
  ENDIF 

  DO i_t = 1,t_size
    DO i = 1,NN(axis)
      DO j = 1,NN(axis)
        F_values(i,j,i_t) = compute_F( t_values(i_t), abs( grid(i) - grid(j) ), hh(axis) )
      ENDDO 
    ENDDO 
  ENDDO 

  DEALLOCATE( grid )

END SUBROUTINE 


! Sundholm (JCP 132, 024102, 2010).
! Within two-divided regions:
!     ([t_i,t_l], [t_l,t_f]),
! num_points1, num_points2 quadrature points are made, respectively.
SUBROUTINE init_t_sampling( &
    num_points1, num_points2, t_i, t_l, t_f, t_values, w_t )

  USE m_constants, ONLY : PI
  IMPLICIT NONE 
  INTEGER :: num_points1, num_points2
  REAL(8) :: t_i, t_l, t_f
  REAL(8) :: t_values(num_points1 + num_points2)
  REAL(8) :: w_t(num_points1 + num_points2)
  !
  INTEGER :: j
  REAL(8) :: s_p, w_p
  REAL(8), ALLOCATABLE :: x_leg(:), w_leg(:)

  ! Linear coord region:  [t_i, t_l]
  ALLOCATE( x_leg(num_points1), w_leg(num_points1) )
  CALL init_gauss_legendre(t_i, t_l, num_points1, x_leg, w_leg)
  !
  DO j=1,num_points1
    t_values(j) = x_leg(j)
    w_t(j)      = w_leg(j)
  ENDDO 

  DEALLOCATE( x_leg, w_leg )

  ! Logarithmic coord region: [t_l, t_f]
  ALLOCATE( x_leg(num_points2), w_leg(num_points2) )
  CALL init_gauss_legendre( log(t_l), log(t_f), num_points2, x_leg, w_leg )

  ! Return the log-coord-partitioned points back to linear t-space.
  s_p = 0.d0
  w_p = 0.d0
  DO j = 1,num_points2
    s_p = x_leg(j)
    w_p = w_leg(j)
    x_leg(j) = exp(s_p)
    w_leg(j) = w_p * exp(s_p)
  ENDDO 

  DO j = 1,num_points2
    t_values(num_points1+j) = x_leg(j)
    w_t(num_points1+j)      = w_leg(j)
  ENDDO 

  DEALLOCATE( x_leg, w_leg )

END SUBROUTINE 



! Generate Gauss-Legendre N-points quadrature formula
SUBROUTINE init_gauss_legendre( x1, x2, N, x, w )
  USE m_constants, ONLY : PI
  IMPLICIT NONE 
  !
  REAL(8) :: x1, x2
  INTEGER :: N
  REAL(8) :: x(N), w(N)
  !
  REAL(8), PARAMETER :: EPS = 3.d-11
  !
  INTEGER :: m, j, i
  REAL(8) :: z1, z, xm, xl, pp, p3, p2, p1

  m = (N + 1)/2        ! The roots are symmetric in the interval, so
  xm = 0.5d0*(x2 + x1) ! we only have to find half of them.
  xl = 0.5d0*(x2 - x1)

  ! FIXME: Need this ?
  z1 = 100.d0  ! some initial number to get the while loop enter the first time
  pp = 0.d0    ! make pp visible outside for loop

  DO i = 1, m ! Loop over the desired roots.
    z = cos( pi*(i-0.25)/(n+0.5) )
    ! Starting with the above approximation to the ith root, we enter the main loop of
    ! refinement by Newton’s method.
    DO 
      p1 = 1.d0
      p2 = 0.d0
      DO j = 1, N     ! Loop up the recurrence relation to get the
        p3 = p2       ! Legendre polynomial evaluated at z.
        p2 = p1
        p1 = ( (2.d0*j - 1.d0)*z*p2 - (j - 1.d0)*p3 )/j
      ENDDO 
      ! p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
      ! by a standard relation involving also p2, the polynomial of one lower order.
      pp = N*(z*p1-p2)/(z*z-1.d0)
      z1 = z
      z  = z1 - p1/pp ! Newton’s method.
      IF( abs(z-z1) < EPS ) EXIT 
    ENDDO
    !
    x(i)     = xm - xl*z                  ! Scale the root to the desired interval,
    x(n-i+1) = xm + xl*z                  ! and put in its symmetric counterpart.
    w(i)     = 2.d0*xl/((1.d0-z*z)*pp*pp) ! Compute the weight
    w(n-i+1) = w(i)                       ! and its symmetric counterpart.
  ENDDO

END SUBROUTINE 


!
SUBROUTINE init_Poisson_solve_DAGE()
!
  USE m_LF3d, ONLY : NN => LF3d_NN
  USE m_Poisson_solve_DAGE
  IMPLICIT NONE 
  INTEGER :: num_points1, num_points2

  !
  ! Explicitly set and allocate variables here
  !
  t_i = 0.d0
  t_l = 1.d0
  t_f = 100000.d0
  !
  num_points1 = 50
  num_points2 = 50
  N_t = num_points1 + num_points2

  ALLOCATE( F_xs(NN(1),NN(1),N_t) )
  ALLOCATE( F_ys(NN(2),NN(2),N_t) )
  ALLOCATE( F_zs(NN(3),NN(3),N_t) )

  ALLOCATE( x_t(N_t), w_t(N_t) )
  CALL init_t_sampling( num_points1, num_points2, t_i, t_l, t_f, &
                               x_t, w_t )
  !
!  WRITE(*,*)
!  WRITE(*,*) 't_sampling:'
!  DO i = 1, 
!    WRITE(*,'(I4,1x,F18.10,1x,F18.10)') i, x_t(i), w_t(i)
!  ENDDO 

  CALL construct_F( 1, N_t, x_t, F_xs )
  CALL construct_F( 2, N_t, x_t, F_ys )
  CALL construct_F( 3, N_t, x_t, F_zs )

END SUBROUTINE 



!---------------------------------------------------
SUBROUTINE Poisson_solve_DAGE( density, potential )
  !
  USE m_constants, ONLY : PI
  USE m_Poisson_solve_DAGE, ONLY : N_t, w_t, F_xs, F_ys, F_zs, t_f
  USE m_LF3d, ONLY : NN => LF3d_NN, &
                     Npoints => LF3d_Npoints, &
                     lin2xyz => LF3d_lin2xyz, &
                     dVol => LF3d_dVol
  IMPLICIT NONE 
  REAL(8) :: density(NN(1),NN(2),NN(3))
  REAL(8) :: potential(Npoints)  ! FIXME: use Nx,Ny,Nz ??
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

  DO i_t = 1,N_t
    
    ! FIXME: use DGEMM 
    DO gg = 1,NN(3)
      T_g(:,:,gg) = matmul( F_xs(:,:,i_t), density(:,:,gg) )
      T_g2(:,:,gg) = matmul( T_g(:,:,gg), F_ys(:,:,i_t) )
    ENDDO 

    ! reorder
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
      potential(ip) = potential(ip) + w_t(i_t)*T_b2(a,b,g)*2.d0/sqrt(PI)
    ENDDO 

  ENDDO 

  DO ip = 1,Npoints
    a = lin2xyz(1,ip)
    b = lin2xyz(2,ip)
    g = lin2xyz(3,ip)
    potential(ip) = (PI/(t_f*t_f))*density(a,b,g) + potential(ip)
  ENDDO 

  potential(:) = potential(:)*sqrt(dVol)

  DEALLOCATE( T_g )
  DEALLOCATE( T_g2 )
  DEALLOCATE( T_b )
  DEALLOCATE( T_b2 )

END SUBROUTINE 

