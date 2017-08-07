PROGRAM HartreeSinc_main

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_constants, ONLY : PI
  IMPLICIT NONE 
  ! LF sinc parameters
  REAL(8) :: scaling(3)
  INTEGER :: NN(3)
  ! t sampling stuffs
  INTEGER :: num_points1, num_points2
  REAL(8) :: t_i, t_l, t_f
  INTEGER :: t_size
  REAL(8), ALLOCATABLE :: t_values(:), w_t(:)
  ! gaussian centers
  INTEGER :: num_gaussian
  REAL(8), ALLOCATABLE :: positions(:,:), exponents(:), coefs(:)
  !
  REAL(8), ALLOCATABLE :: density(:), anal_pot(:), potential(:)
  REAL(8) :: density_norm, anal_energy
  REAL(8), ALLOCATABLE :: F_xs(:,:,:), F_ys(:,:,:), F_zs(:,:,:)
  !
  INTEGER :: i, ip
  CHARACTER(56) :: chars_arg
  INTEGER :: N_in
  INTEGER :: iargc

  !----------------------- Input spec -------------------------------!
  IF( iargc() /= 1 ) THEN 
    WRITE(*,*) 'Exactly one argument must be given:', iargc()
    STOP 
  ENDIF 
  CALL getarg(1, chars_arg )
  READ( chars_arg, *) N_in
  WRITE(*,*) 'N_in = ', N_in

  NN(:) = (/ N_in, N_in, N_in /)
  scaling(:) = (/1.d0, 1.d0, 1.d0/)*(8.d0/(NN(1)-1))
  
  num_points1 = 10
  num_points2 = 100
  t_i = 0.d0
  t_l = 1.d0
  t_f = 10000.d0

  !
  num_gaussian = 1
  ALLOCATE( positions(3,num_gaussian) )
  ALLOCATE( coefs(num_gaussian) )
  ALLOCATE( exponents(num_gaussian) )

  ! First Gaussian
  positions(:,1) = (/ 0.d0, 0.d0, 0.d0 /)
  coefs(1)       = 1.d0
  exponents(1)   = sqrt(2.d0)

  ! Second Gaussian
!  positions(:,2) = (/ -1.d0, 0.d0, 0.d0 /)
!  coefs(2)       = 1.d0
!  exponents(2)   = 1.5d0

  !---------------------- end of input specs ------------------------!


  CALL init_LF3d_sinc( NN, scaling )
  CALL info_LF3d()

  ! t_sampling
  t_size = num_points1 + num_points2
  ALLOCATE( t_values(t_size), w_t(t_size) )
  CALL HartreeSinc_t_sampling( num_points1, num_points2, t_i, t_l, t_f, &
                               t_values, w_t )
  !WRITE(*,*)
  !WRITE(*,*) 't_sampling:'
  !DO i = 1, t_size
  !  WRITE(*,'(I4,1x,F18.10,1x,F18.10)') i, t_values(i), w_t(i)
  !ENDDO 

  !
  ALLOCATE( density(Npoints), anal_pot(Npoints) )
  CALL init_density( num_gaussian, positions, coefs, exponents, density, anal_pot )

  density_norm = 0.d0
  anal_energy = 0.d0
  DO ip = 1, Npoints
    density_norm = density_norm + density(ip)
    anal_energy  = anal_energy + 0.5d0*anal_pot(ip)*density(ip)
  ENDDO 
  density_norm = density_norm*sqrt(dVol)
  anal_energy  = anal_energy*sqrt(dVol)
  WRITE(*,'(1x,A,F18.10)') 'norm of density:', density_norm

  ALLOCATE( F_xs(NN(1),NN(1),t_size ) )
  ALLOCATE( F_ys(NN(2),NN(2),t_size ) )
  ALLOCATE( F_zs(NN(3),NN(3),t_size ) )

  CALL construct_F( 1, t_size, t_values, F_xs )
  CALL construct_F( 2, t_size, t_values, F_ys )
  CALL construct_F( 3, t_size, t_values, F_zs )

!  WRITE(*,*)
!  WRITE(*,*) 'In main:'
!  WRITE(*,*) 'sum(F_xs) = ', sum(F_xs)
!  WRITE(*,*) 'sum(F_ys) = ', sum(F_ys)
!  WRITE(*,*) 'sum(F_zs) = ', sum(F_zs)

!  WRITE(*,*) 'maxval(F_xs) = ', maxval(F_xs)
!  WRITE(*,*) 'minval(F_xs) = ', minval(F_xs)

  ALLOCATE( potential(Npoints) )
  CALL compute_potential( t_size, w_t, F_xs, F_ys, F_zs, &
                          density, potential )
!  DO ip = 1,Npoints
!    potential(ip) = (PI/(t_f*t_f))*density(ip) + potential(ip)
!  ENDDO 
  WRITE(*,*) 'sum(anal_pot)  = ', sum(anal_pot)
  WRITE(*,*) 'sum(potential) = ', sum(potential)
  WRITE(*,'(1x,A,F18.10)') 'numeric        :', 0.5d0*sum( density(:)*potential(:) ) * sqrt(dVol)
  WRITE(*,'(1x,A,F18.10)') 'Analytic energy:', anal_energy
  WRITE(*,'(1x,A,F18.10)') 'analytic: ', exponents(1)/sqrt(2.d0*PI)

  DEALLOCATE( potential )
  DEALLOCATE( F_xs, F_ys, F_zs )
  DEALLOCATE( t_values, w_t )
  DEALLOCATE( density, anal_pot )


END PROGRAM 

