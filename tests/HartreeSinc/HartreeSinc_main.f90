PROGRAM HartreeSinc_main

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
  INTEGER :: i

  !----------------------- Input spec -------------------------------!
  scaling(:) = (/0.2d0, 0.2d0, 0.2d0/)
  NN(:) = (/ 25, 25, 25 /)
  
  num_points1 = 5
  num_points2 = 10
  t_i = 0.0
  t_l = 2.0
  t_f = 10000.0

  !
  num_gaussian = 2
  ALLOCATE( positions(3,num_gaussian) )
  ALLOCATE( coefs(num_gaussian) )
  ALLOCATE( exponents(num_gaussian) )

  ! First Gaussian
  positions(:,1) = (/ 0.d0, 0.d0, 0.d0 /)
  coefs(1)       = 1.d0
  exponents(1)   = 3.d0

  ! Second Gaussian
  positions(:,1) = (/ 1.5d0, 0.d0, 0.d0 /)
  coefs(2)       = 1.d0
  exponents(2)   = 4.d0

  !---------------------- end of input specs ------------------------!


  CALL init_LF3d_sinc( NN, scaling )
  CALL info_LF3d()

  ! t_sampling
  t_size = num_points1 + num_points2
  ALLOCATE( t_values(t_size), w_t(t_size) )
  CALL HartreeSinc_t_sampling( num_points1, num_points2, t_i, t_l, t_f, &
                               t_values, w_t )
  WRITE(*,*)
  WRITE(*,*) 't_sampling:'
  DO i = 1, t_size
    WRITE(*,'(I4,1x,F18.10,1x,F18.10)') i, t_values(i), w_t(i)
  ENDDO 


END PROGRAM 

