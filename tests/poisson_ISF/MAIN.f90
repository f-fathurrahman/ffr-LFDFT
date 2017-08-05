PROGRAM MAIN

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     lin2xyz => LF3d_lin2xyz, &
                     dVol => LF3d_dVol
  USE Poisson_Solver
  IMPLICIT NONE 
  ! LF sinc parameters
  REAL(8) :: hh(3)
  INTEGER :: NN(3)
  ! gaussian centers
  INTEGER :: num_gaussian
  REAL(8), ALLOCATABLE :: positions(:,:), exponents(:), coefs(:)
  !
  REAL(8), ALLOCATABLE :: density(:), anal_pot(:), potential(:)
  REAL(8) :: density_norm, anal_energy
  !
  INTEGER :: nfft1, nfft2, nfft3
  INTEGER :: n1k, n2k, n3k
  INTEGER :: itype_scf
  REAL(8) :: hgrid
  REAL(8) :: ehartree
  !
  REAL(8), ALLOCATABLE :: rhopot(:,:,:)
  REAL(8), POINTER :: kernel(:)
  REAL(8) :: fake_arr(1)
  REAL(8) :: fake_exc
  REAL(8) :: fake_vxc
  !
  INTEGER :: i, j, k, ip
  INTEGER :: N_in
  CHARACTER(56) :: chars_N
  INTEGER :: iargc

  !----------------------- Input spec -------------------------------!

  IF( iargc() /= 1 ) THEN 
    WRITE(*,*) 'ERROR: exactly one argument must be given:'
    WRITE(*,*) 'N'
    STOP 
  ENDIF 

  CALL getarg(1,chars_N)
  READ(chars_N,*) N_in

  NN(:) = N_in
  hh(:) = (/1.d0, 1.d0, 1.d0/)*(16.d0/(NN(1)-1))
  
  !
  num_gaussian = 2
  ALLOCATE( positions(3,num_gaussian) )
  ALLOCATE( coefs(num_gaussian) )
  ALLOCATE( exponents(num_gaussian) )

  ! First Gaussian
  positions(:,1) = (/  1.d0, 0.d0, 0.d0 /)
  coefs(1)       = 1.d0
  exponents(1)   = sqrt(2.d0)

  ! Second Gaussian
  positions(:,2) = (/ -1.d0, 0.d0, 0.d0 /)
  coefs(2)       = 1.d0
  exponents(2)   = 10.d0

  !---------------------- end of input specs ------------------------!


  CALL init_LF3d_sinc( NN, hh )
  CALL info_LF3d()
 
  !
  ALLOCATE( density(Npoints), anal_pot(Npoints) )
  CALL init_density( num_gaussian, positions, coefs, exponents, density, anal_pot )

  density_norm = 0.d0
  anal_energy = 0.d0
  DO ip = 1, Npoints
    density_norm = density_norm + density(ip)
    anal_energy  = anal_energy + 0.5d0*anal_pot(ip)*density(ip)
  ENDDO 
  density_norm = density_norm*dVol
  anal_energy  = anal_energy*dVol

  CALL createKernel( 'F', NN(1),NN(2),NN(3), hh(1),hh(2),hh(3), 14, 0, 1, kernel)

  ALLOCATE( rhopot(NN(1),NN(2),NN(3)) )
  DO ip = 1,Npoints
    i = lin2xyz(1,ip)
    j = lin2xyz(2,ip)
    k = lin2xyz(3,ip)
    rhopot(i,j,k) = density(ip)
  ENDDO 
  
  call PSolver('F','G',0,1, NN(1),NN(2),NN(3),0, hh(1),hh(2),hh(3), &
                rhopot, kernel, fake_arr, ehartree, fake_exc,fake_vxc,0.d0,.false.,1)

  WRITE(*,*)
  WRITE(*,'(1x,A,F18.10)')  'norm of density:', density_norm
  WRITE(*,'(1x,A,F18.10)')  'analytic energy:', anal_energy
  WRITE(*,'(1x,A,F18.10)')  'ehartree = ', ehartree
  WRITE(*,'(1x,A,ES18.10)') 'diff     = ', abs(ehartree-anal_energy)

  ALLOCATE( potential(Npoints) )
  DO ip = 1,Npoints
    i = lin2xyz(1,ip)
    j = lin2xyz(2,ip)
    k = lin2xyz(3,ip)
    potential(ip) = rhopot(i,j,k)
  ENDDO 
  WRITE(*,'(1x,A,F18.10)') 'Calc = ', 0.5d0*sum( density(:)*potential(:) ) *dVol
  WRITE(*,'(1x,A,ES18.10)') 'diff = ', abs(0.5d0*sum( density(:)*potential(:) ) *dVol-anal_energy)

  DEALLOCATE( density, anal_pot )
  CALL dealloc_LF3d()

END PROGRAM 

