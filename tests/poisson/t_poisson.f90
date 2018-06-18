! efefer 30 December 2015
PROGRAM t_poisson
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     lingrid => LF3d_lingrid, &
                     dVol => LF3d_dVol
  IMPLICIT NONE
  !
  INTEGER :: NN(3)
  REAL(8) :: LL(3)
  REAL(8), ALLOCATABLE :: rho(:) ! density
  REAL(8), ALLOCATABLE :: phi(:) ! potential
  !
  REAL(8) :: sigma1, sigma2, dr, x0, y0, z0, dx, dy, dz
  INTEGER :: ip
  REAL(8) :: Uana, Unum
  INTEGER :: N_in
  CHARACTER(56) :: chars_arg
  INTEGER :: iargc
  !
  COMPLEX(8), ALLOCATABLE :: rhoG(:), phiG(:)
  REAL(8) :: UnumG

  IF( iargc() /= 1 ) THEN 
    WRITE(*,*) 'Exactly one argument must be given:', iargc()
    STOP 
  ENDIF 
  CALL getarg(1, chars_arg )
  READ( chars_arg, *) N_in
  WRITE(*,*) 'N_in = ', N_in

  NN(:) = (/ N_in, N_in, N_in /)
  LL = (/ 16.d0, 16.d0, 16.d0 /)
  !
  CALL init_LF3d_p( NN, (/0.d0,0.d0,0.d0/), LL )
  !CALL init_LF3d_c( NN, (/0.d0,0.d0,0.d0/), LL )
  CALL info_LF3d()
 
  ALLOCATE( rho(Npoints) )
  ALLOCATE( phi(Npoints) )

  ! center of the box
  x0 = LL(1)/2.d0
  y0 = LL(2)/2.d0
  z0 = LL(3)/2.d0
  ! Initialize
  sigma1 = 0.75d0
  sigma2 = 0.50d0
  DO ip = 1, Npoints
    dx = lingrid(1,ip) - x0
    dy = lingrid(2,ip) - y0
    dz = lingrid(3,ip) - z0
    dr = sqrt( dx**2 + dy**2 + dz**2 )
    rho(ip) = exp(-dr**2/(2*sigma2**2))/(2*pi*sigma2**2)**1.5d0 - &
              exp(-dr**2/(2*sigma1**2))/(2*pi*sigma1**2)**1.5d0
    !WRITE(*,'(1x,I5,2F18.10)') ip, r, rho(ip)
  ENDDO

  WRITE(*,*)
  WRITE(*,*) 'Integrated rho = ', sum( rho(:) )*dVol

  ! Solve Poisson equation
  CALL Poisson_solve_fft( rho, phi )

  !
  Unum = 0.5d0*sum( rho(:)*phi(:) )*dVol
  Uana = ( (1.d0/sigma1 + 1.d0/sigma2)/2.d0 - sqrt(2.d0)/sqrt(sigma1**2 + sigma2**2) )/sqrt(PI)
  WRITE(*,'(1x,A,F18.10)') 'Unum = ', Unum
  WRITE(*,'(1x,A,F18.10)') 'Uana = ', Uana
  WRITE(*,'(1x,A,E18.10)') 'diff = ', abs(Unum-Uana)

  ALLOCATE( rhoG(Npoints), phiG(Npoints) )
  DO ip = 1,Npoints
    rhoG(ip) = cmplx( rho(ip), 0.d0, kind=8 )
    phiG(ip) = cmplx( phi(ip), 0.d0, kind=8 )
  ENDDO 
  CALL fft_fftw3( rhoG, NN(1), NN(2), NN(3), .false. )
  CALL fft_fftw3( phiG, NN(1), NN(2), NN(3), .false. )

  Unumg = real( sum( rhoG(:)*conjg(phiG(:)) ) )*0.5d0*dVol*Npoints
  WRITE(*,*)
  WRITE(*,*) 'In G space = ', Unumg
  WRITE(*,*) 'Diff = ', abs(Unumg-Uana)

  DEALLOCATE( rho, phi )

  CALL dealloc_ilu0_prec()
  CALL dealloc_nabla2_sparse()
  CALL dealloc_LF3d()
END PROGRAM

