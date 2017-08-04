MODULE m_Poisson_solve_ISF
  IMPLICIT NONE 

  REAL(8), ALLOCATABLE :: karray(:,:,:)
  INTEGER :: ISF_order
  INTEGER :: nfft1, nfft2, nfft3
  INTEGER :: n1k, n2k, n3k
  REAL(8), ALLOCATABLE :: rhopot(:,:,:)
  REAL(8) :: Ehartree
END MODULE 

SUBROUTINE init_Poisson_solve_ISF()
  USE m_LF3d, ONLY : NN => LF3d_NN, &
                     hh => LF3d_hh
  USE m_Poisson_solve_ISF, ONLY : nfft1, nfft2, nfft3, &
                                  n1k, n2k, n3k, &
                                  karray, rhopot, ISF_order
  IMPLICIT NONE 
  REAL(8) :: hgrid

  CALL Dimensions_FFT( NN(1),NN(2),NN(3), nfft1, nfft2, nfft3 )
  n1k = nfft1/2 + 1
  n2k = nfft2/2 + 1
  n3k = nfft3/2 + 1
  WRITE(*,*) 'nfft       : ', nfft1, nfft2, nfft3
  WRITE(*,*) 'Kernel size: ', n1k, n2k, n3k

  ALLOCATE( karray(n1k,n2k,n3k) )
  
  ISF_order = 14 ! 8, 14, 16
  hgrid = hh(1)
  CALL Build_Kernel( NN(1), NN(2), NN(3), nfft1,nfft2,nfft3, hgrid, ISF_order, karray)

  ALLOCATE( rhopot(NN(1),NN(2),NN(3)) )
END SUBROUTINE 


! The subroutine init_Poisson_solve_ISF must be called before
SUBROUTINE Poisson_solve_ISF( rho, phi )
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     lin2xyz => LF3d_lin2xyz, &
                     hh => LF3d_hh, &
                     NN => LF3d_NN
  USE m_Poisson_solve_ISF, ONLY : karray, rhopot, Ehartree, &
                                  nfft1, nfft2, nfft3, karray

  IMPLICIT NONE 
  REAL(8) :: rho(Npoints)
  REAL(8) :: phi(Npoints)
  !
  INTEGER :: ip, i, j, k
  REAL(8) :: hgrid

  DO ip = 1,Npoints
    i = lin2xyz(1,ip)
    j = lin2xyz(2,ip)
    k = lin2xyz(3,ip)
    rhopot(i,j,k) = rho(ip)
  ENDDO 

  hgrid = hh(1)
  CALL PSolver_Kernel( NN(1), NN(2), NN(3), nfft1, nfft2, nfft3, hgrid, karray, &
                       rhopot, Ehartree)

  DO ip = 1,Npoints
    i = lin2xyz(1,ip)
    j = lin2xyz(2,ip)
    k = lin2xyz(3,ip)
    phi(ip) = rhopot(i,j,k)
  ENDDO 

END SUBROUTINE 

