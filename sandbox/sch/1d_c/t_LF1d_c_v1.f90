
! This module contains several global variables and parameter used in this
! program
MODULE m_globals
  USE m_LF1d
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: N = 30
  REAL(8), PARAMETER :: A = -5.d0, B = 5.d0
  !
  TYPE(LF1d_t) :: LF
  !
  ! For plotting purpose
  !
  INTEGER, PARAMETER :: NPTS_PLOT = 400
  !
  REAL(8), ALLOCATABLE :: Vpot(:)
END MODULE


!-------------------------------------
SUBROUTINE init_pot_harmonic(omega, V)
!-------------------------------------
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : N, LF
  IMPLICIT NONE
  ! Arguments
  REAL(8) :: omega
  REAL(8) :: V(N)
  ! local
  INTEGER :: ii

  DO ii=1,N
    V(ii) = 0.5d0*omega**2*LF%grid(ii)**2
  ENDDO
END SUBROUTINE

!------------------------------------
SUBROUTINE plot_pot_harmonic( omega )
!------------------------------------
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : NPTS_PLOT, A, B
  IMPLICIT NONE
  ! Argument
  REAL(8) :: omega
  ! Local
  REAL(8) :: xx, ff
  INTEGER :: ii

  ! Data for plot the potential
  DO ii=1,NPTS_PLOT
    xx = A + (ii-1)*(B-A)/(NPTS_PLOT-1)
    ff = 0.5d0*omega*xx**2
    WRITE(11,*) xx, ff
  ENDDO

END SUBROUTINE


SUBROUTINE test_Etot( evec, N )
  USE m_globals, ONLY : LF, Vpot
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: evec(N,N)
  !
  REAL(8) :: Hpsi(N)
  INTEGER :: i
  REAL(8) :: Ekin, Epot, Etot
  !
  REAL(8) :: ddot
  
  WRITE(*,*)
  WRITE(*,*) 'Calculation of total energy as sum of Ekin and Epot for each states'

  Etot = 0.d0
  DO i = 1, 4
    Hpsi(:) = matmul( LF%D2jl, evec(:,i) )
    Ekin = -0.5d0*ddot( N, evec(:,i),1, Hpsi(:),1 )
    !
    !Epot = sum( Vpot(:)*evec(:,i)**2 ) * LF%h * N / LF%L
    Epot = sum( Vpot(:)*evec(:,i)**2 )  ! h = N/L
    !
    Etot = Etot + Ekin + Epot
    WRITE(*,'(1x,I5,2F18.10)') i, Ekin, Epot
  ENDDO
  WRITE(*,*) 'Etot = ', Etot
END SUBROUTINE


!------------------------------------------------------------------------------
PROGRAM t_LF1d_c
!------------------------------------------------------------------------------
  USE m_constants, ONLY: PI
  USE m_LF1d
  USE m_globals
  IMPLICIT NONE
  !
  INTEGER :: ii, ip, jj
  REAL(8) :: xx, yy
  !
  REAL(8), ALLOCATABLE :: Hamiltonian(:,:)
  REAL(8), ALLOCATABLE :: eval(:)
  ! Parameter for harmonic potential
  REAL(8), PARAMETER :: omega=2.d0

  CALL init_LF1d_c(LF, N, A, B)
  CALL info_LF1d(LF, .FALSE.)
  
  ! Set up potential
  ALLOCATE( Vpot(N) )
  CALL init_pot_harmonic( omega, Vpot )
  CALL plot_pot_harmonic( omega )

  ! Set up the Hamiltonian
  ALLOCATE( Hamiltonian(N,N) )
  !
  Hamiltonian(:,:) = -0.5d0*LF%D2jl(:,:)
  ! Potential, diagonal
  DO ii=1,N
    Hamiltonian(ii,ii) = Hamiltonian(ii,ii) + Vpot(ii)
  ENDDO

  ALLOCATE( eval(N) )
  CALL eig_dsyev( Hamiltonian, eval, N )
  WRITE(*,*)
  WRITE(*,*) 'Eigenvalues:'
  DO ii=1,4
    WRITE(*,'(1x,I5,F18.10)') ii, eval(ii)
  ENDDO
  ! write eigenvectors (coefs)
  DO ii=1,4
    DO ip=1,N
      WRITE(100+ii,*) LF%grid(ip), Hamiltonian(ip,ii)
    ENDDO
  ENDDO

  ! Write to file for plotting, using denser points
  DO ii=1,4
    DO ip = 1, NPTS_PLOT
      !
      xx = A + (ip-1)*LF%L/(NPTS_PLOT-1)
      yy = 0.d0
      DO jj=1,N
        yy = yy + eval_LF1d_c( LF, jj, xx )*Hamiltonian(jj,ii)*sqrt(LF%h)
      ENDDO
      WRITE(200+ii,*) xx, yy
    ENDDO
  ENDDO


  ! The matrix Hamiltonian now carries the eigenvectors
  CALL test_Etot( Hamiltonian, N )

  DEALLOCATE( eval )
  DEALLOCATE( Hamiltonian )
  DEALLOCATE( Vpot )

END PROGRAM



!----------------------------------
SUBROUTINE eig_dsyev(A,eigval,dimA)
!----------------------------------
  IMPLICIT NONE
  ! Arguments
  INTEGER :: dimA
  REAL(8) :: A(dimA,dimA)
  REAL(8) :: eigval(dimA)
  ! Workspace array for ZHEEV
  INTEGER :: lwork
  REAL(8), ALLOCATABLE :: work(:)
  ! Local variables
  INTEGER :: info

  lwork = 3*dimA-1
  ALLOCATE(work(lwork))

  CALL dsyev('v', 'u', dimA, A, dimA, eigval, work, lwork, info)
  IF(info /= 0) THEN
    WRITE(*,*) 'Error on calling dsyev: info=',info
    STOP
  ENDIF

  ! Free memory
  DEALLOCATE(work)
END SUBROUTINE


