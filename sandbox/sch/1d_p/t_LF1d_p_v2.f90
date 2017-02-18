! ffr
! Second implementation: create a more modular program


! This module contains several global variables and parameter used in this
! program
MODULE m_globals
  USE m_LF1d
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: L
  !
  TYPE(LF1d_t) :: LFgrid
  !
  INTEGER :: Nstates
  !
  ! For plotting purpose
  !
  INTEGER :: NPTS_PLOT
END MODULE


!---------------------------------
SUBROUTINE init_pot_mathieu(Vm, V)
!---------------------------------
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : N, L, LFgrid
  IMPLICIT NONE
  ! Arguments
  REAL(8) :: Vm
  REAL(8) :: V(N)
  ! local
  INTEGER :: ii

  DO ii=1,N
    V(ii) = Vm + Vm*cos( 2.d0*PI*LFgrid%grid(ii)/L )
  ENDDO
END SUBROUTINE

!--------------------------------
SUBROUTINE plot_pot_mathieu( Vm )
!--------------------------------
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : NPTS_PLOT, L
  IMPLICIT NONE
  ! Argument
  REAL(8) :: Vm
  ! Local
  REAL(8) :: xx, ff
  INTEGER :: ii
  INTEGER :: iu

  iu = 11
  OPEN(unit=iu,file='V_mathieu',action='write')

  ! Data for plot the potential
  DO ii=0,NPTS_PLOT
    xx = ii*L/NPTS_PLOT
    ff = Vm + Vm*cos( 2.d0*PI*xx/L )
    WRITE(iu,*) xx, ff
  ENDDO
  CLOSE(iu)

END SUBROUTINE


SUBROUTINE plot_evec( evec )
  USE m_globals, ONLY: N, Nstates, LFgrid
  IMPLICIT NONE
  !
  COMPLEX(8) :: evec(N,N)
  !
  INTEGER :: iur, iuc, is, ip

  iur = 12
  iuc = 13
  OPEN(unit=iur,file='evec_real.dat',action='write')
  OPEN(unit=iuc,file='evec_cplx.dat',action='write')
  DO ip=1,N
    WRITE(iur,'(F18.10)',advance='no') LFgrid%grid(ip)
    WRITE(iuc,'(F18.10)',advance='no') LFgrid%grid(ip)
    DO is=1,Nstates
      WRITE(iur,'(F18.10)',advance='no') real( evec(ip,is) )
      WRITE(iuc,'(F18.10)',advance='no') imag( evec(ip,is) )
    ENDDO
    WRITE(iur,*)
    WRITE(iuc,*)
  ENDDO
  CLOSE(iur)
  CLOSE(iuc)
END SUBROUTINE

SUBROUTINE plot_rho( rho )
  USE m_globals, ONLY : N, LFgrid
  IMPLICIT NONE
  !
  REAL(8) :: rho(N)
  !
  INTEGER :: ip, iu

  iu = 15
  OPEN(unit=iu,file='rho.dat',action='write')
  DO ip=1,N
    WRITE(iu,'(2F18.10)') LFgrid%grid(ip), rho(ip)
  ENDDO
  CLOSE(iu)
END SUBROUTINE


!------------------------------------------------------------------------------
PROGRAM t_LF1dp
!------------------------------------------------------------------------------
  USE m_constants, ONLY: PI
  USE m_LF1d
  USE m_globals
  IMPLICIT NONE
  !
  INTEGER :: ii
  !
  REAL(8), ALLOCATABLE :: V(:)
  COMPLEX(8), ALLOCATABLE :: Hamiltonian(:,:)
  REAL(8) :: kpt
  REAL(8), ALLOCATABLE :: eval(:)
  REAL(8), ALLOCATABLE :: rho(:)
  ! Parameter for Mathieu
  REAL(8), PARAMETER :: Vm=1.5d0
  !
  COMPLEX(8) :: zdotc

  N = 51
  L = 5.d0
  NPTS_PLOT = 400
  Nstates = 4

  CALL init_LF1d_p(LFgrid, N, 0.d0, L)
  
  ! Set up potential
  ALLOCATE( V(N) )
  CALL init_pot_mathieu( Vm, V )
  CALL plot_pot_mathieu( Vm )

  ! Set up the Hamiltonian
  ALLOCATE( Hamiltonian(N,N) )
  Hamiltonian(:,:) = cmplx(0.d0,0.d0)
  kpt = 0.2d0
  !
  Hamiltonian(:,:) = -0.5d0*( LFgrid%D2jl(:,:) + cmplx(0.d0,2.d0)*kpt*LFgrid%D1jl(:,:) )
  ! Diagonal
  DO ii=1,N
    Hamiltonian(ii,ii) = Hamiltonian(ii,ii) + 0.5d0*kpt**2 + V(ii)
  ENDDO

  ALLOCATE( eval(N) )
  CALL eig_zheev( Hamiltonian, eval, N )
  WRITE(*,*) 'Eigenvalues:'
  DO ii=1,Nstates
    WRITE(*,'(1x,I5,F18.10)') ii, eval(ii)
  ENDDO
  
  ALLOCATE( rho(N) )
  rho(:) = 0.d0
  ! Test the normalization of eigenvectors and calculate the electron density
  ! Note that the array Hamiltonian now contains eigenvectors
  WRITE(*,*) 'Norm of eigenvectors:'
  DO ii=1,Nstates
    WRITE(*,'(1x,I5,F18.10)') ii, sqrt( real(zdotc( N, Hamiltonian(:,ii),1, Hamiltonian(:,ii),1 )) )
    rho(:) = rho(:) + conjg(Hamiltonian(:,ii))*Hamiltonian(:,ii)
  ENDDO
  WRITE(*,*) 'integrated rho:', sum(rho(:))

  !CALL plot_evec( Hamiltonian )
  !CALL plot_rho( rho )

  DEALLOCATE( eval )
  DEALLOCATE( Hamiltonian )
  DEALLOCATE( V )
  DEALLOCATE( rho )

END PROGRAM



!----------------------------------------
SUBROUTINE eig_zheev(A,eigval,dimA)
!----------------------------------------
  IMPLICIT NONE
  ! Arguments
  INTEGER :: dimA
  COMPLEX(8) :: A(dimA,dimA)
  REAL(8) :: eigval(dimA)
  ! Workspace array for ZHEEV
  INTEGER :: lwork
  COMPLEX(8), ALLOCATABLE :: work(:)
  COMPLEX(8), ALLOCATABLE :: rwork(:)
  ! Local variables
  INTEGER :: info

  lwork = 2*dimA-1
  ALLOCATE(work(lwork))
  ALLOCATE(rwork(3*dimA-2))

  CALL zheev('v', 'u', dimA, A, dimA, eigval, work, lwork, rwork, info)
  IF(info /= 0) THEN
    WRITE(*,*) 'Error on calling zheev: info=',info
    STOP
  ENDIF

  ! Free memory
  DEALLOCATE(work)
  DEALLOCATE(rwork)
END SUBROUTINE


