! ffr
! First implementation, straightforward approach. Only the initialization
! of the grid and the associated basis function are done in another module

PROGRAM test_LF1dp
  USE m_constants, ONLY: PI
  USE m_LF1d
  !USE LinearAlgebra 
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: N=51
  REAL(8), PARAMETER :: L=5.d0
  TYPE(LF1d_t) :: LF
  INTEGER :: ii, jj
  INTEGER, PARAMETER :: NPTS_PLOT=400
  REAL(8) :: xx, ff
  !
  REAL(8), ALLOCATABLE :: V(:)
  COMPLEX(8), ALLOCATABLE :: Hamiltonian(:,:)
  REAL(8) :: kpt
  REAL(8), ALLOCATABLE :: eval(:)
  ! Parameter for Mathieu
  REAL(8), PARAMETER :: Vm=1.5d0

  CALL init_LF1d_p(LF, N, 0.d0, L)

  ! Set up potential
  ALLOCATE( V(N) )
  DO ii=1,N
    V(ii) = Vm + Vm*cos( 2.d0*PI*LF%grid(ii)/L )
    !V(ii) = 0.d0
  ENDDO

  !WRITE(*,*) 'sum(V)       = ', sum(V)
  !WRITE(*,*) 'sum(LF%D1jl) = ', sum(LF%D1jl)
  !WRITE(*,*) 'D1jl = ', LF%D1jl(2,1)
  !WRITE(*,*) 'sum(LF%D2jl) = ', sum(LF%D2jl)
  !WRITE(*,*) 'D2jl = ', LF%D2jl(2,1), LF%D2jl(2,2)

  ! Data for plot the potential
!  DO ii=-NPTS_PLOT,NPTS_PLOT
!    xx = ii*L/NPTS_PLOT
!    !ff = eval_LF1d_p(LFgrid, 1, xx)
!    ff = Vm + Vm*cos( 2.d0*PI*xx/L )
!    WRITE(11,*) xx, ff
!  ENDDO

  ! Set up the Hamiltonian
  ALLOCATE( Hamiltonian(N,N) )
  Hamiltonian(:,:) = cmplx(0.d0,0.d0)
  kpt = 0.2d0
  !
  Hamiltonian(:,:) = -0.5d0*( LF%D2jl(:,:) + cmplx(0.d0,2.d0)*kpt*LF%D1jl(:,:) )
  ! Diagonal
  DO ii=1,N
    Hamiltonian(ii,ii) = Hamiltonian(ii,ii) + 0.5d0*kpt**2 + V(ii)
  ENDDO

  ALLOCATE( eval(N) )
  CALL eig_zheev( Hamiltonian, eval, N ) 
  DO ii=1,4
    WRITE(*,'(1x,I5,F18.10)') ii, eval(ii)
  ENDDO

  DEALLOCATE( eval )
  DEALLOCATE( Hamiltonian )
  DEALLOCATE( V )

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


