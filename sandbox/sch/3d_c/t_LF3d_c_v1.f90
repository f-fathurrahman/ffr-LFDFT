! efefer, 1 January 2016
!
! Solution of Schrodinger equation
! Using full Hamiltonian matrix, it will take a lot of memory for 
! large number of basis set.


! This module contains several global variables and parameter used in this
! program
MODULE m_globals
  USE m_LF3d
  IMPLICIT NONE
  ! These parameters are similar for x, y, and z directions
  INTEGER :: N
  REAL(8) :: A = 0.d0, B = 6.d0
  !
  TYPE(LF3d_t) :: LF
END MODULE


!-------------------------------------
SUBROUTINE init_pot_harmonic(omega, V)
!-------------------------------------
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : N, LF, A, B
  IMPLICIT NONE
  ! Arguments
  REAL(8) :: omega
  REAL(8) :: V(N**3)
  ! local
  INTEGER :: ip
  REAL(8) :: r0(3)

  r0(:) = (B-A)/2.d0

  ! N is assumed to be the same as LF$N
  DO ip=1,N**3
    V(ip) = 0.5d0*omega**2 * norm2( LF%lingrid(:,ip) - r0(:) )**2
  ENDDO
END SUBROUTINE


!------------------------------------------------------------------------------
PROGRAM t_LF3d_c
!------------------------------------------------------------------------------
  USE m_constants, ONLY: PI
  USE m_LF3d
  USE m_globals
  IMPLICIT NONE
  !
  INTEGER :: i1,i2, j1,j2, k1,k2
  INTEGER :: ip1, ip2, ip, ii
  !
  REAL(8), ALLOCATABLE :: V(:)
  REAL(8), ALLOCATABLE :: Hamiltonian(:,:)
  REAL(8), ALLOCATABLE :: eval(:)
  ! Parameter for harmonic potential
  REAL(8), PARAMETER :: omega=2.d0
  !

  N = 11
  A = 0.d0
  B = 6.d0

  CALL init_LF3d_c(LF, (/N,N,N/), (/A,A,A/), (/B,B,B/) )

  CALL write_matrix_hdf5( LF%LFx%D2jl, N, N, 'D2jl.h5', 'D2jl' )
  
  ! Set up potential
  ALLOCATE( V(N**3) )
  CALL init_pot_harmonic( omega, V )

  ! Set up the Hamiltonian
  ! WARNING: the size of the Hamiltonian is N**3
  ALLOCATE( Hamiltonian(N**3,N**3) )
  ! Kinetic part
  DO ip2 = 1, N**3
    DO ip1 = ip2, N**3
      i1 = LF%lin2xyz(1,ip1)
      i2 = LF%lin2xyz(1,ip2)
      !
      j1 = LF%lin2xyz(2,ip1)
      j2 = LF%lin2xyz(2,ip2)
      !
      k1 = LF%lin2xyz(3,ip1)
      k2 = LF%lin2xyz(3,ip2)
      !
      Hamiltonian(ip1,ip2) = 0.d0
      !
      IF( j1 == j2 .AND. k1 == k2 ) THEN
        Hamiltonian(ip1,ip2) = Hamiltonian(ip1,ip2) + LF%LFx%D2jl(i1,i2)
      ENDIF
      !
      IF( i1 == i2 .AND. k1 == k2 ) THEN
        Hamiltonian(ip1,ip2) = Hamiltonian(ip1,ip2) + LF%LFy%D2jl(j1,j2)
      ENDIF
      !
      IF( i1 == i2 .AND. j1 == j2 ) THEN
        Hamiltonian(ip1,ip2) = Hamiltonian(ip1,ip2) + LF%LFz%D2jl(k1,k2)
      ENDIF
      !
      Hamiltonian(ip2,ip1) = Hamiltonian(ip1,ip2)
      !WRITE(999,'(2I8,F18.10)') ip1, ip2, Hamiltonian(ip1,ip2)
    ENDDO
  ENDDO
  !
  Hamiltonian(:,:) = -0.5d0*Hamiltonian(:,:)
  !WRITE(113) Hamiltonian
  CALL write_matrix_hdf5( Hamiltonian, N**3, N**3, 'Kinetic.h5', 'Kinetic' )
  !
  ! Potential, diagonal
  !
  DO ip=1,N**3
    Hamiltonian(ip,ip) = Hamiltonian(ip,ip) + V(ip)
  ENDDO
  !WRITE(114) Hamiltonian
  CALL write_matrix_hdf5( Hamiltonian, N**3, N**3, 'Hamiltonian.h5', 'Hamiltonian' )

  ALLOCATE( eval(N**3) )
  !CALL eig_dsyev( Hamiltonian, eval, N**3 )
  WRITE(*,*)
  WRITE(*,*) 'Eigenvalues:'
  DO ii=1,4
    WRITE(*,'(1x,I5,F18.10)') ii, eval(ii)
  ENDDO

  DEALLOCATE( eval )
  DEALLOCATE( Hamiltonian )
  DEALLOCATE( V )

END PROGRAM



!----------------------------------
SUBROUTINE eig_dsyev(A,eigval,dimA)
!----------------------------------
  IMPLICIT NONE
  ! Arguments
  INTEGER :: dimA
  REAL(8) :: A(dimA,dimA)
  REAL(8) :: eigval(dimA)
  ! Workspace array for DSYEV
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


