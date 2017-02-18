! efefer, 4 May 2016
!
! Solution of Kohn-Sham equation

! Using iterative (partial) diagonalization

!------------------------------------------------------------------------------
PROGRAM t_LF3d_dft
!------------------------------------------------------------------------------
  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : init_LF3d_sinc, init_LF3d_c, init_LF3d_p
  USE m_globals, ONLY : LF, Npoints, Vpsloc, Rho, A,B,h,N, LF_type, Solution_Method, &
     evecs, evals, Nstate, Focc
  IMPLICIT NONE
  
  A = 0.d0
  B = 6.d0
  N = 22

  WRITE(*,'(/,1x,A)') 'Initializing grids and basis functions:'
  WRITE(*,*)          '---------------------------------------'
  CALL init_LF3d_c( LF, (/N,N,N/), (/A,A,A/), (/B,B,B/) )

  Npoints = N**3

  ! Manually set number of states
  Nstate = 4

  ! Set up potential
  ALLOCATE( Vpsloc(Npoints) )
  ALLOCATE( Rho(Npoints) )
  ALLOCATE( evecs(Npoints,Nstate), evals(Nstate) )
  ! Dont't forget to setup Focc manually
  ALLOCATE( Focc(Nstate) )
  Focc(:) = 2.d0

  CALL init_pot_harmonic( 2.d0, Vpsloc )

  CALL solve_diagonalize()

  CALL calc_rho()

  DEALLOCATE( Focc )
  DEALLOCATE( Rho )
  DEALLOCATE( Vpsloc )
  DEALLOCATE( evecs, evals )

END PROGRAM


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

  r0(:) = 0.d0

  WRITE(*,*)
  WRITE(*,*) 'Initializing harmonic potential'
  WRITE(*,*) '-------------------------------'
  WRITE(*,'(1x,A,F10.5)')  'omega  = ', omega
  WRITE(*,'(1x,A,3F10.5)') 'center = ', r0(:)

  ! N is assumed to be the same as LF$N
  DO ip=1,N**3
    V(ip) = 0.5d0*omega**2 * norm2( LF%lingrid(:,ip) - r0(:) )**2
  ENDDO
  WRITE(*,*) 'DEBUG: sum(V) = ', sum(V)
END SUBROUTINE
