! efefer, 15 January 2016
!
! Solution of Schrodinger equation

! Using iterative (partial) diagonalization or total energy minimization
! via conjugate gradient


!------------------------------------------------------------------------------
PROGRAM t_LF3d_sch
!------------------------------------------------------------------------------
  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : init_LF3d_sinc, init_LF3d_c, init_LF3d_p
  USE m_globals, ONLY : LF, Npoints, Vpsloc, Rho, A,B,h,N, LF_type, Solution_Method, &
     evecs, evals, Nstate, Focc
  IMPLICIT NONE

  CALL read_arguments()

  WRITE(*,'(/,1x,A)') 'Initializing grids and basis functions:'
  WRITE(*,*)          '---------------------------------------'
  IF( LF_type == 'sinc' ) THEN
    CALL init_LF3d_sinc(LF, (/N,N,N/), (/h,h,h/) )
    A = LF%LFx%A
    B = LF%LFx%B
  ELSEIF( LF_type == 'box' ) THEN
    CALL init_LF3d_c( LF, (/N,N,N/), (/A,A,A/), (/B,B,B/) )
  ELSEIF( LF_type == 'per' ) THEN
    ! TODO: This is not yet supported
    CALL init_LF3d_p( LF, (/N,N,N/), (/A,A,A/), (/B,B,B/) )
  ENDIF
  
  ! FIXME
  Npoints = N**3

  ! Manually set number of states
  Nstate = 1

  ! Set up potential
  ALLOCATE( Vpsloc(Npoints) )
  ALLOCATE( Rho(Npoints) )
  ALLOCATE( evecs(Npoints,Nstate), evals(Nstate) )
  ! Dont't forget to setup Focc manually
  ALLOCATE( Focc(Nstate) )
  Focc(:) = 2.d0

  CALL init_pspot_H2_hgh( Vpsloc )

  IF( Solution_Method == 'diag') THEN
    CALL solve_diagonalize()
  ELSEIF( Solution_Method == 'cg' .OR. Solution_Method == 'pcg' ) THEN
    CALL solve_minim()
  ELSEIF( Solution_Method == 'MG' ) THEN
    CALL solve_MG() ! FIXME not implemented yet
  ELSE
    WRITE(*,*) 'ERROR: Unrecognized method: ', trim(Solution_Method)
    STOP
  ENDIF

  CALL calc_rho()

  DEALLOCATE( Focc )
  DEALLOCATE( Rho )
  DEALLOCATE( Vpsloc )
  DEALLOCATE( evecs, evals )

END PROGRAM


!-------------------------------------
SUBROUTINE init_pspot_H2_hgh( V )
!-------------------------------------
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : N, LF, A, B
  IMPLICIT NONE
  ! Arguments
  REAL(8) :: V(N**3)
  ! local
  INTEGER, PARAMETER :: NATOM=2
  INTEGER :: ip, ia
  REAL(8) :: r0(3,NATOM)
  !
  REAL(8) :: r, Zion, rloc, C1, C2, rrloc


  r0(:,1) = (/ -0.6991d0, 0.d0, 0.d0 /)
  r0(:,2) = (/ +0.6991d0, 0.d0, 0.d0 /)

  Zion = 1.d0
  rloc = 0.2d0
  C1 = -4.0663326d0
  C2 = 0.6678322d0

  WRITE(*,*)
  WRITE(*,*) 'Initializing HGH pseudopotential for H2 molecule'
  WRITE(*,*) '------------------------------------------------'
  WRITE(*,'(1x,A,3F10.5)') 'pos1 = ', r0(:,1)
  WRITE(*,'(1x,A,3F10.5)') 'pos2 = ', r0(:,2)

  ! TODO Add journal reference
  V(:) = 0.d0
  DO ia = 1, NATOM
    DO ip = 1, N**3
      r = norm2( LF%lingrid(:,ip)-r0(:,ia) ) + 1.d-8
      rrloc = r/rloc
      V(ip) = V(ip) + -Zion/r * erf(r/(sqrt(2.d0)*rloc)) + &
              (C1 + C2*rrloc**2)*exp(-0.5d0*(rrloc)**2)
    ENDDO
  ENDDO
  WRITE(*,*) 'DEBUG: sum(V) = ', sum(V)
END SUBROUTINE


