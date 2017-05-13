!! PURPOSE:
!!
!!   This subroutine solves Kohn-Sham equations by minimizing total energy
!!   functional using conjugate gradient algorithm.
!!   The algorithm is based on T.A. Arias notes.
!!
!! AUTHOR:
!!
!!   Fadjar Fathurrahman
!!
!! MODIFIES:
!! 
!!   Global variables `KS_evecs` and `E_total`
!!
SUBROUTINE KS_solve_Emin_cg( alpha_t, NiterMax, restart )

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_states, ONLY : Nstates, &
                       Focc, &
                       v => KS_evecs
  USE m_energies, ONLY : Etot => E_total

  IMPLICIT NONE

  !! Maximum number of iterations
  INTEGER :: NiterMax

  !! Step size taken when evaluating trial direction
  REAL(8) :: alpha_t

  !! If .TRUE. then starting wavefunction will be read from file
  LOGICAL :: restart

  ! Local
  REAL(8), ALLOCATABLE :: g(:,:), g_old(:,:), g_t(:,:)
  REAL(8), ALLOCATABLE :: d(:,:), d_old(:,:)
  REAL(8), ALLOCATABLE :: tv(:,:)
  REAL(8) :: alpha, beta, denum, Etot_old
  !
  INTEGER :: iter
  REAL(8) :: memGB

  ALLOCATE( g(Npoints,Nstates) )
  ALLOCATE( g_old(Npoints,Nstates) )
  ALLOCATE( g_t(Npoints,Nstates) )
  ALLOCATE( d(Npoints,Nstates) )
  ALLOCATE( d_old(Npoints,Nstates) )

  ALLOCATE( tv(Npoints,Nstates) )

  memGB = Npoints*Nstates*5d0 * 8d0 / (1024d0*1024d0*1024.d0)
  WRITE(*,*)
  WRITE(*,*) 'KS_solve_Emin_cg: memGB = ', memGB

  ! Read starting eigenvectors from file
  IF( restart ) THEN
    READ(112) v   ! FIXME Need to use file name
  ENDIF

  CALL calc_Rhoe( Focc, v )
  CALL update_potentials()
  CALL calc_energies( v )

  Etot_old = Etot

  alpha = 0.d0
  beta  = 0.d0

  g(:,:)     = 0.d0
  g_t(:,:)   = 0.d0
  d(:,:)     = 0.d0
  d_old(:,:) = 0.d0

  DO iter = 1, NiterMax
    !
    ! Evaluate gradient at current trial vectors
    CALL calc_grad( Nstates, v, g )
    !
    ! set search direction
    IF( iter /= 1 ) THEN
      ! Fletcher-Reeves
      beta = sum( g * g ) / sum( g_old * g_old )
    ENDIF
    d(:,:) = -g(:,:) + beta*d_old(:,:)
    !
    ! Evaluate gradient at trial step
    tv(:,:) = v(:,:) + alpha_t * d(:,:)
    CALL orthonormalize( Nstates, tv )
    CALL calc_Rhoe( Focc, tv )
    CALL update_potentials()  ! Now global vars on m_hamiltonian are changed
    CALL calc_grad( Nstates, tv, g_t )
    !
    ! Compute estimate of best step and update current trial vectors
    denum = sum( (g - g_t) * d )
    IF( denum /= 0.d0 ) THEN  ! FIXME: use abs ?
      alpha = abs( alpha_t * sum( g * d )/denum )
    ELSE 
      alpha = 0.d0
    ENDIF
    !WRITE(*,*) 'iter, alpha_t, alpha, beta', iter, alpha_t, alpha, beta

    v(:,:) = v(:,:) + alpha * d(:,:)
    CALL orthonormalize( Nstates, v )
    CALL calc_Rhoe( Focc, v )
    CALL update_potentials()

    CALL calc_energies( v )

    WRITE(*,'(1x,I5,F18.10,E18.10)') iter, Etot, abs(Etot-Etot_old)
    !
    IF( abs(Etot - Etot_old) < 1.d-7 ) THEN
      WRITE(*,*) 'KS_solve_Emin_cg converged in iter', iter
      EXIT
    ENDIF
    !
    Etot_old = Etot
    g_old(:,:) = g(:,:)
    d_old(:,:) = d(:,:)
  ENDDO

  WRITE(111) v

  DEALLOCATE( g, g_old, g_t, d, d_old, tv )
END SUBROUTINE

