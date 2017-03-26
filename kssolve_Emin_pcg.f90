SUBROUTINE kssolve_Emin_pcg( alpha_t, Niter, restart )

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_states, ONLY : Nstates, &
                       Focc, &
                       v => KS_evecs
  USE m_energies, ONLY : Etot => E_total
  IMPLICIT NONE
  !
  INTEGER :: Niter
  REAL(8) :: alpha_t  ! step size
  LOGICAL :: restart
  REAL(8), ALLOCATABLE :: g(:,:), g_old(:,:), g_t(:,:)
  REAL(8), ALLOCATABLE :: d(:,:), d_old(:,:)
  REAL(8), ALLOCATABLE :: Kg(:,:), Kg_old(:,:) ! preconditioned
  REAL(8), ALLOCATABLE :: tv(:,:)
  REAL(8) :: alpha, beta, denum, Etot_old
  !
  INTEGER :: iter, ist
  REAL(8) :: memGB

  WRITE(*,*) 'Pass here 22'

  ALLOCATE( g(Npoints,Nstates) )
  ALLOCATE( g_old(Npoints,Nstates) )
  ALLOCATE( g_t(Npoints,Nstates) )
  ALLOCATE( d(Npoints,Nstates) )
  ALLOCATE( d_old(Npoints,Nstates) )

  ALLOCATE( Kg(Npoints,Nstates) )
  ALLOCATE( Kg_old(Npoints,Nstates) )

  ALLOCATE( tv(Npoints,Nstates) )

  memGB = Npoints*Nstates*8d0 * 8d0 / (1024d0*1024d0*1024.d0)
  WRITE(*,*) 'memGB = ', memGB

  ! Read starting eigenvectors from file
  IF( restart ) THEN
    READ(112) v   ! FIXME Need to use file name
  ENDIF

  CALL calc_rhoe( v, Focc )
  CALL update_potentials()
  CALL calc_energies( v )
  WRITE(*,*) 'Initial Etot = ', Etot

  Etot_old = Etot

  alpha = 0.d0
  beta  = 0.d0

  g(:,:)     = 0.d0
  g_t(:,:)   = 0.d0
  d(:,:)     = 0.d0
  d_old(:,:) = 0.d0
  Kg(:,:)    = 0.d0
  Kg_old(:,:) = 0.d0

  DO iter = 1, Niter
    WRITE(*,*) 'Iter = ', iter
    !
    ! Evaluate gradient at current trial vectors
    CALL calc_grad( Nstates, v, g )
    ! Precondition
    DO ist = 1, Nstates
      !WRITE(*,*) 'ist = ', ist
      CALL linsolve_H( g(:,ist), Kg(:,ist), 100 )
    ENDDO
    WRITE(*,*) 'Pass here 69'
    !
    ! set search direction
    IF( iter /= 1 ) THEN
      ! Fletcher-Reeves
      beta = sum( g * Kg ) / sum( g_old * Kg_old )
    ENDIF
    d(:,:) = -Kg(:,:) + beta*d_old(:,:)
    !
    ! Evaluate gradient at trial step
    tv(:,:) = v(:,:) + alpha_t * d(:,:)
    CALL orthonormalize( Nstates, tv )
    CALL calc_rhoe( tv, Focc )
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
    CALL calc_rhoe( v, Focc )
    CALL update_potentials()

    CALL calc_energies( v )
    !
    !WRITE(*,'(/,1x,A,I5,2F18.10)') 'iter, conv, ||grad||: ', iter, abs(Etot-Etot_old), norm_grad
    !WRITE(*,'(1x,A,3F18.10)') 'Ekin, Epot, Etot: ', Ekin, Epot, Etot
    WRITE(*,'(1x,I5,F18.10,E18.10)') iter, Etot, abs(Etot-Etot_old)
    !
    IF( abs(Etot - Etot_old) < 1.d-7 ) THEN
      WRITE(*,*) 'kssolve_Emin_pcg converged in iter', iter
      EXIT
    ENDIF
    !
    Etot_old = Etot
    g_old(:,:) = g(:,:)
    d_old(:,:) = d(:,:)
    Kg_old(:,:) = Kg(:,:)
  ENDDO

  WRITE(111) v

  DEALLOCATE( g, g_old, g_t, d, d_old, tv, Kg, Kg_old )
END SUBROUTINE

