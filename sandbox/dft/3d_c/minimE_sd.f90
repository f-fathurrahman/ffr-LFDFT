!-------------------------------------------------
SUBROUTINE minimE_sd( alpha, Niter, restart )
!-------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : N, Nstates, Nbasis, deltaV
  USE m_globals, ONLY : evecs, Focc, Rho
  USE m_globals, ONLY : Etot, Ekin, Epot, Ehartree, Exc
  USE m_globals, ONLY : Vhartree, Vxc
  IMPLICIT NONE
  !
  INTEGER :: Niter
  REAL(8) :: alpha  ! step size
  LOGICAL :: restart
  !
  REAL(8), ALLOCATABLE :: grad(:,:)
  REAL(8) :: Etot_old, norm_grad
  REAL(8), ALLOCATABLE :: epsxc(:), depsxc(:)
  !
  INTEGER :: ic, iter

  ALLOCATE( grad(Nbasis,Nstates) )
  ALLOCATE( epsxc(Nbasis) )
  ALLOCATE( depsxc(Nbasis) )

  ! initialize v
  IF( .NOT. restart ) THEN
    DO ic = 1, Nstates
      CALL r8_rand_vec( Nbasis, evecs(:,ic) )
    ENDDO
  ELSE
    READ(112) evecs
  ENDIF
  !
  CALL ortho_gram_schmidt( evecs, Nbasis, Nbasis, Nstates )

  CALL get_Rho( Nbasis, Nstates, Focc, evecs, Rho )

  CALL solve_poisson_cg( Nbasis, -4.d0*PI*Rho(:), Vhartree )
  Ehartree = 0.5d0*sum( Rho(:) * Vhartree(:) ) * deltaV
  CALL excVWN( Nbasis, Rho, epsxc )
  Exc      = sum( Rho(:) * epsxc(:) ) * deltaV
  CALL excpVWN( Nbasis, Rho, depsxc )
  Vxc      = epsxc(:) + Rho(:) * depsxc(:)
  !
  CALL get_Etot( Nbasis, Nstates, evecs, Ekin, Epot, Etot)
  !
  WRITE(*,*) 'Initial Etot = ', Etot
  Etot_old = Etot

  DO iter = 1, Niter
    CALL get_grad( Nstates, evecs, grad, .FALSE. )
    !
    norm_grad = 0.d0
    DO ic = 1, Nstates
      norm_grad = norm_grad + norm2( grad(:,ic) )
    ENDDO
    norm_grad = norm_grad/Nstates
    !
    evecs = evecs - alpha*grad
    !
    CALL ortho_gram_schmidt( evecs, Nbasis, Nbasis, Nstates )
    CALL get_Rho( Nbasis, Nstates, Focc, evecs, Rho )
    !
    CALL solve_poisson_cg( Nbasis, -4.d0*PI*Rho(:), Vhartree )
    Ehartree = 0.5d0*sum( Rho(:) * Vhartree(:) ) * deltaV
    CALL excVWN( Nbasis, Rho, epsxc )
    Exc      = sum( Rho(:) * epsxc(:) ) * deltaV
    CALL excpVWN( Nbasis, Rho, depsxc )
    Vxc      = epsxc(:) + Rho(:) * depsxc(:)
    !
    CALL get_Etot( Nbasis, Nstates, evecs, Ekin, Epot, Etot )
    !
    WRITE(*,'(1x,A,I5,2F18.10)') 'iter, conv, ||grad||: ', iter, abs(Etot-Etot_old), norm_grad
    CALL print_energies()
    ! Convergence criteria in hartree
    IF( abs(Etot - Etot_old) < 1.d-7 ) THEN
      WRITE(*,*) 'SD converged in iter', iter
      EXIT
    ENDIF
    Etot_old = Etot
  ENDDO

  WRITE(111) evecs

  DEALLOCATE( epsxc )
  DEALLOCATE( depsxc )
  DEALLOCATE( grad )

END SUBROUTINE

