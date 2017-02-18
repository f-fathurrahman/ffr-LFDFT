!---------------------------------------------------
SUBROUTINE minimE_pcg( alpha_t, Niter, restart )
!---------------------------------------------------
  USE m_globals, ONLY : N, evecs, Ncol => Nstate, eVtau, eVtau_old
  IMPLICIT NONE
  !
  INTEGER :: Niter
  REAL(8) :: alpha_t  ! step size
  LOGICAL :: restart
  ! TODO: make these variables allocatable instead of automatic in order
  !       to avoid stack error
  REAL(8), ALLOCATABLE :: grad(:,:), grad_old(:,:), grad_t(:,:)
  REAL(8), ALLOCATABLE :: Kgrad(:,:), Kgrad_old(:,:)
  REAL(8), ALLOCATABLE :: dir(:,:), dir_old(:,:)
  REAL(8) :: Etot, Ekin(Ncol), Epot, Etot_old, norm_grad, deltaE
  REAL(8) :: alpha(Ncol), beta(Ncol)
  !
  INTEGER :: ic, iter, iterTrue
  !
  REAL(8) :: ddot

  CALL init_evalsT()

  ALLOCATE( grad(N**3,Ncol), grad_old(N**3,Ncol), grad_t(N**3,Ncol) )
  ALLOCATE( dir(N**3,Ncol), dir_old(N**3,Ncol) )
  ALLOCATE( Kgrad(N**3,Ncol), Kgrad_old(N**3,Ncol) )

  ! initialize v
  IF( .NOT. restart ) THEN
    DO ic = 1, Ncol
      CALL r8_rand_vec( N**3, evecs(:,ic) )
    ENDDO
    CALL ortho_gram_schmidt( evecs, N**3, N**3, Ncol )
  ELSE
    READ(112) evecs
    ! No need to orthonormalize, v should already be ortonormal
  ENDIF

  CALL get_Etot( Ncol, evecs, Ekin, Epot, Etot)
  Etot_old = Etot

  iterTrue = 0
  iter = 0
  cgloop: DO
    iter = iter + 1
    iterTrue = iterTrue + 1
    !
    ! Evaluate gradient at current trial vectors
    !
    CALL get_grad( Ncol, evecs, grad, .FALSE. )
    !
    ! Calculate norm of the gradient, not really used in the minimization
    ! procedure.
    !
    norm_grad = 0.d0
    DO ic = 1, Ncol
      norm_grad = norm_grad + norm2( grad(:,ic) )
    ENDDO
    norm_grad = norm_grad/Ncol
    !
    ! set search direction
    !
    IF( iter == 1 ) THEN
      DO ic=1,Ncol
        ! always do preconditioning at the first step
        CALL apply_PrecCG( N**3, grad(:,ic), Kgrad(:,ic), 50, ic )
        !CALL avgPrec( grad(:,ic), Kgrad(:,ic) )
        !CALL apply_PrecEvals( grad(:,ic), Kgrad(:,ic), ic )
        dir(:,ic) = -Kgrad(:,ic)
      ENDDO
    ELSE
      DO ic=1,Ncol
        !IF( deltaE > 1.d-5 ) THEN
          CALL apply_PrecCG( N**3, grad(:,ic), Kgrad(:,ic), 50, ic )
          !CALL avgPrec( grad(:,ic), Kgrad(:,ic) )
          !CALL apply_PrecEvals( grad(:,ic), Kgrad(:,ic), ic )
        !ELSE
        !  Kgrad = grad
        !ENDIF
        beta(ic) = ddot( N**3, grad(:,ic),1, Kgrad(:,ic),1 ) / &
                   ddot( N**3, grad_old(:,ic),1, Kgrad_old(:,ic),1 )
        !beta(ic) = abs( beta(ic) )
        WRITE(*,'(1x,A,F18.10)') 'beta(ic) = ', beta(ic)
        dir(:,ic) = -Kgrad(:,ic) + beta(ic)*dir_old(:,ic)
      ENDDO
    ENDIF
    !
    ! Evaluate gradient at trial step
    !
    CALL get_grad( Ncol, evecs + alpha_t*dir, grad_t, .TRUE. )
    !
    ! Compute estimate of best step and update current trial vectors
    !
    DO ic = 1, Ncol
      alpha(ic) = alpha_t*ddot( N**3, grad(:,ic),1, dir(:,ic),1 )/&
            ddot( N**3, grad(:,ic)-grad_t(:,ic),1, dir(:,ic), 1 )
      alpha(ic) = abs( alpha(ic) ) ! force alpha to be positive?
      WRITE(*,'(1x,A,F18.10)') 'alpha(ic) = ', alpha(ic)
      evecs(:,ic) = evecs(:,ic) + alpha(ic)*dir(:,ic)
    ENDDO
    !
    CALL ortho_gram_schmidt( evecs, N**3, N**3, Ncol )
    CALL get_Etot_v2( Ncol, evecs, Ekin, Epot, Etot )
    !
    deltaE = Etot_old - Etot
    IF( deltaE < 0.d0 ) THEN
      WRITE(*,*) 
      WRITE(*,*) '*****Negative deltaE, restarting CG:', deltaE
      !iter = 0
      !CYCLE cgloop
    ENDIF
    WRITE(*,*)
    WRITE(*,'(1x,A,I8,2F18.10)') 'iter, Etot, deltaE ', iterTrue, Etot, Etot_old-Etot
    WRITE(*,'(1x,A,F18.10)') 'Epot: ', Epot
    WRITE(*,*) 'Ekin(ic)'
    DO ic = 1, Ncol
      WRITE(*,'(1x,I8,F18.10)') ic, Ekin(ic)
      eVtau(ic) = Ekin(ic)
    ENDDO
    !eVtau_old(:) = eVtau(:)
    !DO ic = 1, Ncol
    !  IF( eVtau(ic) < evTau_old(ic) ) THEN
    !    eVtau(ic) = eVtau_old(ic)
    !  ENDIF
    !ENDDO
    !
    IF( abs(Etot - Etot_old) < 1.d-6 ) THEN
      WRITE(*,*) 'PCG converged in iter', iter
      EXIT
    ENDIF
    !
    Etot_old = Etot
    grad_old(:,:) = grad(:,:)
    dir_old(:,:) = dir(:,:)
    Kgrad_old(:,:) = Kgrad(:,:)
    !
    IF( iter >= Niter ) THEN
      EXIT cgloop
    ENDIF
  ENDDO cgloop

  WRITE(111) evecs

  DEALLOCATE( grad, grad_old, grad_t, dir, dir_old, Kgrad, Kgrad_old )
END SUBROUTINE
