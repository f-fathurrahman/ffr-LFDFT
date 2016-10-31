
!---------------------------------------------------
SUBROUTINE minimE_cg( alpha_t, Niter, restart )
!---------------------------------------------------
  USE m_globals, ONLY : N, evecs, Ncol => Nstate
  IMPLICIT NONE
  !
  INTEGER :: Niter
  REAL(8) :: alpha_t  ! step size
  LOGICAL :: restart
  ! TODO: make these variables allocatable instead of automatic in order
  !       to avoid stack error
  REAL(8), ALLOCATABLE :: grad(:,:), grad_old(:,:), grad_t(:,:)
  REAL(8), ALLOCATABLE :: dir(:,:), dir_old(:,:)
  REAL(8) :: Etot, Ekin(Ncol), Epot, Etot_old, norm_grad, deltaE
  REAL(8) :: alpha(Ncol), beta(Ncol)
  !
  INTEGER :: ic, iter, iterTrue
  !
  REAL(8) :: ddot

  ALLOCATE( grad(N**3,Ncol), grad_old(N**3,Ncol), grad_t(N**3,Ncol) )
  ALLOCATE( dir(N**3,Ncol), dir_old(N**3,Ncol) )

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

  iter = 0
  iterTrue = 0
  cgstep: DO
    iter = iter + 1
    iterTrue = iterTrue + 1
    !
    ! Evaluate gradient at current trial vectors
    !
    CALL get_grad( Ncol, evecs, grad, .FALSE. )
    !CALL write_grad( N**3, Ncol, grad, iter )
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
      dir(:,:) = -grad(:,:)
    ELSE
      DO ic=1,Ncol
        beta(ic) = ddot( N**3, grad(:,ic),1, grad(:,ic),1 ) / &
                   ddot( N**3, grad_old(:,ic),1, grad_old(:,ic),1 )
        !WRITE(*,'(1x,A,F18.10)') 'beta(ic) = ', beta(ic)
        dir(:,ic) = -grad(:,ic) + beta(ic)*dir_old(:,ic)
      ENDDO
    ENDIF
      ! Polak-Ribiere, not as efficient as Fletcher-Reeves?
      !DO ic=1,Ncol
      !  beta(ic) = ddot( N**3, grad(:,ic)-grad_old(:,ic),1, grad(:,ic),1 ) / &
      !             ddot( N**3, grad_old(:,ic),1, grad_old(:,ic),1 )
      !  dir(:,ic) = -grad(:,ic) + beta(ic)*dir_old(:,ic)
      !ENDDO
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
      !WRITE(*,'(1x,A,F18.10)') 'alpha(ic) = ', alpha(ic)
      evecs(:,ic) = evecs(:,ic) + alpha(ic)*dir(:,ic)
    ENDDO
    !
    CALL ortho_gram_schmidt( evecs, N**3, N**3, Ncol )
    CALL get_Etot_v2( Ncol, evecs, Ekin, Epot, Etot )
    !
    deltaE = Etot_old - Etot
    IF( deltaE < 0.d0 ) THEN
      WRITE(*,*) 'Negative deltaE, restarting CG:', deltaE
      iter = 0
      CYCLE cgstep
    ENDIF
    WRITE(*,'(1x,A,I8,2F18.10)') 'iter, Etot, deltaE ', iter, Etot, Etot_old-Etot
    !WRITE(*,'(1x,A,F18.10)') 'Epot: ', Epot
    !WRITE(*,*) 'Ekin(ic)'
    !DO ic = 1, Ncol
    !  WRITE(*,'(1x,I8,F18.10)') ic, Ekin(ic)
    !ENDDO
    !
    IF( abs(Etot - Etot_old) < 1.d-6 ) THEN
      WRITE(*,'(1x,A,I8)') '***CG converged in iter', iterTrue
      EXIT
    ENDIF
    !
    Etot_old = Etot
    grad_old(:,:) = grad(:,:)
    dir_old(:,:) = dir(:,:)
    !
    IF( iter >= Niter ) THEN
      EXIT cgstep
    ENDIF
  ENDDO cgstep

  WRITE(111) evecs

  DEALLOCATE( grad, grad_old, grad_t, dir, dir_old )
END SUBROUTINE

