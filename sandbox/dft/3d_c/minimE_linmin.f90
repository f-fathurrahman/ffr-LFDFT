!-------------------------------------------------------
SUBROUTINE minimE_linmin( Ncol, alpha_t, Niter, restart )
!-------------------------------------------------------
  USE m_globals, ONLY : N
  IMPLICIT NONE
  !
  INTEGER :: Ncol, Niter
  REAL(8) :: alpha_t  ! step size
  LOGICAL :: restart
  ! TODO: make these variables allocatable instead of automatic in order
  !       to avoid stack error
  REAL(8), ALLOCATABLE :: v(:,:), grad(:,:), dir_v(:,:), grad_t(:,:)
  REAL(8) :: Etot, Ekin, Epot, Etot_old, norm_grad
  REAL(8) :: alpha(Ncol)
  !
  INTEGER :: ic, iter
  !
  REAL(8) :: ddot

  ALLOCATE( v(N**3,Ncol), grad(N**3,Ncol), dir_v(N**3,Ncol) )
  ALLOCATE( grad_t(N**3,Ncol) )

  ! initialize v
  IF( .NOT. restart ) THEN
    DO ic = 1, Ncol
      CALL r8_rand_vec( N**3, v(:,ic) )
    ENDDO
    CALL ortho_gram_schmidt( v, N**3, N**3, Ncol )
  ELSE
    READ(112) v
    ! No need to orthonormalize, v should already be ortonormal
  ENDIF

  !FIXME: Restart is not working properly yet. Without proper treatment
  !       it may give false minimum.


  CALL get_Etot( Ncol, v, Ekin, Epot, Etot)
  WRITE(*,*) 'Etot = ', Etot
  Etot_old = Etot

  DO iter = 1, Niter
    !
    ! Evaluate gradient at current trial vectors
    !
    CALL get_grad( Ncol, v, grad, .FALSE. )
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
    dir_v(:,:) = -grad(:,:)
    !
    ! Evaluate gradient at trial step
    !
    CALL get_grad( Ncol, v + alpha_t*dir_v, grad_t, .TRUE. )
    !
    ! Compute estimate of best step
    !
    DO ic = 1, Ncol
      alpha(ic) = alpha_t*ddot( N**3, grad(:,ic),1, dir_v(:,ic),1 )/&
            ddot( N**3, grad(:,ic)-grad_t(:,ic),1, dir_v(:,ic), 1 )
    ENDDO
    ! FIXME: Fuse this loop with the previous loop?
    DO ic = 1, Ncol
      v(:,ic) = v(:,ic) + alpha(ic)*dir_v(:,ic)
    ENDDO
    !
    CALL ortho_gram_schmidt( v, N**3, N**3, Ncol )
    CALL get_Etot( Ncol, v, Ekin, Epot, Etot )
    !
    WRITE(*,'(1x,A,I5,2F18.10)') 'iter, conv, ||grad||: ', iter, abs(Etot-Etot_old), norm_grad
    WRITE(*,'(1x,A,3F18.10)') 'Ekin, Epot, Etot: ', Ekin, Epot, Etot
    !
    IF( abs(Etot - Etot_old) < 1.d-7 ) THEN
      WRITE(*,*) 'linmin converged in iter', iter
      EXIT
    ENDIF
    !
    Etot_old = Etot
  ENDDO

  WRITE(111) v

  DEALLOCATE( v, grad, dir_v, grad_t )

END SUBROUTINE
