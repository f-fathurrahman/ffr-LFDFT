!----------------------------------------------------
SUBROUTINE sch_solve_Emin_cg( alpha_t, Niter, restart )
!----------------------------------------------------
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_states, ONLY : Nstates
  IMPLICIT NONE
  !
  INTEGER :: Niter
  REAL(8) :: alpha_t  ! step size
  LOGICAL :: restart
  REAL(8), ALLOCATABLE :: v(:,:), grad(:,:), grad_old(:,:), grad_t(:,:)
  REAL(8), ALLOCATABLE :: dir(:,:), dir_old(:,:), tv(:,:)
  REAL(8) :: Etot, Ekin, Epot, Etot_old, norm_grad
  REAL(8) :: alpha, beta, denum
  !
  INTEGER :: ist, iter, ip

  ALLOCATE( v(Npoints,Nstates) )
  ALLOCATE( grad(Npoints,Nstates) )
  ALLOCATE( grad_old(Npoints,Nstates) )
  ALLOCATE( grad_t(Npoints,Nstates) )
  ALLOCATE( dir(Npoints,Nstates) )
  ALLOCATE( dir_old(Npoints,Nstates) )

  ALLOCATE( tv(Npoints,Nstates) )

  ! initialize v
  IF( .NOT. restart ) THEN
    DO ist = 1, Nstates
      DO ip = 1, Npoints
        CALL random_number( v(ip,ist) )
      ENDDO
    ENDDO
    CALL ortho_gram_schmidt( v, Npoints, Npoints, Nstates )
    v(:,:) = v(:,:)/sqrt(dVol)
  ELSE
    READ(112) v
  ENDIF

  CALL calc_Energies( v, Ekin, Epot, Etot )
  WRITE(*,*) 'Initial Etot = ', Etot

  Etot_old = Etot

  alpha = 0.d0
  beta  = 0.d0

  grad(:,:)    = 0.d0
  grad_t(:,:)  = 0.d0
  dir(:,:)     = 0.d0
  dir_old(:,:) = 0.d0

  DO iter = 1, Niter
    !
    ! Evaluate gradient at current trial vectors
    CALL calc_grad( Nstates, v, grad )
    !
    ! Calculate norm of the gradient, not really used in the minimization
    ! procedure.
    !norm_grad = 0.d0
    !DO ist = 1, Nstates
    !  norm_grad = norm_grad + sqrt( sum( grad(:,ist)**2 ) )
    !ENDDO
    !norm_grad = norm_grad/Nstates
    !
    ! set search direction
    !
    IF( iter /= 1 ) THEN
      ! Fletcher-Reeves
      beta = sum( grad * grad ) / sum( grad_old * grad_old )
    ENDIF
    dir(:,:) = -grad(:,:) + beta*dir_old(:,:)
    !
    ! Evaluate gradient at trial step
    tv(:,:) = v(:,:) + alpha_t*dir(:,:)
    CALL ortho_gram_schmidt( tv, Npoints, Npoints, Nstates )
    tv(:,:) = tv(:,:)/sqrt(dVol)
    CALL calc_grad( Nstates, tv, grad_t )
    !
    ! Compute estimate of best step and update current trial vectors
    !denum = trace( Nstates, matmul( transpose(grad - grad_t), dir ) )
    denum = sum( (grad-grad_t) * dir )
    IF( denum /= 0.d0 ) THEN  ! FIXME: use abs ?
      alpha = abs( alpha_t * sum( grad * dir )/denum )
    ELSE 
      alpha = 0.d0
    ENDIF
    !WRITE(*,*) 'iter, alpha_t, alpha, beta', iter, alpha_t, alpha, beta

    v(:,:) = v(:,:) + alpha*dir(:,:)
    CALL ortho_gram_schmidt( v, Npoints, Npoints, Nstates )
    v(:,:) = v(:,:)/sqrt(dVol)
    CALL calc_Energies( v, Ekin, Epot, Etot )
    !
    !WRITE(*,'(/,1x,A,I5,2F18.10)') 'iter, conv, ||grad||: ', iter, abs(Etot-Etot_old), norm_grad
    !WRITE(*,'(1x,A,3F18.10)') 'Ekin, Epot, Etot: ', Ekin, Epot, Etot
    WRITE(*,'(1x,I5,F18.10,E18.10)') iter, Etot, abs(Etot-Etot_old)
    !
    IF( abs(Etot - Etot_old) < 1.d-7 ) THEN
      WRITE(*,*) 'sch_solve_Emin_cg converged in iter', iter
      EXIT
    ENDIF
    !
    Etot_old = Etot
    grad_old(:,:) = grad(:,:)
    dir_old(:,:) = dir(:,:)
  ENDDO

  WRITE(111) v

  DEALLOCATE( v, grad, grad_old, grad_t, dir, dir_old, tv )
END SUBROUTINE

