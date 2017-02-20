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
  !
  REAL(8) :: trace

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

  DO iter = 1, Niter
    !
    ! Evaluate gradient at current trial vectors
    CALL calc_grad( Nstates, v, grad )
    WRITE(*,*) 'sum(grad) = ', sum(grad)
    !STOP 'Stop at 54'
    !
    ! Calculate norm of the gradient, not really used in the minimization
    ! procedure.
    norm_grad = 0.d0
    DO ist = 1, Nstates
      norm_grad = norm_grad + sqrt( sum( grad(:,ist)**2 ) )
    ENDDO
    norm_grad = norm_grad/Nstates
    !
    ! set search direction
    !
    IF( iter == 1 ) THEN
      dir(:,:) = -grad(:,:)
    ELSE
      ! Fletcher-Reeves
      beta = trace( Nstates, matmul( transpose(grad), grad ) ) / &
             trace( Nstates, matmul( transpose(grad_old), grad_old ) )
    ENDIF
    dir(:,:) = -grad(:,:) + beta*dir_old(:,:)
    !
    ! Evaluate gradient at trial step
    tv(:,:) = v(:,:) + alpha_t*dir(:,:)
    tv(:,:) = tv(:,:)/sqrt(dVol)
    CALL calc_grad( Nstates, tv, grad_t )
    !
    ! Compute estimate of best step and update current trial vectors
    denum = trace( Nstates, matmul( transpose(grad - grad_t), dir ) )
    IF( abs(denum) >= 1e-9 ) THEN  ! FIXME: use abs ?
      alpha = abs( alpha_t * trace( Nstates, matmul(transpose(grad),dir) )/denum )
    ELSE 
      alpha = 0.d0
    ENDIF
    WRITE(*,*) 'iter, alpha_t, alpha, beta', iter, alpha_t, alpha, beta

    v(:,:) = v(:,:) + alpha*dir(:,:)
    CALL ortho_gram_schmidt( v, Npoints, Npoints, Nstates )
    v(:,:) = v(:,:)/sqrt(dVol)
    CALL calc_Energies( v, Ekin, Epot, Etot )
    !
    WRITE(*,'(1x,A,I5,2F18.10)') 'iter, conv, ||grad||: ', iter, abs(Etot-Etot_old), norm_grad
    WRITE(*,'(1x,A,3F18.10)') 'Ekin, Epot, Etot: ', Ekin, Epot, Etot
    !
    !IF( abs(Etot - Etot_old) < 1.d-7 ) THEN
    !  WRITE(*,*) 'sch_solve_Emin_cg converged in iter', iter
    !  EXIT
    !ENDIF
    !
    Etot_old = Etot
    grad_old(:,:) = grad(:,:)
    dir_old(:,:) = dir(:,:)
  ENDDO

  WRITE(111) v

  DEALLOCATE( v, grad, grad_old, grad_t, dir, dir_old, tv )
END SUBROUTINE

FUNCTION trace( dimA, A ) RESULT( trA )
  IMPLICIT NONE 
  INTEGER :: dimA
  REAL(8) :: A(dimA,dimA)
  REAL(8) :: trA
  INTEGER :: i

  trA = 0.d0
  DO i = 1, dimA
    trA = trA + A(i,i)
  ENDDO
END FUNCTION
