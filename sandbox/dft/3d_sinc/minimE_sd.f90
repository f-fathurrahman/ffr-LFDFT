!-------------------------------------------------
SUBROUTINE minimE_sd( Ncol, alpha, Niter, restart )
!-------------------------------------------------
  USE m_globals, ONLY : N
  IMPLICIT NONE
  !
  INTEGER :: Ncol, Niter
  REAL(8) :: alpha  ! step size
  LOGICAL :: restart
  !
  REAL(8), ALLOCATABLE :: v(:,:), grad(:,:)
  REAL(8) :: Etot, Ekin, Epot, Etot_old, norm_grad
  !
  INTEGER :: ic, iter

  ALLOCATE( v(N**3,Ncol), grad(N**3,Ncol) )

  ! initialize v
  IF( .NOT. restart ) THEN
    DO ic = 1, Ncol
      CALL r8_rand_vec( N**3, v(:,ic) )
    ENDDO
  ELSE
    READ(112) v
  ENDIF

  CALL ortho_gram_schmidt( v, N**3, N**3, Ncol )

  CALL get_Etot( Ncol, v, Ekin, Epot, Etot)
  WRITE(*,*) 'Initial Etot = ', Etot
  Etot_old = Etot

  DO iter = 1, Niter
    CALL get_grad( Ncol, v, grad, .FALSE. )
    !
    norm_grad = 0.d0
    DO ic = 1, Ncol
      norm_grad = norm_grad + norm2( grad(:,ic) )
    ENDDO
    norm_grad = norm_grad/Ncol
    !
    v = v - alpha*grad
    !
    CALL ortho_gram_schmidt( v, N**3, N**3, Ncol )
    CALL get_Etot( Ncol, v, Ekin, Epot, Etot )
    !
    WRITE(*,'(1x,A,I5,2F18.10)') 'iter, conv, ||grad||: ', iter, abs(Etot-Etot_old), norm_grad
    WRITE(*,'(1x,A,3F18.10)') 'Ekin, Epot, Etot: ', Ekin, Epot, Etot
    ! Convergence criteria in hartree
    IF( abs(Etot - Etot_old) < 1.d-7 ) THEN
      WRITE(*,*) 'SD converged in iter', iter
      EXIT
    ENDIF
    Etot_old = Etot
  ENDDO

  WRITE(111) v

  DEALLOCATE( v, grad )

END SUBROUTINE
