! efefer, 1 January 2016
!
! Solution of Schrodinger equation


! Steepest descent solution is VERY slow.
!
! linmin works quite well. Unnecessarily large N, however, will converge
! more slowly.

! This module contains several global variables and parameter used in this
! program
MODULE m_globals
  USE m_LF3d
  IMPLICIT NONE
  ! These parameters are similar for x, y, and z directions
  INTEGER :: N
  REAL(8) :: A, B
  !
  TYPE(LF3d_t) :: LF
  !
  REAL(8), ALLOCATABLE :: Vpot(:)
  REAL(8), ALLOCATABLE :: eval(:)
END MODULE


!-------------------------------------
SUBROUTINE init_pot_harmonic(omega, V)
!-------------------------------------
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : N, LF, A, B
  IMPLICIT NONE
  ! Arguments
  REAL(8) :: omega
  REAL(8) :: V(N**3)
  ! local
  INTEGER :: ip
  REAL(8) :: r0(3)

  r0(:) = (B-A)/2.d0

  ! N is assumed to be the same as LF$N
  DO ip=1,N**3
    V(ip) = 0.5d0*omega**2 * norm2( LF%lingrid(:,ip) - r0(:) )**2
  ENDDO
END SUBROUTINE


! 1-column version
SUBROUTINE apply_Ham(v, Hv)
  USE m_globals, ONLY : N, Vpot
  IMPLICIT NONE
  !
  REAL(8) :: v(N**3)
  REAL(8) :: Hv(N**3)

  CALL apply_laplacian(v, Hv)
  !
  Hv(:) = -0.5d0*Hv(:) + Vpot(:)*v(:)
END SUBROUTINE


! apply Laplacian to input vector v
SUBROUTINE apply_laplacian(v, nabla2_v)
  USE m_globals, ONLY : N, LF
  IMPLICIT NONE
  REAL(8) :: v(N**3)
  REAL(8) :: nabla2_v(N**3)
  !
  INTEGER :: i, j, k, ip, ii, jj, kk
  !
  ! N**3 should be the same as LF3d%N
  DO ip=1,N**3
    i = LF%lin2xyz(1,ip)
    j = LF%lin2xyz(2,ip)
    k = LF%lin2xyz(3,ip)
    !
    nabla2_v(ip) = 0.d0
    !
    DO ii=1,N
      nabla2_v(ip) = nabla2_v(ip) + LF%LFx%D2jl(ii,i)*v(LF%xyz2lin(ii,j,k))
    ENDDO
    !
    DO jj=1,N
      nabla2_v(ip) = nabla2_v(ip) + LF%LFy%D2jl(jj,j)*v(LF%xyz2lin(i,jj,k))
    ENDDO
    !
    DO kk=1,N
      nabla2_v(ip) = nabla2_v(ip) + LF%LFz%D2jl(kk,k)*v(LF%xyz2lin(i,j,kk))
    ENDDO
  ENDDO

END SUBROUTINE



! real(8) version
!------------------------------------------------
SUBROUTINE ortho_gram_schmidt(v, ldv, nrow, ncol)
!------------------------------------------------
  IMPLICIT NONE
  INTEGER :: ldv, nrow, ncol
  REAL(8) :: v(ldv,ncol)
  !
  INTEGER :: ii, jj
  REAL(8) :: zz, puv
  !
  REAL(8) :: ddot

  DO ii = 1, ncol
    zz = ddot( nrow, v(1:nrow,ii),1, v(1:nrow,ii),1 )
    v(1:nrow,ii) = v(1:nrow,ii)/sqrt( zz )
    DO jj = ii+1, ncol
      puv = prj( nrow, v(1:nrow,ii), v(1:nrow,jj) )
      v(1:nrow,jj) = v(1:nrow,jj) - puv*v(1:nrow,ii)
    ENDDO
  ENDDO

  CONTAINS

    ! compute prj = <v|u>/<u|u>
    FUNCTION prj(N,u,v)
      IMPLICIT NONE
      !
      REAL(8) :: prj
      INTEGER :: N
      REAL(8) :: u(N), v(N)
      !
      REAL(8) :: vu, uu
      REAL(8) :: ddot
      !
      ! FIXME: I got the vectors to be orthogonal when I reverse the arguments
      ! for zdotc
      vu = ddot( N, u,1, v,1 )
      uu = ddot( N, u,1, u,1 )
      prj = vu/uu
    END FUNCTION

END SUBROUTINE


!---------------------------
SUBROUTINE r8_rand_vec(N, v)
!---------------------------
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: v(N)
  !
  INTEGER :: i

  DO i=1,N
    CALL random_number( v(i) )
  ENDDO

END SUBROUTINE



! FIXME: Not yet finalized. In this subroutine, we are supposed to test
! the implementation of get_grad() against finite difference result.
!-------------------------
SUBROUTINE test_grad(Ncol)
!-------------------------
  USE m_globals, ONLY : N
  IMPLICIT NONE
  !
  INTEGER :: Ncol
  REAL(8) :: v(N**3,Ncol)
  !
  INTEGER :: ic
  REAL(8) :: Etot, Ekin, Epot
  !
  REAL(8) :: ddot

  WRITE(*,*)
  WRITE(*,*) 'Initializing random vectors:'
  DO ic=1,Ncol
    CALL r8_rand_vec(N**3, v(:,ic))
    WRITE(*,*) ic, ddot( N**3, v(:,ic),1, v(:,ic),1 )
  ENDDO

  CALL ortho_gram_schmidt( v, N**3, N**3, Ncol)
  WRITE(*,*)
  WRITE(*,*) 'Norm:'
  DO ic=1,Ncol
    WRITE(*,*) ic, ddot( N**3, v(:,ic),1, v(:,ic),1 )
  ENDDO
  WRITE(*,*)
  WRITE(*,*) 'Norm wrt to col #1:'
  DO ic=2,Ncol
    WRITE(*,*) ic, ddot( N**3, v(:,ic),1, v(:,1),1 )
  ENDDO

  CALL get_Etot( Ncol, v, Ekin, Epot, Etot )
  WRITE(*,*) 'Ekin = ', Ekin
  WRITE(*,*) 'Epot = ', Epot
  WRITE(*,*) 'Etot = ', Etot
 
END SUBROUTINE


!-------------------------------------------------
SUBROUTINE solve_SD( Ncol, alpha, Niter, restart )
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


!-------------------------------------------------------
SUBROUTINE solve_linmin( Ncol, alpha_t, Niter, restart )
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




!------------------------------------------------------------------------------
PROGRAM t_LF3d_c
!------------------------------------------------------------------------------
  USE m_constants, ONLY: PI
  USE m_LF3d
  USE m_globals
  IMPLICIT NONE
  !
  INTEGER :: ii
  ! Parameter for harmonic potential
  REAL(8), PARAMETER :: omega=2.d0
  !
  
  N = 20
  A = 0.d0
  B = 6.d0

  CALL init_LF3d_c(LF, (/N,N,N/), (/A,A,A/), (/B,B,B/) )
  
  ! Set up potential
  ALLOCATE( Vpot(N**3) )
  CALL init_pot_harmonic( omega, Vpot )


  !CALL test_grad(4)

  !CALL solve_SD( 4, 3.d-4, 3000, .FALSE. )
  CALL solve_linmin( 4, 3.d-5, 1000, .FALSE. )
  !CALL solve_pcg( 4, 3.d-5, 1000, .FALSE. )

  !ALLOCATE( eval(N**3) )
  !
  !WRITE(*,*)
  !WRITE(*,*) 'Eigenvalues:'
  !DO ii=1,4
  !  WRITE(*,'(1x,I5,F18.10)') ii, eval(ii)
  !ENDDO

  !DEALLOCATE( eval )
  DEALLOCATE( Vpot )

END PROGRAM



!----------------------------------
SUBROUTINE eig_dsyev(A,eigval,dimA)
!----------------------------------
  IMPLICIT NONE
  ! Arguments
  INTEGER :: dimA
  REAL(8) :: A(dimA,dimA)
  REAL(8) :: eigval(dimA)
  ! Workspace array for DSYEV
  INTEGER :: lwork
  REAL(8), ALLOCATABLE :: work(:)
  ! Local variables
  INTEGER :: info

  lwork = 3*dimA-1
  ALLOCATE(work(lwork))

  CALL dsyev('v', 'u', dimA, A, dimA, eigval, work, lwork, info)
  IF(info /= 0) THEN
    WRITE(*,*) 'Error on calling dsyev: info=',info
    STOP
  ENDIF

  ! Free memory
  DEALLOCATE(work)
END SUBROUTINE


