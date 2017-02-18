! efefer, 16 January 2016
!
! Solution of Schrodinger equation using conjugate gradient method
! with preconditioner as proposed in PhysRevB.51.14057



! This module contains several global variables and parameter used in this
! program
MODULE m_globals
  USE m_LF3d
  IMPLICIT NONE
  ! These parameters are similar for x, y, and z directions
  INTEGER :: N
  REAL(8) :: A, B
  REAL(8) :: h
  !
  TYPE(LF3d_t) :: LF
  !
  REAL(8), ALLOCATABLE :: Vpot(:)
  REAL(8), ALLOCATABLE :: eval(:)
  !
  REAL(8), ALLOCATABLE :: Kprec(:)
END MODULE


!-------------------------------------
SUBROUTINE init_pot_harmonic(omega, V)
!-------------------------------------
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : N, LF
  IMPLICIT NONE
  ! Arguments
  REAL(8) :: omega
  REAL(8) :: V(N**3)
  ! local
  INTEGER :: ip
  REAL(8) :: r0(3)

  r0(:) = 0.d0

  ! N is assumed to be the same as LF$N
  DO ip=1,N**3
    V(ip) = 0.5d0*omega**2 * norm2( LF%lingrid(:,ip) - r0(:) )**2
  ENDDO
END SUBROUTINE

!-------------------------------------
SUBROUTINE init_pspot_H( V )
!-------------------------------------
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : N, LF, A, B
  IMPLICIT NONE
  ! Arguments
  REAL(8) :: V(N**3)
  ! local
  INTEGER :: ip
  REAL(8) :: r0(3)
  !
  REAL(8) :: rc1, rc2, aa, bb, r

  r0(:) = 0.d0 ! position of H atom

  rc1 = 0.25d0
  rc2 = 0.284d0
  aa = -1.9287d0
  bb = 0.3374d0

  ! FIXME Only pure radial potential
  DO ip=1,N**3
    r = norm2( LF%lingrid(:,ip)-r0(:) )
    !WRITE(*,*)
    V(ip) = -1.d0/r * erf( r/rc1 ) + (aa + bb*r**2)*exp(-(r/rc2)**2)
    !WRITE(113,*) r, V(ip)
  ENDDO
END SUBROUTINE


!--------------------------------
SUBROUTINE init_pot_coulomb(Z, V)
!--------------------------------
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : N, LF, A, B
  IMPLICIT NONE
  ! Arguments
  REAL(8) :: Z
  REAL(8) :: V(N**3)
  ! local
  INTEGER :: ip
  REAL(8) :: r0(3)

  r0(:) = 0.d0

  ! N is assumed to be the same as LF$N
  DO ip=1,N**3
    V(ip) = -Z/norm2( LF%lingrid(:,ip) - r0(:) )
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


! 1-column version
SUBROUTINE apply_Ham_v2(v, Hv, Kin, lambda)
  USE m_globals, ONLY : N, Vpot
  IMPLICIT NONE
  !
  REAL(8) :: v(N**3)
  REAL(8) :: Hv(N**3)
  REAL(8) :: Kin, lambda
  !
  REAL(8) :: ddot

  CALL apply_laplacian(v, Hv)
  !
  Kin = -0.5d0*ddot( N**3, v,1, Hv,1 )
  !
  Hv(:) = -0.5d0*Hv(:) + Vpot(:)*v(:)
  !
  lambda = ddot( N**3, v,1, Hv,1 )
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



! multicolumn version
! v are asssumed to be orthogonal
SUBROUTINE get_grad( Ncol, v, grad, do_ortho)
  USE m_globals, ONLY : N
  IMPLICIT NONE
  !
  INTEGER :: Ncol
  REAL(8) :: v(N**3,Ncol)
  REAL(8) :: grad(N**3,Ncol)
  LOGICAL, OPTIONAL :: do_ortho  ! FIXME The keyword OPTIONAL is not working (?)
  !
  REAL(8) :: Hv(N**3)
  INTEGER :: ic, icc
  !
  REAL(8) :: ddot

  IF( .NOT. present( do_ortho ) ) do_ortho = .FALSE.

  IF( do_ortho ) THEN
    CALL ortho_gram_schmidt( v, N**3, N**3, ncol )
  ENDIF

  !WRITE(*,*) 'Calling get_grad'
  !
  DO ic = 1, Ncol
    CALL apply_Ham( v(:,ic), Hv(:) )
    grad(:,ic) = Hv(:)
    DO icc = 1, Ncol
      grad(:,ic) = grad(:,ic) - ddot( N**3, v(:,icc),1, Hv(:),1 )*v(:,icc)
    ENDDO
  ENDDO
END SUBROUTINE



! multicolumn version
! v are asssumed to be orthogonal
SUBROUTINE get_grad_v2( Ncol, v, grad, do_ortho, Kin, lambda)
  USE m_globals, ONLY : N
  IMPLICIT NONE
  !
  INTEGER :: Ncol
  REAL(8) :: v(N**3,Ncol)
  REAL(8) :: grad(N**3,Ncol)
  REAL(8) :: lambda(Ncol)
  REAL(8) :: Kin(Ncol)
  LOGICAL, OPTIONAL :: do_ortho  ! FIXME The keyword OPTIONAL is not working (?)
  !
  REAL(8) :: Hv(N**3)
  INTEGER :: ic, icc
  !
  REAL(8) :: ddot

  IF( .NOT. present( do_ortho ) ) do_ortho = .FALSE.

  IF( do_ortho ) THEN
    CALL ortho_gram_schmidt( v, N**3, N**3, ncol )
  ENDIF

  !WRITE(*,*) 'Calling get_grad'
  !
  DO ic = 1, Ncol
    CALL apply_Ham_v2( v(:,ic), Hv(:), Kin(ic), lambda(ic) )
    grad(:,ic) = Hv(:)
    DO icc = 1, Ncol
      grad(:,ic) = grad(:,ic) - ddot( N**3, v(:,icc),1, Hv(:),1 )*v(:,icc)
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


! psi are assumed to be orthogonalized
!-----------------------------------------------
SUBROUTINE get_Etot(Ncol, psi, Ekin, Epot, Etot)
!-----------------------------------------------
  USE m_globals, ONLY : N, Vpot, LF
  IMPLICIT NONE
  !
  INTEGER :: Ncol
  REAL(8) :: psi(N**3,Ncol)
  REAL(8) :: nabla2_psi(N**3)
  REAL(8) :: Etot, Ekin, Epot
  !
  INTEGER :: ic
  REAL(8) :: deltaV
  !
  REAL(8) :: ddot

  deltaV = LF%LFx%h * LF%LFy%h * LF%LFz%h
  !
  Etot = 0.d0
  Ekin = 0.d0
  Epot = 0.d0
  ! assume all occupancies are 1.d0
  DO ic=1,Ncol
    CALL apply_laplacian( psi(:,ic), nabla2_psi(:) )
    Ekin = Ekin + -0.5d0*ddot( N**3, psi(:,ic),1, nabla2_psi(:),1 )
    !
    Epot = Epot + sum( Vpot(:)*psi(:,ic)**2 ) ! FIXME: We don't need to multiply to deltaV here
  ENDDO
  Etot = Ekin + Epot
END SUBROUTINE

! CG algorithm, no preconditioning
!------------------------------------------------------------
SUBROUTINE solve_cg( Ncol, alpha_t, Niter, restart )
!------------------------------------------------------------
  USE m_globals, ONLY : N
  IMPLICIT NONE
  !
  INTEGER :: Ncol, Niter
  REAL(8) :: alpha_t  ! step size
  LOGICAL :: restart
  ! TODO: make these variables allocatable instead of automatic in order
  !       to avoid stack error
  REAL(8), ALLOCATABLE :: v(:,:), grad(:,:), grad_old(:,:), grad_t(:,:)
  REAL(8), ALLOCATABLE :: dir(:,:), dir_old(:,:)
  REAL(8) :: Etot, Ekin, Epot, Etot_old, norm_grad
  REAL(8) :: alpha(Ncol), beta(Ncol)
  !
  INTEGER :: ic, iter
  !
  REAL(8) :: ddot

  ALLOCATE( v(N**3,Ncol), grad(N**3,Ncol), grad_old(N**3,Ncol), grad_t(N**3,Ncol) )
  ALLOCATE( dir(N**3,Ncol), dir_old(N**3,Ncol) )

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

  CALL get_Etot( Ncol, v, Ekin, Epot, Etot)
  WRITE(*,*) 'Initial Etot = ', Etot
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
    IF( iter == 1 ) THEN
      dir(:,:) = -grad(:,:)
    ELSE
      ! Fletcher-Reeves
      DO ic=1,Ncol
        beta(ic) = ddot( N**3, grad(:,ic),1, grad(:,ic),1 ) / &
                   ddot( N**3, grad_old(:,ic),1, grad_old(:,ic),1 )
        dir(:,ic) = -grad(:,ic) + beta(ic)*dir_old(:,ic)
      ENDDO
      ! Polak-Ribiere, not as efficient as Fletcher-Reeves?
      !DO ic=1,Ncol
      !  beta(ic) = ddot( N**3, grad(:,ic)-grad_old(:,ic),1, grad(:,ic),1 ) / &
      !             ddot( N**3, grad_old(:,ic),1, grad_old(:,ic),1 )
      !  dir(:,ic) = -grad(:,ic) + beta(ic)*dir_old(:,ic)
      !ENDDO
    ENDIF
    !
    ! Evaluate gradient at trial step
    !
    CALL get_grad( Ncol, v + alpha_t*dir, grad_t, .TRUE. )
    !
    ! Compute estimate of best step and update current trial vectors
    !
    DO ic = 1, Ncol
      alpha(ic) = alpha_t*ddot( N**3, grad(:,ic),1, dir(:,ic),1 )/&
            ddot( N**3, grad(:,ic)-grad_t(:,ic),1, dir(:,ic), 1 )
      v(:,ic) = v(:,ic) + alpha(ic)*dir(:,ic)
    ENDDO
    !
    CALL ortho_gram_schmidt( v, N**3, N**3, Ncol )
    CALL get_Etot( Ncol, v, Ekin, Epot, Etot )
    !
    WRITE(*,'(1x,A,I5,2F18.10)') 'iter, conv, ||grad||: ', iter, abs(Etot-Etot_old), norm_grad
    WRITE(*,'(1x,A,3F18.10)') 'Ekin, Epot, Etot: ', Ekin, Epot, Etot
    !
    IF( abs(Etot - Etot_old) < 1.d-7 ) THEN
      WRITE(*,*) 'cg converged in iter', iter
      EXIT
    ENDIF
    !
    Etot_old = Etot
    grad_old(:,:) = grad(:,:)
    dir_old(:,:) = dir(:,:)
  ENDDO

  WRITE(111) v

  DEALLOCATE( v, grad, grad_old, grad_t, dir, dir_old )
END SUBROUTINE


! Seitsonen-Puska-Nieminen (PRB-51-14057)
!---------------------------------------------------------------
SUBROUTINE solve_pcg_SPN( Ncol, alpha_t, Niter, restart, Aprec )
!---------------------------------------------------------------
  USE m_globals, ONLY : N, Vpot
  IMPLICIT NONE
  !
  INTEGER :: Ncol, Niter
  REAL(8) :: alpha_t  ! step size
  LOGICAL :: restart
  ! TODO: make these variables allocatable instead of automatic in order
  !       to avoid stack error
  REAL(8), ALLOCATABLE :: v(:,:), grad(:,:), grad_old(:,:), grad_t(:,:)
  REAL(8), ALLOCATABLE :: dir(:,:), dir_old(:,:)
  REAL(8) :: Etot, Ekin, Epot, Etot_old, norm_grad
  REAL(8) :: alpha(Ncol), beta(Ncol), lambda(Ncol), Kin(Ncol)
  REAL(8), ALLOCATABLE :: Kprec(:), x(:)
  REAL(8) :: Aprec
  !
  INTEGER :: ic, iter
  !
  REAL(8) :: ddot

  ALLOCATE( v(N**3,Ncol), grad(N**3,Ncol), grad_old(N**3,Ncol), grad_t(N**3,Ncol) )
  ALLOCATE( dir(N**3,Ncol), dir_old(N**3,Ncol) )
  ALLOCATE( Kprec(N**3), x(N**3) )

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

  CALL get_Etot( Ncol, v, Ekin, Epot, Etot)
  WRITE(*,*) 'Initial Etot = ', Etot
  Etot_old = Etot

  DO iter = 1, Niter
    !
    ! Evaluate gradient at current trial vectors
    !
    CALL get_grad_v2( Ncol, v, grad, .FALSE., Kin, lambda )
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
        x(:) = Aprec*abs( lambda(ic) - Vpot(:) )/Kin(ic)
        Kprec(:) = ( 27d0 + 18d0*x(:) + 12d0*x(:)**2 + 8d0*x(:)**3 ) / &
                   ( 27d0 + 18d0*x(:) + 12d0*x(:)**2 + 8d0*x(:)**3 + 16d0*x(:)**4 )
        dir(:,ic) = -Kprec(:)*grad(:,ic)
      ENDDO
    ELSE
      ! Fletcher-Reeves
      DO ic=1,Ncol
        x(:) = Aprec*abs( lambda(ic) - Vpot(:) )/Kin(ic)
        Kprec(:) = ( 27d0 + 18d0*x(:) + 12d0*x(:)**2 + 8d0*x(:)**3 ) / &
                   ( 27d0 + 18d0*x(:) + 12d0*x(:)**2 + 8d0*x(:)**3 + 16d0*x(:)**4 )
        beta(ic) = ddot( N**3, grad(:,ic),1, Kprec(:)*grad(:,ic),1 ) / &
                   ddot( N**3, grad_old(:,ic),1, Kprec(:)*grad_old(:,ic),1 )
        dir(:,ic) = -Kprec(:)*grad(:,ic) + beta(ic)*dir_old(:,ic)
      ENDDO
      ! Polak-Ribiere, not as efficient as Fletcher-Reeves?
      !DO ic=1,Ncol
      !  beta(ic) = ddot( N**3, grad(:,ic)-grad_old(:,ic),1, grad(:,ic),1 ) / &
      !             ddot( N**3, grad_old(:,ic),1, grad_old(:,ic),1 )
      !  dir(:,ic) = -grad(:,ic) + beta(ic)*dir_old(:,ic)
      !ENDDO
    ENDIF
    !
    ! Evaluate gradient at trial step
    !
    CALL get_grad( Ncol, v + alpha_t*dir, grad_t, .TRUE. )
    !
    ! Compute estimate of best step and update current trial vectors
    !
    DO ic = 1, Ncol
      alpha(ic) = alpha_t*ddot( N**3, grad(:,ic),1, dir(:,ic),1 )/&
            ddot( N**3, grad(:,ic)-grad_t(:,ic),1, dir(:,ic), 1 )
      v(:,ic) = v(:,ic) + alpha(ic)*dir(:,ic)
    ENDDO
    !
    CALL ortho_gram_schmidt( v, N**3, N**3, Ncol )
    CALL get_Etot( Ncol, v, Ekin, Epot, Etot )
    !
    WRITE(*,'(1x,A,I5,2F18.10)') 'iter, conv, ||grad||: ', iter, abs(Etot-Etot_old), norm_grad
    WRITE(*,'(1x,A,3F18.10)') 'Ekin, Epot, Etot: ', Ekin, Epot, Etot
    !
    IF( abs(Etot - Etot_old) < 1.d-7 ) THEN
      WRITE(*,*) 'pcg converged in iter', iter
      EXIT
    ENDIF
    !
    Etot_old = Etot
    grad_old(:,:) = grad(:,:)
    dir_old(:,:) = dir(:,:)
  ENDDO

  WRITE(111) v

  DEALLOCATE( v, grad, grad_old, grad_t, dir, dir_old, Kprec, x )
END SUBROUTINE


!------------------------------------------------------------------------------
PROGRAM t_LF3d_c
!------------------------------------------------------------------------------
  USE m_constants, ONLY: PI
  USE m_LF3d
  USE m_globals
  IMPLICIT NONE
  !
  INTEGER :: ip
  REAL(8), ALLOCATABLE :: vec(:)
  

  N = 50
  h = 0.25d0

  CALL init_LF3d_sinc(LF, (/N,N,N/), (/h,h,h/) )
  
  ! A and B are the same in x,y,z directions
  A = LF%LFx%A
  B = LF%LFx%B

  ! Set up potential
  ALLOCATE( Vpot(N**3) )
  CALL init_pspot_H( Vpot )
  !CALL init_pot_coulomb(1.d0, Vpot)

  CALL solve_pcg_SPN( 2, 3.d-5, 1000, .FALSE., 0.05 )
  !CALL solve_cg( 2, 3.d-5, 1000, .FALSE. )

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


