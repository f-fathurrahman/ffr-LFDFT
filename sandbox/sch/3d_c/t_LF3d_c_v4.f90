! efefer, 15 January 2016
!
! Solution of Schrodinger equation

! Using iterative (partial) diagonalization

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
  !
  REAL(8), ALLOCATABLE :: Kprec(:)
END MODULE

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

  r0(:) = (B-A)/2.d0 + 1d-8 ! position of H atom

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

  r0(:) = (B-A)/2.d0 + 1.d-8

  ! N is assumed to be the same as LF$N
  DO ip=1,N**3
    V(ip) = -Z/norm2( LF%lingrid(:,ip) - r0(:) )
  ENDDO
END SUBROUTINE


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


SUBROUTINE h_psi( npwx, npw, nvec, psi, hpsi )
  IMPLICIT NONE
  INTEGER :: npwx, npw, nvec
  REAL(8) :: psi(npwx,nvec)
  REAL(8) :: hpsi(npwx,nvec)
  !
  INTEGER :: iv

  DO iv=1,nvec
    CALL apply_Ham( psi(:,iv), hpsi(:,iv) )
  ENDDO

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


SUBROUTINE diag_lanczos( Ncol, evals, evecs )
  USE m_globals, ONLY : N
  IMPLICIT NONE
  !
  INTEGER :: Ncol
  REAL(8) :: evals(Ncol)
  REAL(8) :: evecs(N**3,Ncol)
  !
  REAL(8), ALLOCATABLE :: Vm(:,:), W(:)  ! Lanczos vectors
  INTEGER :: j
  REAL(8) :: r1
  REAL(8) :: alpha(Ncol), beta(Ncol)  ! alpha and beta of the triangular matrix
  REAL(8) :: evecsT(Ncol,Ncol)  ! the eigenvector of triangular matrix
  REAL(8), ALLOCATABLE :: work(:)
  INTEGER :: info
  !
  REAL(8) :: ddot

  ALLOCATE( Vm(N**3,Ncol), W(N**3) )
  ALLOCATE( work(2*Ncol-1) )

  CALL r8_rand_vec( N**3,Vm(:,1) )
  r1 = sqrt( ddot( N**3, Vm(:,1),1, Vm(:,1),1 ) )
  Vm(:,1) = Vm(:,1)/r1

  r1 = sqrt( ddot(N**3,Vm(:,1),1,Vm(:,1),1) )

  WRITE(*,*) 'r1 = ', r1
  
  DO j = 1, Ncol-1
    CALL apply_Ham( Vm(:,j), W(:) )
    alpha(j) = ddot( N**3, W(:),1, Vm(:,j),1 )
    !WRITE(*,*) 'j, alpha = ', j, alpha(j)
    IF(j>1) THEN
      W(:) = W(:) - alpha(j)*Vm(:,j) - beta(j-1)*Vm(:,j-1)
    ELSE
      W(:) = W(:) - alpha(j)*Vm(:,j) 
    ENDIF
    beta(j) = sqrt( ddot( N**3, W(:),1, W(:),1 ) )
    Vm(:,j+1) = W(:)/beta(j)
    r1 = sqrt( ddot(N**3, Vm(:,j+1),1, Vm(:,j+1),1 ) )
    WRITE(*,*) 'j, r1 = ', r1
  ENDDO

  CALL apply_Ham( Vm(:,Ncol), W(:) )
  alpha(Ncol) = ddot( N**3, W(:),1, Vm(:,ncol),1 )
  !WRITE(*,*) 'Ncol, alpha = ', Ncol,alpha(Ncol)

  IF( Ncol > 1 ) THEN
    CALL dstev('V',Ncol, alpha, beta, evecsT,Ncol, work, info )
    WRITE(*,*) 'info = ', info
    evals(:) = alpha(:)
    WRITE(*,*) 'Pass here'
  ELSE
    evals(1) = alpha(1)
  ENDIF

  WRITE(*,*) evals(:)

  DEALLOCATE( Vm, W )
  DEALLOCATE( work )


END SUBROUTINE


!------------------------------------------------------------------------------
PROGRAM t_LF3d_c
!------------------------------------------------------------------------------
  USE m_constants, ONLY: PI
  USE m_LF3d
  USE m_globals
  IMPLICIT NONE
  !
  REAL(8), ALLOCATABLE :: evals(:), evecs(:,:)
  INTEGER :: Nstates
  INTEGER, ALLOCATABLE :: btype(:)
  INTEGER :: dav_iter, gstart, notcnv
  INTEGER :: lrot
  REAL(8) :: ethr
  INTEGER :: ist
  
  N = 31
  A = 0.d0
  B = 6.d0

  CALL init_LF3d_c(LF, (/N,N,N/), (/A,A,A/), (/B,B,B/) )

  ! Set up potential
  ALLOCATE( Vpot(N**3) )
  !CALL init_pot_harmonic( 2.d0, Vpot )
  CALL init_pspot_H( Vpot )
  !CALL init_pot_coulomb(1.d0, Vpot)

  Nstates = 1
  ALLOCATE( evecs(N**3,Nstates), evals(Nstates) )

  !CALL diag_lanczos( Nstates, evals, evecs )

  gstart = 2
  ethr = 1.d-5
  ALLOCATE( btype(Nstates) )
  btype(:) = 1 ! all bands are occupied
  lrot = .FALSE.
 
  WRITE(*,*) 'Pass here ...'
  DO ist = 1, Nstates
    CALL r8_rand_vec( N**3, evecs(:,ist) )
  ENDDO
  CALL ortho_gram_schmidt( evecs, N**3, N**3, Nstates )

  WRITE(*,*) 'Calling regterg ...'
  CALL regterg( N**3, N**3, Nstates, 4*Nstates, evecs, ethr, &
                    gstart, evals, btype, notcnv, lrot, dav_iter )
  WRITE(*,*) 'dav_iter = ', dav_iter
  WRITE(*,*) evals

  DEALLOCATE( Vpot )
  DEALLOCATE( evecs, evals )

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


