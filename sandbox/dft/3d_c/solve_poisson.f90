!--------------------------------------------
SUBROUTINE solve_poisson_cg_v2( N, rho, phi )
!--------------------------------------------
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: rho(N), phi(N)
  !
  REAL(8), ALLOCATABLE :: nabla2vec(:), r(:), z(:), p(:)
  REAL(8), ALLOCATABLE :: r_old(:), z_old(:), p_old(:)
  INTEGER :: ip, iter
  REAL(8) :: alpha_k, norm_res, beta_k
  LOGICAL :: conv
  !
  REAL(8) :: ddot

  ALLOCATE( nabla2vec(N), r(N), z(N), p(N) )
  ALLOCATE( r_old(N), z_old(N), p_old(N) )
  
  ! guess of phi
  DO ip=1,N
    phi(ip) = 0.d0
  ENDDO

  CALL apply_Laplacian( phi, nabla2vec )
  r(:) = rho(:) - nabla2vec(:)
  z(:) = r(:)
  p(:) = z(:)

  r_old(:) = r(:)
  z_old(:) = z(:)
  p_old(:) = p(:)

  conv = .FALSE.
  DO iter = 1, 1000
    CALL apply_Laplacian( p, nabla2vec )
    alpha_k = ddot( N, r,1, z,1 )/ddot( N, p,1, nabla2vec,1 )
    !
    phi(:) = phi(:) + alpha_k*p(:)
    !
    r(:) = r(:) - alpha_k*nabla2vec(:)
    norm_res = sqrt( ddot( N, r,1, r,1 ) )
    !
    !WRITE(*,*) 'iter, norm_res = ', iter, norm_res
    IF( norm_res < 1.d-9 ) THEN
      WRITE(*,*) 'Convergence achieved in solve_poisson_cg_v2: N, iter', N, iter
      conv = .TRUE.
      EXIT
    ENDIF
    !
    !CALL precond( N, r, z )
    z(:) = r(:)
    ! Fletcher-Reeves
    !beta_k = ddot( N, z,1, r,1 ) / ddot( N, z_old,1, r_old,1 )
    ! Polak-Ribiere
    beta_k = ddot( N, z,1, r-r_old,1 ) / ddot( N, z_old,1, r_old,1 )
    p(:) = z(:) + beta_k*p_old(:)
    !
    z_old(:) = z(:)
    r_old(:) = r(:)
    p_old(:) = p(:)

  ENDDO

  IF( .NOT. conv ) WRITE(*,*) 'No convergence in solve_poisson_pcg, N', N

  DEALLOCATE( nabla2vec, r, z, p )
  DEALLOCATE( r_old, z_old, p_old )
  
END SUBROUTINE



SUBROUTINE solve_poisson_cg( Nbasis, rho, phi )
  IMPLICIT NONE
  !
  INTEGER :: Nbasis
  REAL(8) :: rho(Nbasis), phi(Nbasis)
  !
  REAL(8), ALLOCATABLE :: r(:), p(:) ! residual
  REAL(8), ALLOCATABLE :: nabla2_phi(:)
  INTEGER :: ip, iter
  REAL(8) :: rsold, rsnew, alpha
  REAL(8) :: ddot

  ALLOCATE( r(Nbasis), p(Nbasis), nabla2_phi(Nbasis) )

  !
  DO ip=1,Nbasis
    phi(ip) = 0.d0
  ENDDO

  CALL apply_laplacian( phi, nabla2_phi )
  r(:) = rho(:) - nabla2_phi(:)
  p(:) = r(:)

  !
  rsold = ddot(Nbasis,r,1,r,1)
  !WRITE(*,*) 'rsold = ', rsold

  DO iter=1,1000
    CALL apply_laplacian( p, nabla2_phi )
    !
    !alpha = rsold/dot_product(p,nabla2_phi)
    alpha = rsold/ddot( Nbasis,p,1, nabla2_phi,1 )
    !
    phi(:) = phi(:) + alpha*p(:)
    !
    r(:) = r(:) - alpha*nabla2_phi(:)
    !
    !rsnew = dot_product(r,r)
    rsnew = ddot( Nbasis, r,1, r,1 )
    !WRITE(*,*) 'cg-poisson = ', sqrt(rsnew)
    !
    IF(sqrt(rsnew) < 1.d-9) THEN
      WRITE(*,'(/1x,A,I5,E18.10)') 'Convergence in solve_poisson_cg in iter:', iter, sqrt(rsnew)
      EXIT
    ENDIF
    p(:) = r(:) + (rsnew/rsold)*p(:)
    rsold = rsnew
  ENDDO

  DEALLOCATE( r, p, nabla2_phi )
END SUBROUTINE


