! efefer 5 January 2016

MODULE gbl_poisson
  USE m_LF3d
  IMPLICIT NONE
  INTEGER :: Nx, Ny, Nz
  REAL(8) :: Lx, Ly, Lz
  TYPE(LF3d_t) :: LF3d
END MODULE

! apply preconditioner to input vector v
SUBROUTINE precond(N, v, Kv)
  USE gbl_poisson, ONLY : LF3d
  IMPLICIT NONE
  INTEGER :: N
  REAL(8) :: v(N)
  REAL(8) :: Kv(N)
  !
  INTEGER :: i, j, k, ip, Nx,Ny,Nz

  Nx = LF3d%LFx%N
  Ny = LF3d%LFy%N
  Nz = LF3d%LFz%N
  !
  ! N should be the same as LF3d%N
  DO ip=1,N
    i = LF3d%lin2xyz(1,ip)
    j = LF3d%lin2xyz(2,ip)
    k = LF3d%lin2xyz(3,ip)
    !
    !Kv(ip) = v(ip) / LF3d%LFx%D2jl(i,i) &
    !                + v(ip) / LF3d%LFy%D2jl(j,j) &
    !                + v(ip) / LF3d%LFz%D2jl(k,k)
    ! probably this is the correct one
    Kv(ip) = v(ip) / ( LF3d%LFx%D2jl(i,i) + LF3d%LFy%D2jl(j,j) + LF3d%LFz%D2jl(k,k)  )
  ENDDO

END SUBROUTINE


! apply Laplacian to input vector v
SUBROUTINE apply_laplacian(N, v, nabla2_v)
  USE gbl_poisson
  IMPLICIT NONE
  INTEGER :: N
  REAL(8) :: v(N)
  REAL(8) :: nabla2_v(N)
  !
  INTEGER :: i, j, k, ip, ii, jj, kk
  !
  ! N should be the same as LF3d%N
  DO ip=1,N
    i = LF3d%lin2xyz(1,ip)
    j = LF3d%lin2xyz(2,ip)
    k = LF3d%lin2xyz(3,ip)
    !
    nabla2_v(ip) = 0.d0
    !
    DO ii=1,Nx
      nabla2_v(ip) = nabla2_v(ip) + LF3d%LFx%D2jl(ii,i)*v(LF3d%xyz2lin(ii,j,k))
    ENDDO
    !
    DO jj=1,Ny
      nabla2_v(ip) = nabla2_v(ip) + LF3d%LFy%D2jl(jj,j)*v(LF3d%xyz2lin(i,jj,k))
    ENDDO
    !
    DO kk=1,Nz
      nabla2_v(ip) = nabla2_v(ip) + LF3d%LFz%D2jl(kk,k)*v(LF3d%xyz2lin(i,j,kk))
    ENDDO
  ENDDO

END SUBROUTINE

!------------------------------------------
SUBROUTINE solve_poisson_pcg( N, rho, phi )
!------------------------------------------
  USE gbl_poisson
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
  
  ! random guess of phi
  DO ip=1,N
    !CALL random_number( phi(ip) )
    phi(ip) = 0.d0
  ENDDO

  CALL apply_laplacian( N, phi, nabla2vec )
  r(:) = rho(:) - nabla2vec(:)
  !CALL precond( N, r, z )
  z(:) = r(:)
  p(:) = z(:)

  r_old(:) = r(:)
  z_old(:) = z(:)
  p_old(:) = p(:)

  conv = .FALSE.
  DO iter = 1, 1000
    CALL apply_laplacian( N, p, nabla2vec )
    alpha_k = ddot( N, r,1, z,1 )/ddot( N, p,1, nabla2vec,1 )
    !
    phi(:) = phi(:) + alpha_k*p(:)
    !
    r(:) = r(:) - alpha_k*nabla2vec(:)
    norm_res = sqrt( ddot( N, r,1, r,1 ) )
    !
    WRITE(*,*) 'iter, norm_res = ', iter, norm_res
    IF( norm_res < 5.d-10 ) THEN
      WRITE(*,*) 'Convergence achieved in solve_poisson_pcg: N, iter', N, iter
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

!-----------------------------------------
SUBROUTINE solve_poisson_cg( N, rho, phi )
!-----------------------------------------
  USE gbl_poisson
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: rho(N), phi(N)
  !
  REAL(8), ALLOCATABLE :: r(:), p(:) ! residual
  REAL(8), ALLOCATABLE :: nabla2_phi(:)
  INTEGER :: ip, iter, iterTrue
  REAL(8) :: rsold, rsnew, alpha, deltars
  LOGICAL :: conv
  !
  REAL(8) :: ddot

  ALLOCATE( r(N), p(N), nabla2_phi(N) )

  !
  DO ip=1,N
    !CALL random_number( phi(ip) )
    phi(ip) = 0.d0
  ENDDO

  CALL apply_laplacian( N, phi, nabla2_phi )
  r(:) = rho(:) - nabla2_phi(:)
  p(:) = r(:)

  !
  rsold = ddot(N, r,1, r,1)
  
  conv = .FALSE.
  iter = 0
  iterTrue = 0
  cgloop: DO
    iter = iter + 1
    iterTrue = iterTrue + 1
    CALL apply_laplacian( N, p, nabla2_phi )
    !
    alpha = rsold/ddot(N, p,1, nabla2_phi,1)
    !
    phi(:) = phi(:) + alpha*p(:)
    !
    r(:) = r(:) - alpha*nabla2_phi(:)
    !
    rsnew = ddot(N, r,1, r,1)
    deltars = rsold - rsnew
    WRITE(*,*) 'conv = ', sqrt(rsnew)
    !IF( deltars < 0.d0 ) THEN
    !  WRITE(*,*) 'Bad convergence, restarting CG', deltars
    !  iter = 0
    !  CYCLE cgloop
    !ENDIF
    !
    IF(sqrt(rsnew) < 5.d-10) THEN
    !IF(rsnew < 1.d-10) THEN
      WRITE(*,*) 'Convergence achieved in solve_poisson_cg: N, iter', N, iterTrue
      conv = .TRUE.
      EXIT
    ENDIF
    p(:) = r(:) + (rsnew/rsold)*p(:)
    rsold = rsnew
    IF( iter >= 1000 ) THEN
      EXIT cgloop
    ENDIF
  ENDDO cgloop

  IF( .NOT. conv ) WRITE(*,*) 'No convergence in solve_poisson_cg, N', N

  DEALLOCATE( r, p, nabla2_phi )
END SUBROUTINE


!------------------------------------------------------------------------------
PROGRAM t_poisson
!------------------------------------------------------------------------------
  USE m_constants
  USE m_LF3d
  USE gbl_poisson
  IMPLICIT NONE
  !

  REAL(8), ALLOCATABLE :: rho(:) ! density
  REAL(8), ALLOCATABLE :: phi(:) ! potential
  !
  REAL(8) :: sigma1, sigma2, r, x0, y0, z0, deltaV
  INTEGER :: ip
  REAL(8) :: Uana, Unum
  !
  INTEGER :: nread
  CHARACTER(20) :: method
  CHARACTER(20) :: buffer
  !
  REAL(8) :: error
  REAL(8), ALLOCATABLE :: nabla2phi(:)

  CALL getarg(1,buffer)
  READ(buffer,*) nread
  !
  CALL getarg(2,method)

  Nx = Nread
  Ny = Nread
  Nz = Nread
  !
  Lx = 16.d0
  Ly = 16.d0
  Lz = 16.d0
  !
  CALL init_LF3d_p( LF3d, (/Nx,Ny,Nz/), (/0.d0,0.d0,0.d0/), (/Lx,Ly,Lz/) )
  
  ALLOCATE( rho(Nx*Ny*Nz) )
  ALLOCATE( phi(Nx*Ny*Nz) )

  ! center of the box
  x0 = Lx/2.d0
  y0 = Ly/2.d0
  z0 = Lz/2.d0
  ! Initialize
  sigma1 = 0.75d0
  sigma2 = 0.50d0
  DO ip = 1, LF3d%N
    !FIXME add some small positive value to avoid zero distance
    r = norm2( LF3d%lingrid(:,ip) - (/x0,y0,z0/) )
    rho(ip) = exp(-r**2/(2*sigma2**2))/(2*pi*sigma2**2)**1.5d0 !- &
              !exp(-r**2/(2*sigma1**2))/(2*pi*sigma1**2)**1.5d0
    !WRITE(*,'(1x,I5,2F18.10)') ip, r, rho(ip)
  ENDDO

  deltaV = LF3d%LFx%h * LF3d%LFy%h * LF3d%LFz%h
  WRITE(*,*) 'deltaV = ', deltaV
  WRITE(*,*) sum( rho(:) )*deltaV

  STOP

  ! Solve Poisson equation
  IF( trim(method) == 'pcg' ) THEN
    CALL solve_poisson_pcg( LF3d%N, -4.d0*PI*rho, phi )
  ELSEIF( trim(method) == 'cg' ) THEN
    CALL solve_poisson_cg( LF3d%N, -4.d0*PI*rho, phi )
  ELSE
    WRITE(*,*) 'ERROR: unrecognized method', method
    STOP
  ENDIF

  ! Test error of the solution for Poisson equation
  ALLOCATE( nabla2phi(LF3d%N) )
  CALL apply_laplacian( LF3d%N, phi, nabla2phi )
  error = 0.d0
  DO ip=1,LF3d%N
    error = error + abs( nabla2phi(ip) + 4.d0*PI*rho(ip) )
  ENDDO
  error = error/LF3d%N
  WRITE(*,'(1x,A,E18.10)') '>>> Average error of Poisson solution:', error

  !
  ! Calculation of Hartree energy
  !
  Unum = 0.5d0*sum( rho(:)*phi(:) )*deltaV
  Uana = ( (1.d0/sigma1 + 1.d0/sigma2)/2.d0 - sqrt(2.d0)/sqrt(sigma1**2 + sigma2**2) )/sqrt(PI)
  WRITE(*,'(1x,A,I10,E18.10)') 'N, Unum = ', LF3d%N, Unum
  WRITE(*,'(1x,A,I10,E18.10)') 'N, Uana = ', LF3d%N, Uana
  WRITE(*,'(1x,A,I10,E18.10)') 'N, Error = ', LF3d%N, abs(Unum-Uana) 

  DEALLOCATE( rho, phi, nabla2phi )
END PROGRAM

