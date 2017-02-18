! rho will be multiplied by -4*pi in this subroutine

SUBROUTINE solve_poisson_cg( rho, phi )
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  IMPLICIT NONE
  !
  REAL(8) :: rho(Npoints), phi(Npoints)
  REAL(8), ALLOCATABLE :: r(:), p(:) ! residual
  REAL(8), ALLOCATABLE :: nabla2_phi(:)
  INTEGER :: ip, iter
  REAL(8) :: rsold, rsnew, alpha
  LOGICAL :: conv

  ALLOCATE( r(Npoints), p(Npoints), nabla2_phi(Npoints) )

  !
  DO ip=1,Npoints
    !CALL random_number( phi(ip) )
    phi(ip) = 0.d0
  ENDDO

  CALL apply_Laplacian( phi, nabla2_phi )
  r(:) = -4.d0*PI*rho(:) - nabla2_phi(:)  ! NOTICE that we multiply rho by -4*pi
  p(:) = r(:)

  !
  rsold = dot_product(r,r)
  WRITE(*,*) 'Initial rsold = ', rsold
  
  conv = .FALSE.
  DO iter = 1,Npoints
    CALL apply_Laplacian( p, nabla2_phi )
    !
    alpha = rsold/dot_product(p,nabla2_phi)
    !
    phi(:) = phi(:) + alpha*p(:)
    !
    r(:) = r(:) - alpha*nabla2_phi(:)
    !
    rsnew = dot_product(r,r)
    WRITE(*,'(1x,A,E18.10)') 'rconv = ', sqrt(rsnew)
    !
    IF(sqrt(rsnew) < 1.d-10) THEN
    !IF(rsnew < 1.d-10) THEN
      WRITE(*,*) 'Convergence in solve_poisson_cg: iter', iter
      conv = .TRUE.
      EXIT
    ENDIF
    p(:) = r(:) + (rsnew/rsold)*p(:)
    rsold = rsnew
  ENDDO

  IF( .NOT. conv ) WRITE(*,*) 'No convergence in solve_poisson_cg'

  DEALLOCATE( r, p, nabla2_phi )
END SUBROUTINE

