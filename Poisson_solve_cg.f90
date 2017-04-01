!! PURPOSE
!!
!!   This subroutine solves Poisson equation using conjugate gradient
!!   algorithm.
!!
!! AUTHOR
!!
!!   Fadjar Fathurrahman
!!
!! NOTES
!!
!!   The input `rho` will be multiplied by -4*pi.
!!   The output is given in `phi`.

SUBROUTINE Poisson_solve_cg( rho, phi )
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  IMPLICIT NONE
  !
  REAL(8) :: rho(Npoints), phi(Npoints)
  REAL(8), ALLOCATABLE :: r(:), p(:) ! residual
  REAL(8), ALLOCATABLE :: nabla2_phi(:)
  INTEGER :: iter, NmaxIter
  REAL(8) :: rsold, rsnew, alpha
  LOGICAL :: conv
  !
  REAL(8) :: ddot

  ALLOCATE( r(Npoints), p(Npoints), nabla2_phi(Npoints) )

  NmaxIter = Npoints/100

  !
  !DO ip=1,Npoints
    !CALL random_number( phi(ip) )
  !  phi(ip) = 0.d0
  !ENDDO

  CALL op_nabla2( phi, nabla2_phi )
  r(:) = -4.d0*PI*rho(:) - nabla2_phi(:)  ! NOTICE that we multiply rho by -4*pi
  p(:) = r(:)

  !
  rsold = ddot( Npoints, r,1, r, 1 )
  
  conv = .FALSE.

  DO iter = 1,NmaxIter
    CALL op_nabla2( p, nabla2_phi )
    !
    alpha = rsold/ddot( Npoints, p,1, nabla2_phi,1 )
    !
    phi(:) = phi(:) + alpha*p(:)
    !
    r(:) = r(:) - alpha*nabla2_phi(:)
    !
    rsnew = ddot( Npoints, r,1, r,1 )
    !WRITE(*,'(1x,A,E18.10)') 'rconv = ', sqrt(rsnew)
    !
    IF(sqrt(rsnew) < 1.d-10) THEN
    !IF(rsnew < 1.d-10) THEN
      WRITE(*,*) 'Convergence in Poisson_solve_cg: iter, sqrt(rsnew):', iter, sqrt(rsnew)
      conv = .TRUE.
      EXIT
    ENDIF
    p(:) = r(:) + (rsnew/rsold)*p(:)
    rsold = rsnew
  ENDDO

  IF( .NOT. conv ) WRITE(*,*) 'No convergence in Poisson_solve_cg'

  DEALLOCATE( r, p, nabla2_phi )
END SUBROUTINE

