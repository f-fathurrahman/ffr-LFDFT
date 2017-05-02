!! PURPOSE
!!
!!   This subroutine solves Poisson equation using preconditioned
!!   conjugate gradient algorithm.
!!
!! AUTHOR
!!
!!   Fadjar Fathurrahman
!!
!! NOTES
!!
!!   The input `rho` will be multiplied by -4*pi.
!!   The output is given in `phi`.

SUBROUTINE Poisson_solve_pcg( rho, phi )
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  IMPLICIT NONE
  !
  REAL(8) :: rho(Npoints), phi(Npoints)
  REAL(8), ALLOCATABLE :: r(:), v(:), z(:) ! residual
  REAL(8), ALLOCATABLE :: nabla2_phi(:)
  INTEGER :: iter, NmaxIter
  REAL(8) :: c, d, omega
  LOGICAL :: conv
  !
  REAL(8) :: ddot

  ALLOCATE( r(Npoints), v(Npoints), nabla2_phi(Npoints), z(Npoints) )

  NmaxIter = Npoints/100

  CALL op_nabla2( phi, nabla2_phi )
  r(:) = -4.d0*PI*rho(:) - nabla2_phi(:)  ! NOTICE that we multiply rho by -4*pi

  CALL prec_ilu0( r, z )
  !z(:) = 2.d0*z(:)

  v(:) = z(:)

  !
  c = ddot( Npoints, r,1, z, 1 )
  
  conv = .FALSE.

  DO iter = 1,NmaxIter
    CALL op_nabla2( v, z )
    !
    omega = c/ddot( Npoints, v,1, z,1 )
    !
    phi(:) = phi(:) + omega*v(:)
    !
    r(:) = r(:) - omega*z(:)
    !
    CALL prec_ilu0( r, z )
    !z(:) = 2.d0*z(:)
    d = ddot( Npoints, r,1, z,1 )
    !
    WRITE(*,'(1x,A,I8,E18.10)') 'rconv = ', iter, sqrt(d)
    !
    IF(sqrt(d) < 1.d-10) THEN
    !IF(rsnew < 1.d-10) THEN
      !WRITE(*,*) 'Convergence in Poisson_solve_cg: iter, sqrt(rsnew):', iter, sqrt(rsnew)
      conv = .TRUE.
      EXIT
    ENDIF
    v(:) = z(:) + (d/c)*v(:)
    c = d
  ENDDO

  IF( .NOT. conv ) WRITE(*,*) 'No convergence in Poisson_solve_cg'

  DEALLOCATE( r, v, nabla2_phi, z )
END SUBROUTINE

