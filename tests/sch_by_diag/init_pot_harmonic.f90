!--------------------------------------------
SUBROUTINE init_pot_harmonic( omega, center )
!--------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     lingrid => LF3d_lingrid
  USE m_hamiltonian, ONLY : V_ps_loc
  IMPLICIT NONE
  ! Arguments
  REAL(8) :: omega
  REAL(8) :: center(3)
  ! local
  INTEGER :: ip
  REAL(8) :: dr(3)

  WRITE(*,*) 'center = ', center(:)
  DO ip = 1, Npoints
    dr(:) = lingrid(:,ip) - center(:)
    V_ps_loc(ip) = 0.5d0*omega**2 * ( dr(1)**2 + dr(2)**2 + dr(3)**2 )
  ENDDO
END SUBROUTINE
