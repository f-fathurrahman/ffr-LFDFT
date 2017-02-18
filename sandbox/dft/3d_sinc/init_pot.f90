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

  r0(:) = (/0.d0,0.d0,0.d0/)

  ! N is assumed to be the same as LF$N
  DO ip=1,N**3
    V(ip) = 0.5d0*omega**2 * norm2( LF%lingrid(:,ip) - r0(:) )**2
  ENDDO
END SUBROUTINE
