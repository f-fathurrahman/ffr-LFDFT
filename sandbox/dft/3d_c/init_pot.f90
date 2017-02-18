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

  r0(:) = A + (B-A)/2.d0 ! position of H atom

  rc1 = 0.25d0
  rc2 = 0.284d0
  aa  = -1.9287d0
  bb  = 0.3374d0

  WRITE(*,*)
  WRITE(*,*) 'Initializing pseudopotential for H atom'
  WRITE(*,*) '---------------------------------------'
  WRITE(*,'(1x,A,3F10.5)') 'pos = ', r0(:)

  ! TODO Add journal reference for this pseudopotential
  DO ip=1,N**3
    r = norm2( LF%lingrid(:,ip)-r0(:) )
    V(ip) = -1.d0/r * erf( r/rc1 ) + (aa + bb*r**2)*exp(-(r/rc2)**2)
  ENDDO
  WRITE(*,*) 'DEBUG: sum(V) = ', sum(V)
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

  r0(:) = A + (B-A)/2.d0
  WRITE(*,*) 'center = ', r0

  ! N is assumed to be the same as LF$N
  DO ip=1,N**3
    V(ip) = 0.5d0*omega**2 * norm2( LF%lingrid(:,ip) - r0(:) )**2
  ENDDO
  WRITE(*,*) 'DEBUG: sum(Vpot) = ', sum(V)
END SUBROUTINE
