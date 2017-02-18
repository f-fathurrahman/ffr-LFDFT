SUBROUTINE init_potential(pot_type)
  USE m_globals, ONLY : Vpot
  IMPLICIT NONE
  CHARACTER(8) :: pot_type

  CALL init_pspot_H( Vpot )

END SUBROUTINE


! Calculate Hartree potential and energy
SUBROUTINE calc_hartree()
  USE m_globals, ONLY : Rho, Vhartree, Ehartree, PI, LF
  IMPLICIT NONE

  CALL solve_poisson_cg( -4.d0*PI*Rho, Vhartree)

  Ehartree = 0.5d0*sum(Rho*Vhartree)*LF%LFx%h*LF%LFy%h*LF%LFz%h
  WRITE(*,*) 'Ehartree = ', Ehartree
END SUBROUTINE



SUBROUTINE init_pot_softCoulomb( Z, V )
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

  ! N is assumed to be the same as LF%N
  DO ip=1,N**3
    V(ip) = -Z/sqrt( 1 + norm2( LF%lingrid(:,ip) - r0(:) )**2 )
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

  r0(:) = 1.d-8 ! position of H atom

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



!-------------------------------------
SUBROUTINE init_pspot_H2( V )
!-------------------------------------
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : N, LF, A, B
  IMPLICIT NONE
  ! Arguments
  REAL(8) :: V(N**3)
  ! local
  INTEGER, PARAMETER :: NATOM=2
  INTEGER :: ip, ia
  REAL(8) :: r0(3,NATOM)
  !
  REAL(8) :: rc1, rc2, aa, bb, r


  r0(:,1) = (/ -0.5d0, 0.d0, 0.d0 /)
  r0(:,2) = (/ +0.5d0, 0.d0, 0.d0 /)

  rc1 = 0.25d0
  rc2 = 0.284d0
  aa = -1.9287d0
  bb = 0.3374d0

  ! FIXME Only pure radial potential
  V(:) = 0.d0
  DO ia = 1, NATOM
    DO ip = 1, N**3
      r = norm2( LF%lingrid(:,ip)-r0(:,ia) ) + 1.d-8
      V(ip) = V(ip) + -1.d0/r * erf( r/rc1 ) + (aa + bb*r**2)*exp(-(r/rc2)**2)
    ENDDO
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

  r0(:) = 1.d-8

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

  r0(:) = 0.d0

  ! N is assumed to be the same as LF$N
  DO ip=1,N**3
    V(ip) = 0.5d0*omega**2 * norm2( LF%lingrid(:,ip) - r0(:) )**2
  ENDDO
END SUBROUTINE


