!! PURPOSE:
!!
!!   This subroutine calculates total energy components.
!!
!! AUTHOR:
!!
!!   Fadjar Fathurrahman
!!
!! MODIFIES:
!!   
!!   Global variables defined in module `m_energies`
!!
!! IMPORTANT:
!!
!!   The input `psi` is assumed to be orthonormalized.
!!
SUBROUTINE calc_energies( psi )

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_states, ONLY : Nstates, Focc
  USE m_hamiltonian, ONLY : V_ps_loc, V_Hartree, Rhoe
  USE m_energies
  IMPLICIT NONE
  !
  REAL(8) :: psi(Npoints, Nstates)
  !
  REAL(8), ALLOCATABLE  :: nabla2_psi(:)
  INTEGER :: ist
  REAL(8), ALLOCATABLE :: epsxc(:)
  !
  REAL(8) :: ddot

  ALLOCATE( nabla2_psi(Npoints) )
  ALLOCATE( epsxc(Npoints) )

  E_total   = 0.d0
  E_kinetic = 0.d0
  E_ps_loc  = 0.d0
  E_Hartree = 0.d0
  E_xc      = 0.d0

  ! assume all states are occupied
  DO ist = 1, Nstates
    CALL op_nabla2( psi(:,ist), nabla2_psi(:) )
    E_kinetic = E_kinetic + Focc(ist) * (-0.5d0) * ddot( Npoints, psi(:,ist),1, nabla2_psi(:),1 ) * dVol
  ENDDO

  E_ps_loc = sum( Rhoe(:) * V_ps_loc(:) )*dVol

  E_Hartree = 0.5d0*sum( Rhoe(:) * V_Hartree(:) )*dVol

  CALL excVWN( Npoints, Rhoe, epsxc )
  E_xc = sum( Rhoe(:) * epsxc(:) )*dVol

  E_total = E_kinetic + E_ps_loc + E_Hartree + E_xc + E_nn

  DEALLOCATE( epsxc )
  DEALLOCATE( nabla2_psi )
END SUBROUTINE
