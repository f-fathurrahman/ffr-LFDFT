!!
!! PURPOSE:
!!
!!   This subroutine calculates electronic density, given `psi`
!!   (which need not to be Kohn-Sham states) and occupation number `Focc`
!! 
!! AUTHOR:
!!
!!   Fadjar Fathurrahman
!!
!! MODIFY:
!!
!!   Global variable `rhoe`
!!
SUBROUTINE calc_rhoe( Focc, psi )

  USE m_options, ONLY : T_PRINT_INTEG_RHO
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints  , dVol => LF3d_dVol
  USE m_states, ONLY : Nstates_occ, Nstates
  USE m_hamiltonian, ONLY : Rhoe
  IMPLICIT NONE
  REAL(8) :: psi(Npoints,Nstates)
  REAL(8) :: Focc(Nstates)
  INTEGER :: ist

  Rhoe(:) = 0.d0
  ! FIXME: only for FIXED occ ? Not working for metals with fractional occupation ??
  DO ist = 1, Nstates_occ
    Rhoe(:) = Rhoe(:) + Focc(ist) * psi(:,ist) * psi(:,ist)
  ENDDO
  
  IF( T_PRINT_INTEG_RHO ) THEN 
    WRITE(*,*)
    WRITE(*,*) 'Calculating electron density'
    WRITE(*,'(1x,A,F18.10)') 'Integrated electron density:', sum( Rhoe(:) )*dVol
  ENDIF 

END SUBROUTINE 

