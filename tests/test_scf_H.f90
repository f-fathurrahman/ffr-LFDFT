! Do an SCF calculation using Davidson method to diagonalize
! the Kohn-Sham Hamiltonian
! 
! This will eventually be organized into one subroutine.
PROGRAM test_scf
  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol, &
                     lingrid => LF3d_lingrid
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  USE m_hamiltonian, ONLY : V_ps_loc, Rhoe
  USE m_energies, ONLY : Etot => E_total
  USE ps_hgh_m, ONLY : hgh_t, &
                       hgh_init, &
                       hgh_process, &
                       hgh_end, &
                       vlocalr_scalar
  IMPLICIT NONE
  !
  INTEGER :: ist, ip, iterSCF
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  REAL(8) :: Etot_old, dEtot
  REAL(8), ALLOCATABLE :: Rhoe_old(:)
  REAL(8), PARAMETER :: mixing_beta = 0.1d0
  REAL(8), PARAMETER :: LL(3) = (/ 16.d0, 16.d0, 16.d0 /)
  REAL(8) :: center(3)
  REAL(8) :: dr
  TYPE(hgh_t) :: ps
  
  NN = (/ 55, 55, 55 /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ LL(1), LL(2), LL(3) /)

  CALL init_LF3d_p( NN, AA, BB )

  CALL info_LF3d()

  ! Set up potential
  CALL alloc_hamiltonian()

  CALL init_V_ps_loc_H_hgh_G( Npoints, V_ps_loc )

  WRITE(*,*) 'sum(V_ps_loc) = ', sum(V_ps_loc)

  ! Initialize electronic states variables
  Nstates = 1

  ALLOCATE( evecs(Npoints,Nstates), evals(Nstates) )
  ALLOCATE( Focc(Nstates) )

  Focc(:) = 1.d0

  DO ist = 1, Nstates
    DO ip = 1, Npoints
      CALL random_number( evecs(ip,ist) )
    ENDDO
  ENDDO
  CALL orthonormalize( Nstates, evecs )
  CALL ortho_check( Npoints, Nstates, dVol, evecs )

  CALL calc_rhoe( evecs, Focc )
  CALL update_potentials()

  ALLOCATE( Rhoe_old(Npoints) )

  Etot_old = 0.d0
  Rhoe_old(:) = Rhoe(:)

  DO iterSCF = 1, 100

    CALL Sch_solve_diag()
    CALL calc_energies( evecs ) ! not updating potentials

    dEtot = abs(Etot - Etot_old)

    WRITE(*,*)
    WRITE(*,*) 'SCF iter', iterSCF, Etot, dEtot

    IF( dEtot < 1d-6) THEN 
      WRITE(*,*)
      WRITE(*,*) 'SCF converged!!!'
      EXIT 
    ENDIF 

    CALL calc_rhoe( evecs, Focc )

    Rhoe(:) = 0.7d0*Rhoe(:) + 0.3d0*Rhoe_old(:)

    WRITE(*,'(1x,A,F18.10)') 'After mix: integRho = ', sum(Rhoe)*dVol

    CALL update_potentials()

    Etot_old = Etot
    Rhoe_old(:) = Rhoe(:)
  ENDDO 

  CALL info_energies()

  DEALLOCATE( evecs, evals )
  DEALLOCATE( Focc )
  
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()

END PROGRAM
