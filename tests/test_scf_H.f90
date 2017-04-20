! Do an SCF calculation using Davidson method to diagonalize
! the Kohn-Sham Hamiltonian
! 
! This will eventually be organized into one subroutine.
PROGRAM test_scf
  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  USE m_hamiltonian, ONLY : V_ps_loc, Rhoe
  USE m_energies, ONLY : Etot => E_total
  USE m_PsPot, ONLY : PsPot_Dir
  USE m_options, ONLY : ethr => DIAG_DAVIDSON_QE_ETHR
  USE m_states, ONLY : Nelectrons
  IMPLICIT NONE
  !
  INTEGER :: ist, ip, iterSCF
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  REAL(8) :: Etot_old, dEtot
  REAL(8) :: dr2
  REAL(8), ALLOCATABLE :: Rhoe_old(:)
  REAL(8), PARAMETER :: mixing_beta = 0.1d0
  REAL(8), PARAMETER :: LL(3) = (/ 16.d0, 16.d0, 16.d0 /)
  REAL(8) :: ddot
  REAL(8) :: integRho
  
  CALL init_atoms_xyz('../structures/LiH.xyz')

  PsPot_Dir = '../HGH/'
  CALL init_PsPot()

  NN = (/ 55, 55, 55 /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ LL(1), LL(2), LL(3) /)

  CALL init_LF3d_p( NN, AA, BB )

  CALL info_LF3d()

  CALL init_states()

  CALL init_strfact()
  CALL calc_Ewald()

  ! Set up potential
  CALL alloc_hamiltonian()

  CALL init_nabla2_sparse()
  CALL init_ilu0_prec()

  CALL init_V_ps_loc_G( )

  WRITE(*,*) 'sum(V_ps_loc) = ', sum(V_ps_loc)

  ALLOCATE( evecs(Npoints,Nstates), evals(Nstates) )

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

  dr2 = 1.d0
  DO iterSCF = 1, 100

    IF( iterSCF==1 ) THEN
      ethr = 1.d-1
    ELSE 
      IF( iterSCF == 2 ) ethr = 1.d-2
      !ethr = min( ethr, 1.d-2*dEtot / max(1.d0,Nelectrons) )
      ethr = ethr/5.d0
      ethr = max( ethr, 1d-13 )
      WRITE(*,'(1x,A,ES18.10)') 'ethr = ', ethr
    ENDIF 

    CALL Sch_solve_diag()
    CALL calc_energies( evecs ) ! not updating potentials

    dEtot = abs(Etot - Etot_old)

    IF( dEtot < 1d-6) THEN 
      WRITE(*,*)
      WRITE(*,*) 'SCF converged!!!'
      EXIT 
    ENDIF 

    WRITE(*,*)
    WRITE(*,'(1x,A,I5,F18.10,2ES18.10)') 'SCF iter', iterSCF, Etot, dEtot, dr2

    CALL calc_rhoe( evecs, Focc )

    Rhoe(:) = 0.5d0*Rhoe(:) + 0.5d0*Rhoe_old(:)
    IF( iterSCF > 2 ) THEN 
      dr2 = sqrt( ddot( Npoints, Rhoe(:)-Rhoe_old(:), 1, Rhoe(:)-Rhoe_old(:), 1 ) )
    ENDIF

    integRho = sum(Rhoe)*dVol
    WRITE(*,'(1x,A,F18.10)') 'After mix: integRho = ', integRho
    IF( abs(integRho - Nelectrons) > 1.0d-6 ) THEN
      WRITE(*,*) 'Rescaling Rho'
      Rhoe(:) = Nelectrons/integRho * Rhoe(:)
      integRho = sum(Rhoe)*dVol
      WRITE(*,'(1x,A,F18.10)') 'After rescaling: integRho = ', integRho
    ENDIF 

    CALL update_potentials()

    Etot_old = Etot
    Rhoe_old(:) = Rhoe(:)
  ENDDO 

  CALL info_energies()

  DEALLOCATE( evecs, evals )
  DEALLOCATE( Focc )
 
  CALL dealloc_atoms()
  CALL dealloc_PsPot()
  CALL dealloc_nabla2_sparse()
  CALL dealloc_ilu0_prec()
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()

END PROGRAM

