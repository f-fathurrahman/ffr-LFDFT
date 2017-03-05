! Do an SCF calculation using CG method to minimize
! Kohn-Sham energy functional
! 
! This will eventually be organized into one subroutine.
PROGRAM test_Emin_cg

  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  USE m_hamiltonian, ONLY : V_ps_loc
  IMPLICIT NONE
  !
  INTEGER :: ist, ip
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  REAL(8), PARAMETER :: mixing_beta = 0.1d0
  
  NN = (/ 25, 25, 25 /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 6.d0, 6.d0, 6.d0 /)

  CALL init_LF3d_p( NN, AA, BB )

  CALL info_LF3d()

  ! Set up potential
  CALL alloc_hamiltonian()

  CALL init_V_ps_loc_harmonic( 2.d0, 0.5*(BB-AA) )

  WRITE(*,*) 'sum(V_ps_loc) = ', sum(V_ps_loc)

  ! Initialize electronic states variables
  Nstates = 4

  ALLOCATE( evecs(Npoints,Nstates), evals(Nstates) )
  ALLOCATE( Focc(Nstates) )

  Focc(:) = 2.d0

  DO ist = 1, Nstates
    DO ip = 1, Npoints
      CALL random_number( evecs(ip,ist) )
    ENDDO
  ENDDO
  CALL orthonormalize( Nstates, evecs )
  CALL test_orthonormal( Npoints, Nstates, dVol, evecs )

  CALL kssolve_Emin_cg( 3.d-5, 100, .FALSE. )

  CALL info_energies()

  DEALLOCATE( evecs, evals )
  DEALLOCATE( Focc )
  
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()

END PROGRAM
