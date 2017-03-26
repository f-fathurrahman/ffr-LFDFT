! Do an SCF calculation using CG method to minimize
! Kohn-Sham energy functional

PROGRAM test_Emin_cg_H

  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol, &
                     lingrid => LF3d_lingrid
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  USE m_hamiltonian, ONLY : V_ps_loc
  USE ps_hgh_m, ONLY : hgh_t, &
                       hgh_init, &
                       hgh_process, &
                       hgh_end, &
                       vlocalr_scalar
  IMPLICIT NONE
  !
  INTEGER :: ist, ip
  INTEGER :: NN(3)
  REAL(8), PARAMETER :: LL(3) = (/ 16.d0, 16.d0, 16.d0 /)
  REAL(8) :: center(3), AA(3), BB(3)
  REAL(8), ALLOCATABLE :: dr(:)
  TYPE(hgh_t) :: ps
  
  NN = (/ 55, 55, 55 /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ LL(1), LL(2), LL(3) /)

  CALL init_LF3d_p( NN, AA, BB )

  CALL info_LF3d()

  ! Set up potential
  CALL alloc_hamiltonian()

  CALL hgh_init( ps, 'tests/pseudo_HGH/HGH/H.hgh' )
  !CALL hgh_process( ps )
  !CALL hgh_info( ps )

  ALLOCATE( dr(Npoints) )
  center(:) = 0.5d0*LL(:)
  !center(:) = 0.d0  ! to test calc_dr_periodic

  CALL calc_dr_periodic( LL, center, Npoints, lingrid, dr )
  CALL init_V_ps_loc_H_hgh_G( Npoints, V_ps_loc )

  !DO ip = 1, Npoints
    !dr(ip) = sqrt( (lingrid(1,ip) - center(1))**2 + &
    !           (lingrid(2,ip) - center(2))**2 + &
    !           (lingrid(3,ip) - center(3))**2 )
    !V_ps_loc(ip) = vlocalr_scalar( dr(ip), ps )
  !ENDDO 

  WRITE(*,*) 'sum(V_ps_loc) = ', sum(V_ps_loc)
  !STOP

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
  CALL test_orthonormal( Npoints, Nstates, dVol, evecs )

  CALL kssolve_Emin_cg( 3.d-5, 100, .FALSE. )

  CALL info_energies()

  DEALLOCATE( evecs, evals )
  DEALLOCATE( Focc )
  
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()

END PROGRAM
