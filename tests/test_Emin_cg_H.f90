! Do an SCF calculation using CG method to minimize
! Kohn-Sham energy functional

PROGRAM test_Emin_cg_H

  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol, &
                     lingrid => LF3d_lingrid, &
                     LF3d_TYPE, LF3d_PERIODIC
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  USE m_hamiltonian, ONLY : V_ps_loc
  !USE ps_hgh_m, ONLY : hgh_t, &
  !                     hgh_init, &
  !                     hgh_process, &
  !                     hgh_end, &
  !                     vlocalr_scalar
  IMPLICIT NONE
  !
  INTEGER :: ist, ip
  INTEGER :: NN(3)
  REAL(8) :: center(3), AA(3), BB(3)
  REAL(8), ALLOCATABLE :: dr(:)
  !TYPE(hgh_t) :: ps
  REAL(8) :: LL(3)
  
  NN = (/ 60, 60, 60 /)
  
  LL(:) = (/ 20.d0, 20.d0, 20.d0 /)
  AA(:) = (/ 0.d0, 0.d0, 0.d0 /)
  BB(:) = LL(:)
  !CALL init_LF3d_p( NN, AA, BB )
  center(:) = 0.5d0*LL

  !AA = (/ -8.d0, -8.d0, -8.d0 /)
  !BB = (/  8.d0,  8.d0,  8.d0 /)
  CALL init_LF3d_c( NN, AA, BB )
  !center(:) = (/ 8.d0, 8.d0, 8.d0 /)

  CALL info_LF3d()

  ! Set up potential
  CALL alloc_hamiltonian()
  CALL init_K_diag()

  !CALL hgh_init( ps, 'pseudo_HGH/HGH/H.hgh' )
  !CALL hgh_process( ps )
  !CALL hgh_info( ps )

  ALLOCATE( dr(Npoints) )

  IF ( LF3d_TYPE == LF3d_PERIODIC ) THEN
    CALL init_V_ps_loc_H_hgh_G( Npoints, V_ps_loc )
  ELSE
    CALL calc_dr( center, Npoints, lingrid, dr )
    CALL init_V_ps_loc_H_hgh( Npoints, dr, V_ps_loc )
    !DO ip = 1, Npoints
    !  V_ps_loc(ip) = vlocalr_scalar( dr(ip), ps )
    !ENDDO 
  ENDIF 

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

  CALL KS_solve_Emin_cg( 3.d-5, 200, .FALSE. )
  !CALL KS_solve_Emin_cg( 3.d-4, 200, .FALSE. )

  CALL info_energies()

  DEALLOCATE( evecs, evals )
  DEALLOCATE( Focc )
  
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()

END PROGRAM
