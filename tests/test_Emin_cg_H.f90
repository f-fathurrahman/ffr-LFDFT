! Do an SCF calculation using CG method to minimize
! Kohn-Sham energy functional

PROGRAM test_Emin_cg_H

  USE m_constants, ONLY: PI

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol, &
                     Gv => LF3d_Gv

  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs

  USE m_hamiltonian, ONLY : V_ps_loc

  USE m_atoms, ONLY : Nspecies, Natoms, atm2species, &
                      atpos => AtomicCoords, &
                      strf => StructureFactor, &
                      Zv => AtomicValences
  IMPLICIT NONE
  !
  INTEGER :: ist, ip
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  !TYPE(hgh_t) :: ps
  REAL(8) :: LL(3)
  
  NN = (/ 55, 55, 55 /)
  
  LL(:) = (/ 16.d0, 16.d0, 16.d0 /)
  AA(:) = (/ 0.d0, 0.d0, 0.d0 /)
  BB(:) = LL(:)

  CALL init_LF3d_p( NN, AA, BB )
  CALL info_LF3d()

  CALL init_atoms_xyz( '../structures/H.xyz' )
  Zv(1) = 1.d0
  CALL info_atoms()

  CALL init_strfact()

  ! At this point we can already calculate E_nn
  CALL calc_Ewald()

  ! Set up potential
  CALL alloc_hamiltonian()

  CALL init_nabla2_sparse()
  CALL init_ilu0_prec()

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

  !CALL KS_solve_Emin_cg( 3.d-5, 200, .FALSE. )
  !CALL KS_solve_Emin_cg( 3.d-4, 200, .FALSE. )
  CALL KS_solve_Emin_pcg( 3.d-5, 200, .TRUE. )

  CALL info_energies()

  DEALLOCATE( evecs, evals )
  DEALLOCATE( Focc )
  
  CALL dealloc_nabla2_sparse()
  CALL dealloc_ilu0_prec()
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()

END PROGRAM
