PROGRAM test_sch

  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_states, ONLY : Nstates, &
                       Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  USE m_hamiltonian, ONLY : V_ps_loc
  IMPLICIT NONE
  !
  INTEGER :: ist, ip
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  !
  REAL(8) :: ddot
  
  NN = (/ 35, 35, 35 /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 6.d0, 6.d0, 6.d0 /)

  CALL init_LF3d_c( NN, AA, BB )

  CALL info_LF3d()

  ! Set up potential
  CALL alloc_hamiltonian()

  CALL init_nabla2_sparse()
  CALL init_ilu0_prec()

  CALL init_V_ps_loc_harmonic( 2.d0, 0.5*(BB-AA) )
  WRITE(*,*) 'sum(V_ps_loc) = ', sum(V_ps_loc)

  Nstates = 4
  ALLOCATE( Focc(Nstates) )
  Focc(:) = 1.d0

  ALLOCATE( evecs(Npoints,Nstates), evals(Nstates) )

  DO ist = 1, Nstates
    DO ip = 1, Npoints
      CALL random_number( evecs(ip,ist) )
    ENDDO
  ENDDO
  CALL ortho_gram_schmidt( evecs, Npoints, Npoints, Nstates )

  CALL diag_lobpcg( evals, evecs, 1.0d-4, .TRUE. )

  WRITE(*,*)
  WRITE(*,*) 'Final eigenvalues:'
  DO ist = 1, Nstates
    WRITE(*,'(1x,I4,F18.10)') ist, evals(ist)
  ENDDO

  ! remember to renormalize eigenvectors
  evecs(:,:) = evecs(:,:)/sqrt(dVol)

  WRITE(*,*)
  WRITE(*,*) 'Check normalization:'
  DO ist = 1, Nstates
    WRITE(*,'(1x,I4,F18.10)') ist, ddot( Npoints, evecs(:,ist), 1, evecs(:,ist), 1 )*dVol
  ENDDO 

  DEALLOCATE( evecs, evals )
  DEALLOCATE( Focc )
  CALL dealloc_nabla2_sparse()
  CALL dealloc_ilu0_prec()
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()

END PROGRAM


