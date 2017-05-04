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
  INTEGER, ALLOCATABLE :: btype(:)
  INTEGER :: dav_iter, notcnv
  REAL(8) :: ethr
  INTEGER :: ist, ip
  INTEGER :: NN(3)
  REAL(8) :: hh(3)
  REAL(8) :: Ekin, Epot, Etot
  !
  REAL(8) :: ddot
  
  NN = (/ 25, 25, 25 /)
  hh = (/ 0.3d0, 0.3d0, 0.3d0 /)

  CALL init_LF3d_sinc( NN, hh )

  CALL info_LF3d()

  ! Set up potential
  CALL alloc_hamiltonian()

  CALL init_nabla2_sparse()
  CALL init_ilu0_prec()

  CALL init_V_ps_loc_harmonic( 2.d0, (/0.d0,0.d0,0.d0/) )
  WRITE(*,*) 'sum(V_ps_loc) = ', sum(V_ps_loc)

  Nstates = 4
  ALLOCATE( Focc(Nstates) )
  Focc(:) = 1.d0

  ALLOCATE( evecs(Npoints,Nstates), evals(Nstates) )

  ALLOCATE( btype(Nstates) )
  btype(:) = 1 ! all bands are occupied
  DO ist = 1, Nstates
    DO ip = 1, Npoints
      CALL random_number( evecs(ip,ist) )
    ENDDO
  ENDDO
  CALL ortho_gram_schmidt( evecs, Npoints, Npoints, Nstates )
  DO ist = 1, Nstates
    WRITE(*,*) ist, ddot( Npoints, evecs(:,ist), 1, evecs(:,ist), 1 )
  ENDDO 
  !evecs(:,:) = evecs(:,:)/sqrt(dVol)

  !ethr = 1.d-1
  !CALL diag_davidson_qe( Npoints, Nstates, 4*Nstates, evecs, ethr, &
  !                       evals, btype, notcnv, dav_iter )
  !WRITE(*,*) 'dav_iter = ', dav_iter
 
  !ethr = 1.d-3
  !CALL diag_davidson_qe( Npoints, Nstates, 4*Nstates, evecs, ethr, &
  !                       evals, btype, notcnv, dav_iter )
  !WRITE(*,*) 'dav_iter = ', dav_iter
  
  ethr = 1.d-6
  !CALL diag_davidson_qe( Npoints, Nstates, 4*Nstates, evecs, ethr, &
  !                       evals, btype, notcnv, dav_iter )
  !CALL diag_davidson( evals, evecs, ethr )
  CALL diag_lobpcg( Nstates, evals, evecs )

  WRITE(*,*) 'dav_iter = ', dav_iter
  
  DO ist = 1, Nstates
    WRITE(*,'(1x,I4,F18.10)') ist, evals(ist)
  ENDDO

  DO ist = 1, Nstates
    WRITE(*,*) ist, ddot( Npoints, evecs(:,ist), 1, evecs(:,ist), 1 )
  ENDDO 

  evecs(:,:) = evecs(:,:)/sqrt(dVol)

  !CALL calc_Energies( evecs, Ekin, Epot, Etot )
  !WRITE(*,*) 'Ekin = ', Ekin
  !WRITE(*,*) 'Epot = ', Epot
  !WRITE(*,*) 'Etot = ', Etot

  DEALLOCATE( evecs, evals )
  
  111 WRITE(*,*) 'Deallocating ...'
  
  DEALLOCATE( Focc )
  CALL dealloc_nabla2_sparse()
  CALL dealloc_ilu0_prec()
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()

END PROGRAM


