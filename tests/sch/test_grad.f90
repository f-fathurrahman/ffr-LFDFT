PROGRAM test_grad
  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_states, ONLY : Nstates
  USE m_hamiltonian, ONLY : V_ps_loc
  IMPLICIT NONE
  !
  INTEGER :: ist
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  REAL(8), ALLOCATABLE :: v(:,:), g(:,:), Kg(:,:)
  REAL(8) :: Etot, Ekin, Epot
  
  NN = (/ 25, 25, 25 /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 6.d0, 6.d0, 6.d0 /)

  CALL init_LF3d_p( NN, AA, BB )

  CALL info_LF3d()

  ! Set up potential
  ALLOCATE( V_ps_loc(Npoints) )
  CALL init_pot_harmonic( 2.d0, 0.5*(BB-AA) )
  WRITE(*,*) 'sum(V_ps_loc) = ', sum(V_ps_loc)

  Nstates = 4
  ALLOCATE( v(Npoints,Nstates) )
  ALLOCATE( g(Npoints,Nstates) )
  ALLOCATE( Kg(Npoints,Nstates) )

  v(:,:) = 0.d0
  DO ist = 1, Nstates
    v(ist,ist) = 1.d0/sqrt(dVol)
  ENDDO 

  CALL calc_grad( Nstates, v, g )
  CALL prec_G2( Nstates, g, Kg )
  WRITE(*,*) 'sum(g) = ', sum(g)
  WRITE(*,*) 'sum(Kg) = ', sum(Kg)

  CALL calc_Energies( v, Ekin, Epot, Etot )
  WRITE(*,*) Etot, Ekin, Epot

  DEALLOCATE( v )
  DEALLOCATE( g )
  DEALLOCATE( Kg )
  
  deallocate(V_ps_loc)
  CALL dealloc_LF3d()
END PROGRAM 
