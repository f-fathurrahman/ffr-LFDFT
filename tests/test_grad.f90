PROGRAM test_grad

  USE m_PsPot, ONLY : PsPot_Dir
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, dVol => LF3d_dVol
  USE m_states, ONLY : Nstates, Focc
  USE m_energies, ONLY : E_total
  IMPLICIT NONE 
  INTEGER :: Narg
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  CHARACTER(64) :: filexyz, arg_N
  INTEGER :: ip, ist, N_in
  INTEGER :: iargc  ! pgf90 
  REAL(8), ALLOCATABLE :: psi0(:,:), psi(:,:), dW(:,:)
  REAL(8), ALLOCATABLE :: g0(:,:)
  REAL(8) :: E0, delta, E1, dE
  INTEGER :: iexp

  Narg = iargc()
  IF( Narg /= 2 ) THEN 
    WRITE(*,*) 'ERROR: exactly two arguments must be given:'
    WRITE(*,*) '       N and path to structure file'
    STOP 
  ENDIF 

  CALL getarg( 1, arg_N )
  READ(arg_N, *) N_in

  CALL getarg( 2, filexyz )

  CALL init_atoms_xyz(filexyz)

  ! Override PsPot_Dir
  PsPot_Dir = '../HGH/'
  CALL init_PsPot()

  !
  NN = (/ N_in, N_in, N_in /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)
  !AA = (/ -8.d0, -8.d0, -8.d0 /)
  !BB = (/  8.d0,  8.d0,  8.d0 /)
  CALL init_LF3d_p( NN, AA, BB )

  CALL info_atoms()
  CALL info_PsPot()
  CALL info_LF3d()

  CALL init_betaNL()

  ! Initialize occupation numbers
  CALL init_states()

  ! Structure factor, shifted to FFT grid
  CALL init_strfact_shifted()

  ! Ewald energy
  CALL calc_Ewald()

  ! Memory for potentials
  CALL alloc_hamiltonian()

  ! Local pseudopotential
  CALL init_V_ps_loc_G()

  ! Laplacian matrix
  CALL init_nabla2_sparse()


  ALLOCATE( psi(Npoints,Nstates) )
  ALLOCATE( psi0(Npoints,Nstates) )
  ALLOCATE( dW(Npoints,Nstates) )
  ALLOCATE( g0(Npoints,Nstates) )

  ! Initialize to random wavefunction
  DO ist = 1, Nstates
    DO ip = 1, Npoints
      CALL random_number( psi0(ip,ist) )
    ENDDO
  ENDDO
  CALL orthonormalize( Nstates, psi0 )

  CALL calc_rhoe( psi0, Focc )
  CALL update_potentials()
  CALL calc_betaNL_psi( Nstates, psi0 )
  CALL calc_energies( psi0 )
  E0 = E_total
  CALL calc_grad( Nstates, psi0, g0 )

  WRITE(*,*) 'E0 = ', E0

  ! random direction
  DO ist = 1, Nstates
    DO ip = 1, Npoints
      CALL random_number( dW(ip,ist) )
    ENDDO
  ENDDO

  DO iexp = 1,9

    delta = 10.d0**(-iexp)

    WRITE(*,*)
    WRITE(*,'(1x,A,ES18.10)') 'delta = ', delta

    dE = sum( g0 * delta * dW )*dVol*sum(Focc)
    psi = psi0 + delta*dW

    CALL orthonormalize( Nstates, psi )
    CALL calc_rhoe( psi, Focc )
    CALL update_potentials()
    CALL calc_betaNL_psi( Nstates, psi )

    CALL calc_energies( psi )
    E1 = E_total

    WRITE(*,'(1x,A,4F18.10)') 'E0, E1, diff, dE = ', E0, E1, E1-E0, dE
    WRITE(*,*) 'ratio = ', (E1-E0)/dE
  ENDDO 

  !
  DEALLOCATE( Focc )
  CALL dealloc_nabla2_sparse()
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()
  CALL dealloc_PsPot()
  CALL dealloc_atoms()

END PROGRAM

