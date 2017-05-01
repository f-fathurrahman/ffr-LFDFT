PROGRAM do_Emin_pcg

  USE m_options, ONLY : FREE_NABLA2
  USE m_PsPot, ONLY : PsPot_Dir
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol, &
                     lingrid => LF3d_lingrid, &
                     xyz2lin => LF3d_xyz2lin
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  USE m_hamiltonian, ONLY : V_ps_loc, V_Hartree

  IMPLICIT NONE 
  INTEGER :: Narg
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  CHARACTER(64) :: filexyz, arg_N
  INTEGER :: ip, ist, N_in, ix, iy, iz

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

  !CALL shift_atoms()

  CALL info_atoms()
  CALL info_PsPot()
  CALL info_LF3d()

  ! Initialize occupation numbers
  CALL init_states()

  CALL init_strfact()

  CALL calc_Ewald()

  CALL alloc_hamiltonian()

  !CALL init_V_ps_loc_G()
  !CALL init_V_coul_G()
  CALL init_V_ps_loc_G_interp()
  !CALL init_V_coul_G_interp()

  CALL init_nabla2_sparse()
  CALL init_ilu0_prec()

  IF( FREE_NABLA2 ) THEN 
    CALL dealloc_nabla2_sparse()
  ENDIF 

  ALLOCATE( evecs(Npoints,Nstates), evals(Nstates) )

  DO ist = 1, Nstates
    DO ip = 1, Npoints
      CALL random_number( evecs(ip,ist) )
    ENDDO
  ENDDO
  CALL orthonormalize( Nstates, evecs )
  CALL ortho_check( Npoints, Nstates, dVol, evecs )

  CALL KS_solve_Emin_pcg( 3.d-5, 1000, .FALSE. )
  !CALL KS_solve_Emin_pcg( 3.d-5, 1000, .TRUE. )

  CALL info_energies()

  iy = 1
  iz = 1
  WRITE(*,*) 'iy iz = ', iy, iz
  DO ix = 1,NN(1)
    ip = xyz2lin(ix,iy,iz)
    WRITE(N_in,'(3F22.12)') lingrid(1,ip), V_ps_loc(ip), V_Hartree(ip)
  ENDDO 

  !
  DEALLOCATE( evecs, evals )
  DEALLOCATE( Focc )
  CALL dealloc_nabla2_sparse()
  CALL dealloc_ilu0_prec()
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()
  CALL dealloc_PsPot()
  CALL dealloc_atoms()

END PROGRAM

