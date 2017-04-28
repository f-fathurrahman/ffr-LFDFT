PROGRAM do_print_V_ps_loc

  USE m_PsPot, ONLY : PsPot_Dir
  USE m_LF3d, ONLY : dVol => LF3d_dVol, &
                     grid_x => LF3d_grid_x, &
                     xyz2lin => LF3d_xyz2lin
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

  CALL init_V_ps_loc_G()

  iy = 28
  iz = 28
  DO ix = 1,NN(1)
    ip = xyz2lin(ix,iy,iz)
    WRITE(*,*) lingrid(:,ip)
  ENDDO 

  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()
  CALL dealloc_PsPot()
  CALL dealloc_atoms()

END PROGRAM

