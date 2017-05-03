PROGRAM test_Ylm_real

  USE m_PsPot, ONLY : PsPot_Dir
  USE m_LF3d, ONLY : lingrid => LF3d_lingrid, &
                     grid_x => LF3d_grid_x, &
                     grid_y => LF3d_grid_y, &
                     grid_z => LF3d_grid_z, &
                     Npoints => LF3d_Npoints, &
                     LL => LF3d_LL
  USE m_atoms, ONLY : atpos => AtomicCoords, &
                      SpeciesSymbols, atm2species, Natoms
  USE m_constants, ONLY : ANG2BOHR
  IMPLICIT NONE 
  INTEGER :: Narg
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  CHARACTER(64) :: filexyz, arg_N
  INTEGER :: ip, N_in, ia
  REAL(8) :: dr, dr_vec(3)
  REAL(8) :: LatVecs(3,3), origin(3)
  REAL(8), ALLOCATABLE :: fun(:)
  REAL(8) :: Ylm_real

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
  ! so that coord given in xyz file is in bohr
  !atpos(:,:) = atpos(:,:)/ANG2BOHR  

  ! Override PsPot_Dir
  PsPot_Dir = '../HGH/'
  CALL init_PsPot()

  !
  NN = (/ N_in, N_in, N_in /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)
  CALL init_LF3d_p( NN, AA, BB )

  CALL info_atoms()
  CALL info_PsPot()
  CALL info_LF3d()

  ! Initialize occupation numbers
  CALL init_states()

  CALL init_strfact_shifted()

  CALL calc_Ewald()

  CALL alloc_hamiltonian()

  CALL init_V_ps_loc_G()

  LatVecs(:,:) = 0.d0
  LatVecs(1,1) = BB(1) - AA(1)
  LatVecs(2,2) = BB(2) - AA(2)
  LatVecs(3,3) = BB(3) - AA(3)
  CALL xsf_struct( LatVecs, Natoms, atpos, SpeciesSymbols, atm2species, 444 )
  origin(1) = 0.5d0*( grid_x(2) - grid_x(1) )
  origin(2) = 0.5d0*( grid_y(2) - grid_y(1) )
  origin(3) = 0.5d0*( grid_z(2) - grid_z(1) )
  !CALL xsf_fast_datagrid_3d(V_ps_loc, NN(1), NN(2), NN(3), NN(1), NN(2), NN(3), &
  !          origin, LatVecs, 444)

  ALLOCATE( fun(Npoints) )
  ia = 1
  DO ip = 1,Npoints
    CALL calc_dr_periodic_1pnt( LL, atpos(:,ia), lingrid(:,ip), dr_vec )
    dr = sqrt( dr_vec(1)**2 + dr_vec(2)**2 + dr_vec(3)**2 )
    fun(ip) = exp(-0.5d0*dr) * Ylm_real( 3,1, dr_vec )
  ENDDO 
  CALL xsf_fast_datagrid_3d(fun, NN(1), NN(2), NN(3), NN(1), NN(2), NN(3), &
            origin, LatVecs, 444)


  DEALLOCATE( fun )

END PROGRAM 

