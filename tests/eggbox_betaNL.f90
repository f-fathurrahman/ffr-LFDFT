#include "init_betaNL_a.f90"

PROGRAM eggbox_betaNL

  USE m_constants, ONLY : Ry2eV
  USE m_options, ONLY : FREE_NABLA2
  USE m_PsPot, ONLY : PsPot_Dir, NbetaNL
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  !
  USE m_atoms, ONLY : atpos => AtomicCoords
  USE m_hamiltonian, ONLY : V_ps_loc, Rhoe, V_ps_loc_long
  USE m_energies, ONLY : E_ps_loc
  USE m_LF3d, ONLY : dVol => LF3d_dVol
  USE m_options, ONLY : T_PRINT_INTEG_RHO, I_ALG_DIAG, ETHR_EVALS
  USE m_constants, ONLY : ANG2BOHR
  USE m_grid_atom_cube, ONLY : Npoints_a, dVol_a
  !
  IMPLICIT NONE 
  INTEGER :: Narg
  INTEGER :: NN(3)
  REAL(8) :: hh(3), AA(3), BB(3)
  CHARACTER(64) :: filexyz, arg_N
  INTEGER :: ip, ist, N_in
  INTEGER :: iargc  ! pgf90 
  INTEGER :: tstart, counts_per_second, tstop
  CHARACTER(1) :: typ
  REAL(8) :: center(3)
  REAL(8), ALLOCATABLE :: V_short_a(:), Rhoe_a(:)

  CALL system_clock( tstart, counts_per_second )

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
  atpos(:,:) = atpos(:,:)/ANG2BOHR  

  ! Override PsPot_Dir
  PsPot_Dir = '../HGH/'
  CALL init_PsPot()

  typ = 'p'
  IF( typ == 's' ) THEN 
    NN = (/ N_in, N_in, N_in /)
    hh(:) = (/1.d0, 1.d0, 1.d0/)*(16.d0/(NN(1)-1))
    CALL init_LF3d_sinc( NN, hh )
  ELSE ! typ == 'p'
    NN = (/ N_in, N_in, N_in /)
    AA = (/ 0.d0, 0.d0, 0.d0 /)
    BB = (/ 16.d0, 16.d0, 16.d0 /)
    CALL init_LF3d_p( NN, AA, BB )
  ENDIF 

  CALL info_atoms()
  CALL info_PsPot()
  CALL info_LF3d()

  center(:) = atpos(:,1)
  CALL init_grid_atom_cube( center, 3.0d0, 55 ) 

  CALL init_betaNL_a()

  CALL dealloc_LF3d()
  CALL dealloc_PsPot()
  CALL dealloc_atoms()

  CALL system_clock( tstop )

  WRITE(*,*)
  WRITE(*,*) 'Total elapsed time: ', dble(tstop - tstart)/counts_per_second, ' second.'
  WRITE(*,*)

END PROGRAM 

