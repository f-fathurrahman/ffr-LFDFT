PROGRAM visualize_gridpts

  USE m_constants, ONLY : Ry2eV, PI
  USE m_options, ONLY : FREE_NABLA2
  USE m_PsPot, ONLY : PsPot_Dir, NbetaNL
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     lingrid => LF3d_lingrid
  !
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  !
  USE m_atoms, ONLY : atpos => AtomicCoords
  USE m_hamiltonian, ONLY : V_ps_loc, Rhoe, V_ps_loc_long
  USE m_energies, ONLY : E_ps_loc
  USE m_LF3d, ONLY : dVol => LF3d_dVol, lingrid => LF3d_lingrid, lin2xyz => LF3d_lin2xyz
  USE m_options, ONLY : T_PRINT_INTEG_RHO, I_ALG_DIAG, ETHR_EVALS
  USE m_constants, ONLY : ANG2BOHR
  USE m_grid_atom_cube, ONLY : Npoints_a, dVol_a, grid_a
  !
  IMPLICIT NONE 
  INTEGER :: Narg
  INTEGER :: NN(3)
  REAL(8) :: hh(3), AA(3), BB(3)
  CHARACTER(64) :: filexyz, arg_in
  INTEGER :: N_a
  INTEGER :: ip, ist, N_in
  INTEGER :: ix, iy, iz, idx_center
  INTEGER :: iargc  ! pgf90 
  INTEGER :: tstart, counts_per_second, tstop
  CHARACTER(1) :: typ
  REAL(8) :: center(3), xpos, ypos, r_cut
  REAL(8), ALLOCATABLE :: V_short_a(:), Rhoe_a(:)

  CALL system_clock( tstart, counts_per_second )

  Narg = iargc()
  IF( Narg /= 5 ) THEN 
    WRITE(*,*) 'ERROR:'
    WRITE(*,*) 'exactly 5 arguments must be given:'
    WRITE(*,*) 'N, xpos, ypos, r_cut (in bohr), and N_a'
    STOP 
  ENDIF 

  CALL getarg( 1, arg_in )
  READ(arg_in, *) N_in

  CALL getarg( 2, arg_in )
  READ(arg_in, *) xpos

  CALL getarg( 3, arg_in )
  READ(arg_in, *) ypos

  CALL getarg( 4, arg_in )
  READ(arg_in, *) r_cut

  CALL getarg( 5, arg_in )
  READ(arg_in, *) N_a

  ALLOCATE(atpos(3,1))
  atpos(:,1) = (/ xpos, ypos, 8.d0 /)

  typ = 'p'
  IF( typ == 's' ) THEN  ! sinc LF
    NN = (/ N_in, N_in, N_in /)
    hh(:) = (/1.d0, 1.d0, 1.d0/)*(16.d0/(NN(1)-1))
    CALL init_LF3d_sinc( NN, hh )
  ELSEIF(typ == 'p') THEN ! periodic LF
    NN = (/ N_in, N_in, N_in /)
    AA = (/ 0.d0, 0.d0, 0.d0 /)
    BB = (/ 16.d0, 16.d0, 16.d0 /)
    CALL init_LF3d_p( NN, AA, BB )
  ELSE 
    WRITE(*,*)
    WRITE(*,*) 'ERROR: Unknown typ:', typ
    STOP 
  ENDIF 

  CALL info_LF3d()

  center(:) = atpos(:,1)

  CALL init_grid_atom_cube( center, r_cut, N_a, .FALSE. ) 

  idx_center = N_in/2 + 1

  DO ip = 1, Npoints
    ix = lin2xyz(1,ip)
    iy = lin2xyz(2,ip)
    iz = lin2xyz(3,ip)
    IF( iz == idx_center ) THEN 
      WRITE(100,'(1x,3F18.10)') lingrid(1:2,ip), sin(lingrid(1,ip)*PI/8.d0) * cos(lingrid(2,ip)*PI/8.d0)
    ENDIF 
  ENDDO 


  WRITE(101,*) '#', Npoints_a
  DO ip = 1,Npoints_a
    WRITE(101,'(1x,3F18.10)') grid_a(:,ip)
  ENDDO 
 
  CALL dealloc_LF3d()
  CALL dealloc_atoms()

  CALL system_clock( tstop )

  WRITE(*,*)
  WRITE(*,*) 'Total elapsed time: ', dble(tstop - tstart)/counts_per_second, ' second.'
  WRITE(*,*)

END PROGRAM 

#include "interp_Rhoe_a.f90"
#include "interp_Rhoe_a_sinc.f90"

