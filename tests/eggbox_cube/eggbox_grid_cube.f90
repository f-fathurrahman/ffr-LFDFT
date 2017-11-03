PROGRAM eggbox_grid_cube

  USE m_constants, ONLY : Ry2eV
  USE m_options, ONLY : FREE_NABLA2
  USE m_PsPot, ONLY : NbetaNL, PsPot_FilePath
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  !
  USE m_atoms, ONLY : atpos => AtomicCoords, atm2species, Zv => AtomicValences, &
                      Natoms, Nspecies
  USE m_hamiltonian, ONLY : V_ps_loc, Rhoe, V_ps_loc_long
  USE m_energies, ONLY : E_ps_loc
  USE m_LF3d, ONLY : dVol => LF3d_dVol
  USE m_options, ONLY : T_PRINT_INTEG_RHO, I_ALG_DIAG, ETHR_EVALS
  USE m_constants, ONLY : ANG2BOHR
  USE m_grid_atom_cube, ONLY : Npoints_a, dVol_a
  !
  USE m_atoms, ONLY : SpeciesSymbols
  USE m_PsPot, ONLY : Ps_HGH_Params, NbetaNL
  USE m_Ps_HGH, ONLY : init_Ps_HGH_Params
  !
  IMPLICIT NONE 
  INTEGER :: Narg
  INTEGER :: NN(3)
  INTEGER :: N_a
  REAL(8) :: hh(3), AA(3), BB(3)
  CHARACTER(64) :: arg_in
  INTEGER :: ip, ist, N_in
  INTEGER :: iargc  ! pgf90 
  INTEGER :: tstart, counts_per_second, tstop
  CHARACTER(1) :: typ
  REAL(8) :: center(3), xpos, ypos, r_cut
  REAL(8), ALLOCATABLE :: V_short_a(:), Rhoe_a(:)

  CALL system_clock( tstart, counts_per_second )

  Narg = iargc()
  IF( Narg /= 6 ) THEN 
    WRITE(*,*) 'ERROR:'
    WRITE(*,*) 'exactly r arguments must be given:'
    WRITE(*,*) 'N, path to PsPot file, xpos, ypos, r_cut, N_a'
    STOP 
  ENDIF 

  CALL getarg( 1, arg_in )
  READ(arg_in, *) N_in

  CALL getarg( 2, arg_in )
  ALLOCATE( PsPot_FilePath(1) )
  PsPot_FilePath(1) = arg_in
  
  CALL getarg( 3, arg_in )
  READ(arg_in, *) xpos
  
  CALL getarg( 4, arg_in )
  READ(arg_in, *) ypos

  CALL getarg( 5, arg_in )
  READ(arg_in, *) r_cut

  CALL getarg( 6, arg_in )
  READ(arg_in, *) N_a

  ALLOCATE( atpos(3,1) )

  atpos(:,1) = (/ xpos, ypos, 8.d0 /)

  ALLOCATE( atm2species(1) )
  Natoms = 1
  Nspecies = 1
  atm2species(1) = 1

  !CALL init_PsPot()
  ALLOCATE( SpeciesSymbols(1) )
  ALLOCATE( Ps_HGH_Params(1) )
  !
  CALL init_Ps_HGH_Params( Ps_HGH_Params(1), PsPot_FilePath(1) )

  ALLOCATE( Zv(1) )
  Zv = Ps_HGH_Params(1)%Zval

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

  CALL info_PsPot()
  CALL info_LF3d()

  ! Only use local pseudopotential
  NbetaNL = 0

  ! Initialize occupation numbers
  CALL init_states()

  IF( typ == 'p' ) THEN 
    CALL init_strfact_shifted()
  ENDIF 

  ! Memory for potentials
  CALL alloc_hamiltonian()

  ! Local pseudopotential
  IF( typ == 'p' ) THEN 
    CALL init_V_ps_loc_G()
  ELSEIF( typ == 's' ) THEN  
    CALL init_V_ps_loc()
  ELSE 
    WRITE(*,*) 'ERROR: unknown typ = ', typ
    STOP 
  ENDIF 

  ! Local pseudopotential (long part)
  ! FIXME: only for periodic CASE
  ALLOCATE( V_ps_loc_long(Npoints) )
  IF( typ=='s') THEN 
    CALL init_V_ps_loc_long()
  ELSEIF( typ=='p') THEN 
    CALL init_V_ps_loc_G_long()
  ELSE 
    WRITE(*,*) 'ERROR: Unknown typ:', typ
    STOP 
  ENDIF 

  !
  center(:) = atpos(:,1)
  CALL init_grid_atom_cube( center, r_cut, N_a, .FALSE. ) 
  !
  ALLOCATE( V_short_a(Npoints_a) )
  CALL init_V_ps_loc_short( center, V_short_a, typ )
 
  ! Laplacian matrix
  CALL init_nabla2_sparse()
  ! ILU0 preconditioner based on kinetic matrix
  CALL init_ilu0_prec()

  IF( FREE_NABLA2 ) THEN 
    CALL dealloc_nabla2_sparse()
  ENDIF 

  ! Manually allocate KS eigenvectors and eigenvalues
  ALLOCATE( evecs(Npoints,Nstates), evals(Nstates) )

  ! Initialize to random wavefunction
  DO ist = 1, Nstates
    DO ip = 1, Npoints
      CALL random_number( evecs(ip,ist) )
    ENDDO
  ENDDO
  CALL orthonormalize( Nstates, evecs )

  
  ! Diagonalize
  I_ALG_DIAG = 3
  ETHR_EVALS = 1d-5;
  CALL Sch_solve_diag()

  ! Calculate electron density
  T_PRINT_INTEG_RHO = .TRUE.
  CALL calc_rhoe( Focc, evecs )


  ! Calculate local pseudopotential energy
  E_ps_loc = sum( Rhoe(:) * V_ps_loc(:) )*dVol

  !
  ALLOCATE( Rhoe_a(Npoints_a) )
  !
  IF( typ == 'p' ) THEN 
    ! periodic
    CALL interp_Rhoe_a( Rhoe, Rhoe_a )
  ELSE
    ! non periodic
    CALL interp_Rhoe_a_sinc( Rhoe, Rhoe_a )
  ENDIF 

  WRITE(*,*)
  WRITE(*,*) 'These two values should be close:'
  WRITE(*,*) 'integ(V_short a) = ', sum(V_short_a)*dVol_a
  WRITE(*,*) 'integ(V_short t) = ', sum(V_ps_loc(:) - V_ps_loc_long(:))*dVol

  WRITE(*,*)
  WRITE(*,*) 'These two values should be close:'
  WRITE(*,'(1x,A,F18.10)') 'Ps short a:', sum(V_short_a(:)*Rhoe_a(:))*dVol_a
  WRITE(*,'(1x,A,F18.10)') 'Ps short t:', sum( (V_ps_loc(:)-V_ps_loc_long(:)) *Rhoe(:))*dVol

  WRITE(*,*)
  WRITE(*,'(1x,A,F18.10)') 'Ps long  a/t:', sum(V_ps_loc_long(:)*Rhoe(:))*dVol

  WRITE(*,*)
  WRITE(*,*) 'These two values should be close:'
  WRITE(*,'(1x,A,F18.10)') 'E_ps_loc = ', E_ps_loc
  WRITE(*,'(1x,A,F18.10)') 'Ps total a:', sum(V_short_a(:)*Rhoe_a(:))*dVol_a + sum(V_ps_loc_long(:)*Rhoe(:))*dVol

  CALL dump_data()


  CALL dealloc_nabla2_sparse()
  CALL dealloc_ilu0_prec()
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()
  CALL dealloc_PsPot()
  CALL dealloc_atoms()

  CALL system_clock( tstop )

  WRITE(*,*)
  WRITE(*,*) 'Total elapsed time: ', dble(tstop - tstart)/counts_per_second, ' second.'
  WRITE(*,*)

CONTAINS 

SUBROUTINE dump_data()
  USE m_LF3d, ONLY : grid_x => LF3d_grid_x, &
                     lingrid => LF3d_lingrid, &
                     lin2xyz => LF3d_lin2xyz
  USE m_grid_atom_cube
  INTEGER :: idx_center
  INTEGER :: iz

  idx_center = N_in/2 + 1
  WRITE(*,*) 'idx_center = ', idx_center
  WRITE(*,*) 'x = ', grid_x(idx_center)

  DO ip = 1, Npoints
    iz = lin2xyz(3,ip)
    IF( iz == idx_center ) THEN 
      WRITE(100,'(1x,3F18.10)') lingrid(1:2,ip), Rhoe(ip)
    ENDIF 
  ENDDO 

  idx_center = N_a/2 + 1
  DO ip = 1, Npoints_a
    iz = lin2xyz_a(3,ip)
    IF( iz == idx_center ) THEN 
      WRITE(101,'(1x,4F18.10)') lingrid_a(1:2,ip), Rhoe_a(ip), V_short_a(ip)
    ENDIF 
  ENDDO 

END SUBROUTINE 



END PROGRAM 

SUBROUTINE init_V_ps_loc_short( center, V_short_a, typ )
  USE m_grid_atom_cube, ONLY : Npoints_a, lingrid_a
  USE m_LF3d, ONLY : LL => LF3d_LL
  USE m_PsPot, ONLY : Ps => Ps_HGH_Params
  USE m_Ps_HGH, ONLY : hgh_eval_Vloc_R_short
  IMPLICIT NONE 
  REAL(8) :: center(3)
  REAL(8) :: V_short_a(Npoints_a)
  CHARACTER(1) :: typ
  INTEGER :: ip, isp
  REAL(8) :: dr_vec(3)
  REAL(8) :: dr

  isp = 1
  IF( typ == 'p' ) THEN 
    ! periodic case
    DO ip = 1,Npoints_a
      CALL calc_dr_periodic_1pnt( LL, center, lingrid_a(:,ip), dr_vec )
      dr = sqrt( dr_vec(1)**2 + dr_vec(2)**2 + dr_vec(3)**2 )
      V_short_a(ip) = hgh_eval_Vloc_R_short( Ps(isp), dr ) 
    ENDDO 
  ELSE 
    ! non-periodic
    DO ip = 1,Npoints_a
      CALL calc_dr_1pnt( center, lingrid_a(:,ip), dr )
      V_short_a(ip) = hgh_eval_Vloc_R_short( Ps(isp), dr ) 
    ENDDO 
  ENDIF 
  WRITE(*,*) 'sum(V_short_a) = ', sum(V_short_a)

END SUBROUTINE 

#include "interp_Rhoe_a.f90"
#include "interp_Rhoe_a_sinc.f90"

