#include "m_supersample.f90"
#include "m_LF3d_ss.f90"
#include "m_grid_atom.f90"
#include "m_grid_atom_ss.f90"
#include "init_V_ps_loc_short_ss.f90"
#include "sinc.f90"

PROGRAM test
  USE m_LF3d, ONLY : grid_x => LF3d_grid_x, &
                     lingrid => LF3d_lingrid, &
                     lin2xyz => LF3d_lin2xyz, &
                     Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_grid_atom
  USE m_LF3d_ss, ONLY : lingrid_ss => LF3d_lingrid_ss, &
                        lin2xyz_ss => LF3d_lin2xyz_ss, &
                        Npoints_ss => LF3d_Npoints_ss
  USE m_grid_atom_ss
  USE m_PsPot, ONLY : NbetaNL, PsPot_FilePath, Ps_HGH_Params
  USE m_atoms, ONLY : atpos => AtomicCoords, atm2species, Zv => AtomicValences, &
                      Natoms, Nspecies, SpeciesSymbols
  USE m_Ps_HGH, ONLY : init_Ps_HGH_Params
  USE m_options, ONLY : FREE_NABLA2, T_PRINT_INTEG_RHO, I_ALG_DIAG, ETHR_EVALS
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  USE m_hamiltonian, ONLY : Rhoe, V_ps_loc
  USE m_energies, ONLY : E_ps_loc

  IMPLICIT NONE 
  REAL(8) :: center(3)
  REAL(8) :: xpos, ypos
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  REAL(8) :: cutoff
  !
  INTEGER :: Nargs
  CHARACTER(32) :: arg_in
  INTEGER :: N_in, idx_center, idx_center_ss
  INTEGER :: ip, ip_a, ix, iy, iz
  INTEGER :: ist
  !
  REAL(8), ALLOCATABLE :: V_ps_loc_short_ss(:)

  CALL parse_arguments()
  
  CALL setup_dummy_atoms()

  NN = (/ N_in, N_in, N_in /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)
  CALL init_LF3d_p( NN, AA, BB )
  CALL info_LF3d()

  CALL init_LF3d_p_ss( 3 )
  CALL init_grid_atom( center, cutoff )
  CALL init_grid_atom_ss( center, cutoff )

  CALL output_grid()
  CALL output_grid_ss()


  CALL setup_potentials()


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

  WRITE(*,*)
  WRITE(*,*) 'E_ps_loc = ', E_ps_loc


CONTAINS 


SUBROUTINE setup_potentials()
  USE m_hamiltonian, ONLY : V_ps_loc_long, V_ps_loc

  LOGICAL :: T_DO_SS

  CALL init_strfact_shifted()

  ! Memory for potentials
  CALL alloc_hamiltonian()

  T_DO_SS = .FALSE.
  IF( T_DO_SS ) THEN 
  ! Do supersampling stuffs
    ALLOCATE( V_ps_loc_long(Npoints) )
    ! Local pseudopotential (long part)
    CALL init_V_ps_loc_G_long()
    ALLOCATE( V_ps_loc_short_ss(Npoints) )
    CALL init_V_ps_loc_short_ss( V_ps_loc_short_ss )
    V_ps_loc(:) = V_ps_loc_long(:) + V_ps_loc_short_ss(:)

  ELSE 
  ! no supersampling
    WRITE(*,*)
    WRITE(*,*) 'No supersampling'
    CALL init_V_ps_loc_G()
  ENDIF 
  
  ! Only use local pseudopotential
  NbetaNL = 0

  ! Initialize occupation numbers
  CALL init_states()

END SUBROUTINE 


SUBROUTINE setup_dummy_atoms()

  ALLOCATE( atpos(3,1) )
  atpos(:,1) = (/ xpos, ypos, 8.d0 /)

  ALLOCATE( atm2species(1) )
  Natoms = 1
  Nspecies = 1
  atm2species(1) = 1

  ALLOCATE( SpeciesSymbols(1) )
  ALLOCATE( Ps_HGH_Params(1) )
  !
  CALL init_Ps_HGH_Params( Ps_HGH_Params(1), PsPot_FilePath(1) )

  ALLOCATE( Zv(1) )
  Zv = Ps_HGH_Params(1)%Zval

  center(:) = atpos(:,1)
END SUBROUTINE 



SUBROUTINE parse_arguments()
  INTEGER :: iargc
  !
  Nargs = iargc()
  IF( Nargs /= 5 ) THEN 
    WRITE(*,*) 'ERROR: need 5 arguments'
    STOP 
  ENDIF 

  !
  ! Parse arguments
  !
  CALL getarg( 1, arg_in )
  READ(arg_in,*) N_in

  CALL getarg( 2, arg_in )
  ALLOCATE( PsPot_FilePath(1) )
  PsPot_FilePath(1) = arg_in

  CALL getarg( 3, arg_in )
  READ(arg_in,*) xpos
  
  CALL getarg( 4, arg_in )
  READ(arg_in,*) ypos

  CALL getarg( 5, arg_in )
  READ(arg_in,*) cutoff

  WRITE(*,*)
  WRITE(*,*) 'N_in   = ', N_in
  WRITE(*,*) 'xpos   = ', xpos
  WRITE(*,*) 'ypos   = ', ypos
  WRITE(*,*) 'cutoff = ', cutoff

  WRITE(99,'(1x,3F18.10)') xpos, ypos, cutoff

END SUBROUTINE 


SUBROUTINE output_grid()

  idx_center = N_in/2 + 1
  WRITE(*,*) 'idx_center = ', idx_center
  WRITE(*,*) 'x = ', grid_x(idx_center)

  DO ip = 1, Npoints
    ix = lin2xyz(1,ip)
    iy = lin2xyz(2,ip)
    iz = lin2xyz(3,ip)
    IF( iz == idx_center ) THEN 
      WRITE(100,'(1x,3F18.10)') lingrid(1:3,ip)
    ENDIF 
  ENDDO 

  DO ip_a = 1, Ngrid_atom
    ip = idx_grid_atom(ip_a)
    ix = lin2xyz(1,ip)
    iy = lin2xyz(2,ip)
    iz = lin2xyz(3,ip)
    IF( iz == idx_center ) THEN 
      WRITE(101,'(1x,3F18.10)') lingrid(1:3,ip)
    ENDIF 
  ENDDO 

END SUBROUTINE 


SUBROUTINE output_grid_ss()

  idx_center_ss = 3*N_in/2 + 1
  WRITE(*,*) 'idx_center_ss = ', idx_center_ss

  DO ip = 1, Npoints_ss
    ix = lin2xyz_ss(1,ip)
    iy = lin2xyz_ss(2,ip)
    iz = lin2xyz_ss(3,ip)
    IF( iz == idx_center_ss ) THEN 
      WRITE(102,'(1x,3F18.10)') lingrid_ss(1:3,ip)
    ENDIF 
  ENDDO 

  DO ip_a = 1, Ngrid_atom_ss
    ip = idx_grid_atom_ss(ip_a)
    ix = lin2xyz_ss(1,ip)
    iy = lin2xyz_ss(2,ip)
    iz = lin2xyz_ss(3,ip)
    IF( iz == idx_center_ss ) THEN 
      WRITE(103,'(1x,3F18.10)') lingrid_ss(1:3,ip)
    ENDIF 
  ENDDO 

END SUBROUTINE 


END PROGRAM 

