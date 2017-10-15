PROGRAM test_interp_Rhoe

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

  CALL info_atoms()
  CALL info_PsPot()
  CALL info_LF3d()

! Only use local pseudopotential
  ! CALL init_betaNL()
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
  ELSE 
    CALL init_V_ps_loc()
  ENDIF 

  !
  center(:) = atpos(:,1)
  CALL init_grid_atom_cube( center, 1.0d0, 55 )  ! 1.0 is quite reasonable for hydrogen atom
  !
  ALLOCATE( V_short_a(Npoints_a) )
  CALL init_V_ps_loc_short( center, V_short_a )

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
  WRITE(*,'(1x,A,F18.10)') 'E_ps_loc = ', E_ps_loc

  !
  ALLOCATE( Rhoe_a(Npoints) )
  CALL interp_Rhoe( Rhoe, Rhoe_a )

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

END PROGRAM 


SUBROUTINE interp_Rhoe( Rhoe, Rhoe_a )
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     NN => LF3d_NN, &
                     lingrid => LF3d_lingrid, &
                     xyz2lin => LF3d_xyz2lin, &
                     grid_x => LF3d_grid_x, &
                     grid_y => LF3d_grid_y, &
                     grid_z => LF3d_grid_z, &
                     dVol => LF3d_dVol
  USE bspline
  IMPLICIT NONE 
  REAL(8) :: Rhoe(Npoints)
  REAL(8) :: Rhoe_a(Npoints)
  !
  INTEGER :: Nx, Ny, Nz, kx, ky, kz
  INTEGER :: iknot
  REAL(8), ALLOCATABLE :: x(:), y(:), z(:)
  REAL(8), ALLOCATABLE :: tx(:), ty(:), tz(:)
  REAL(8), ALLOCATABLE :: Rhoe_tmp(:,:,:), bcoef(:,:,:)
  INTEGER :: iflag
  INTEGER :: ip
  !
  INTEGER :: i,j,k, ii,jj,kk, ip_a
  REAL(8) :: shiftx, shifty, shiftz
  INTEGER :: idx, idy, idz, iloy, iloz, inbvx, inbvy, inbvz
  REAL(8) :: val, dx, dy, dz

  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  iknot = 0

  kx = 4
  ky = 4
  kz = 4
  ALLOCATE( tx(Nx+1+kz), ty(Ny+1+ky), tz(Nz+1+kz) )
  ALLOCATE( x(Nx+1), y(Ny+1), z(Nz+1) )

  shiftx = 0.5d0*( grid_x(2) - grid_x(1) )
  !shiftx = 0.d0
  DO i = 1,Nx
    x(i) = grid_x(i) - shiftx
  ENDDO 
  x(Nx+1) = x(Nx) + 2*shiftx

  shifty = 0.5d0*( grid_y(2) - grid_y(1) )
  !shifty = 0.d0
  DO j = 1,Ny
    y(j) = grid_y(j) - shifty
  ENDDO 
  y(Ny+1) = y(Ny) + 2*shifty

  shiftz = 0.5d0*( grid_z(2) - grid_z(1) )
  !shiftz = 0.d0
  DO k = 1,Nz
    z(k) = grid_z(k) - shiftz
  ENDDO 
  z(Nz+1) = z(Nz) + 2*shiftz

  Rhoe_a(:) = 0.d0

  ALLOCATE( Rhoe_tmp(Nx+1,Ny+1,Nz+1) )
  ALLOCATE( bcoef(Nx+1,Ny+1,Nz+1) )

  iloy = 1
  iloz = 1
  inbvx = 1
  inbvy = 1
  inbvz = 1
  idx = 0
  idy = 0
  idz = 0

  ! copy to interp_Rhoe
  DO kk = 1,Nz+1
  DO jj = 1,Ny+1
  DO ii = 1,Nx+1
    i = ii
    j = jj
    k = kk
    IF( kk == Nz+1 ) k = 1
    IF( jj == Ny+1 ) j = 1
    IF( ii == Nx+1 ) i = 1
    ip = xyz2lin(i,j,k)
    !
    Rhoe_tmp(ii,jj,kk) = Rhoe(ip)
  ENDDO 
  ENDDO 
  ENDDO 

  CALL db3ink( x, Nx+1, y, Ny+1, z, Nz+1, Rhoe_tmp, kx,ky,kz, iknot, tx,ty,tz, bcoef, iflag )
  WRITE(*,*) 'db3ink iflag = ', iflag

  DO ip = 1, Npoints
    dx = lingrid(1,ip) - shiftx
    dy = lingrid(2,ip) - shifty
    dz = lingrid(3,ip) - shiftz
    CALL db3val( dx, dy, dz, idx,idy,idz, tx,ty,tz, Nx+1,Ny+1,Nz+1,kx,ky,kz, bcoef,&
                 val, iflag, inbvx, inbvy, inbvz, iloy, iloz )
    IF( iflag /= 0 ) THEN 
      WRITE(*,*) 'ERROR in calling db3val: iflag = ', iflag
      STOP 
    ENDIF 
    Rhoe_a(ip_a) = val
  ENDDO 

  WRITE(*,*) 'integ(Rhoe_a) = ', sum(Rhoe_a)*dVol
  WRITE(*,*) 'integ(Rhoe) = ', sum(Rhoe)*dVol

  DEALLOCATE( x, y, z )
  DEALLOCATE( tx, ty, tz )
  DEALLOCATE( Rhoe_tmp )
  DEALLOCATE( bcoef )

END SUBROUTINE 

SUBROUTINE init_V_ps_loc_short( center, V_short_a )
  USE m_grid_atom_cube, ONLY : Npoints_a, grid_a
  USE m_LF3d, ONLY : LL => LF3d_LL
  USE m_PsPot, ONLY : Ps => Ps_HGH_Params
  USE m_Ps_HGH, ONLY : hgh_eval_Vloc_R_short
  IMPLICIT NONE 
  REAL(8) :: center(3)
  REAL(8) :: V_short_a(Npoints_a)
  INTEGER :: ip, isp
  REAL(8) :: dr_vec(3)
  REAL(8) :: dr

  isp = 1
  DO ip = 1,Npoints_a
    CALL calc_dr_periodic_1pnt( LL, center, grid_a(:,ip), dr_vec )
    dr = sqrt( dr_vec(1)**2 + dr_vec(2)**2 + dr_vec(3)**2 )
    V_short_a(ip) = hgh_eval_Vloc_R_short( Ps(isp), dr ) 
  ENDDO 
  WRITE(*,*) 'sum(V_short_a) = ', sum(V_short_a)

END SUBROUTINE 

#include "interp_Rhoe_a.f90"

