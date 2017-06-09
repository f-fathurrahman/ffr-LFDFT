#include "gen_rho_G.f90"

PROGRAM test_rho_gaussian_G

  USE m_LF3d, ONLY : xyz2lin => LF3d_xyz2lin, &
                     lingrid => LF3d_lingrid, &
                     Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_atoms, ONLY : atpos => AtomicCoords, Natoms, SpeciesSymbols, atm2species
  USE m_constants, ONLY : ANG2BOHR, PI, EPS_SMALL
  USE m_LF3d, ONLY : grid_x => LF3d_grid_x, &
                     grid_y => LF3d_grid_y, &
                     grid_z => LF3d_grid_z
  !
  IMPLICIT NONE 
  !
  INTEGER :: Narg
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  CHARACTER(64) :: arg_N
  INTEGER :: ip, N_in, ix, iy, iz, unitxsf
  REAL(8), ALLOCATABLE :: rho_gauss(:)
  REAL(8) :: center(3), origin(3)
  REAL(8) :: LatVecs(3,3)
  REAL(8) :: integRho
  REAL(8) :: length, znucl, zion, densty

  Narg = iargc()
  IF( Narg /= 1 ) THEN 
    WRITE(*,*) 'ERROR: exactly one argument must be given: N'
    STOP 
  ENDIF 

  CALL getarg( 1, arg_N )
  READ(arg_N, *) N_in

  CALL init_atoms_xyz('ATOM.xyz')
  ! so that coord given in xyz file is in bohr
  atpos(:,:) = atpos(:,:)/ANG2BOHR  

  center(:) = atpos(:,1)

  !
  NN = (/ N_in, N_in, N_in /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)
  CALL init_LF3d_p( NN, AA, BB )

  ! Generate Gaussian density
  ALLOCATE( rho_gauss(Npoints) )
  densty = 0.d0
  zion = 4.d0
  znucl = 6.d0
  CALL atmlength( densty, length, zion, znucl )
  WRITE(*,*) 'length = ', length
  !
  CALL gen_rho_G( center, length, rho_gauss )
  !
  integRho = sum(rho_gauss)*dVol
  WRITE(*,*) 'integRho = ', integRho, sum(rho_gauss)
  ! normalize
  rho_gauss(:) = rho_gauss(:)/integRho*zion
  integRho = sum(rho_gauss)*dVol
  WRITE(*,*) 'renormalized: integRho = ', integRho, sum(rho_gauss)

  LatVecs(:,:) = 0.d0
  LatVecs(1,1) = BB(1) - AA(1)
  LatVecs(2,2) = BB(2) - AA(2)
  LatVecs(3,3) = BB(3) - AA(3)

  unitxsf = 444
  origin(1) = 0.5d0*( grid_x(2) - grid_x(1) )
  origin(2) = 0.5d0*( grid_y(2) - grid_y(1) )
  origin(3) = 0.5d0*( grid_z(2) - grid_z(1) )

  CALL xsf_struct( LatVecs, Natoms, atpos, SpeciesSymbols, atm2species, unitxsf )
  CALL xsf_fast_datagrid_3d(rho_gauss, NN(1), NN(2), NN(3), NN(1), NN(2), NN(3), &
            origin, LatVecs, unitxsf)
  
  !
  ix = 1
  iz = 1
  WRITE(*,*) 'ix iz = ', ix, iz
  DO iy = 1,NN(2)
    ip = xyz2lin(ix,iy,iz)
    WRITE(N_in,'(2F22.12)') lingrid(2,ip), rho_gauss(ip)
  ENDDO 

  ! Free memory
  DEALLOCATE( rho_gauss )

  CALL dealloc_LF3d()
  CALL dealloc_atoms()

END PROGRAM


