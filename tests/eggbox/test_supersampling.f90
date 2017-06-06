#include "m_LF3d_supersample.f90"
#include "init_LF3d_supersample.f90"
#include "deallocate_LF3d_supersample.f90"
#include "sinc.f90"
#include "m_grid_ss_atom.f90"
#include "init_grid_ss_atom.f90"
#include "m_grid_atom.f90"
#include "init_grid_atom.f90"
#include "supersample.f90"

PROGRAM test_eggbox

  USE m_PsPot, ONLY : PsPot_Dir
  USE m_LF3d, ONLY : xyz2lin => LF3d_xyz2lin, &
                     lingrid => LF3d_lingrid, &
                     Npoints => LF3d_Npoints, &
                     LL => LF3d_LL
  USE m_atoms, ONLY : atpos => AtomicCoords
  USE m_constants, ONLY : ANG2BOHR, PI, EPS_SMALL
  USE m_PsPot, ONLY : Ps => Ps_HGH_Params
  USE m_Ps_HGH, ONLY : hgh_eval_Vloc_R_short
  !
  IMPLICIT NONE 
  !
  INTEGER :: Narg
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  CHARACTER(64) :: arg_N
  INTEGER :: ip, N_in, ix, iy, iz
  REAL(8), ALLOCATABLE :: V_short(:), V_short_ss(:)
  REAL(8) :: dr_vec(3), center(3), dr
  INTEGER :: l, m, isp, iprj
  REAL(8) :: Ylm_real

  Narg = iargc()
  IF( Narg /= 1 ) THEN 
    WRITE(*,*) 'ERROR: exactly one arguments must be given: N'
    STOP 
  ENDIF 

  CALL getarg( 1, arg_N )
  READ(arg_N, *) N_in

  CALL init_atoms_xyz('ATOM.xyz')
  ! so that coord given in xyz file is in bohr
  atpos(:,:) = atpos(:,:)/ANG2BOHR  

  center(:) = atpos(:,1)

  ! Override PsPot_Dir
  PsPot_Dir = '../../HGH/'
  CALL init_PsPot()

  !
  NN = (/ N_in, N_in, N_in /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)
  CALL init_LF3d_p( NN, AA, BB )

  CALL init_LF3d_supersample( 3 )

  CALL init_strfact_shifted()

  CALL alloc_hamiltonian()

  !
  l = 0
  m = 0
  isp = 1
  !
  ALLOCATE( V_short(Npoints) )
  V_short(:) = 0.d0
  DO ip = 1, Npoints
    CALL calc_dr_periodic_1pnt( LL, center, lingrid(:,ip), dr_vec )
    dr = sqrt( dr_vec(1)**2 + dr_vec(2)**2 + dr_vec(3)**2 )
    V_short(ip) = hgh_eval_Vloc_R_short( Ps(isp), dr ) 
  ENDDO 
  WRITE(*,*) 'sum(V_short) = ', sum(V_short)
  
  !
  ALLOCATE( V_short_ss(Npoints) )
  
  CALL init_grid_atom( center, 2.5d0 )
  CALL init_grid_ss_atom( center, 2.5d0 )
  !
  CALL supersample( V_short, V_short_ss )

  ix = 1
  iz = 1
  WRITE(*,*) 'ix iz = ', ix, iz
  DO iy = 1,NN(2)
    ip = xyz2lin(ix,iy,iz)
    WRITE(N_in,'(4F22.12)') lingrid(2,ip), V_short(ip), V_short_ss(ip), 
  ENDDO 

  ! Free memory
  DEALLOCATE( V_short, V_short_ss )

  CALL dealloc_LF3d_supersample()
  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()
  CALL dealloc_PsPot()
  CALL dealloc_atoms()

END PROGRAM

