#include "init_grid_atom_cube.f90"

PROGRAM test

  USE m_PsPot, ONLY : PsPot_Dir
  USE m_LF3d, ONLY : lingrid => LF3d_lingrid, &
                     Npoints => LF3d_Npoints, &
                     LL => LF3d_LL, &
                     dVol => LF3d_dVol
  USE m_atoms, ONLY : atpos => AtomicCoords
  USE m_constants, ONLY : ANG2BOHR, PI, EPS_SMALL
  USE m_PsPot, ONLY : Ps => Ps_HGH_Params
  USE m_Ps_HGH
  USE m_grid_atom_cube
  !
  IMPLICIT NONE 
  !
  INTEGER :: Narg
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  CHARACTER(64) :: arg_N, arg_N_a
  INTEGER :: ip, N_in, N_a
  REAL(8), ALLOCATABLE :: V_short(:), V_long(:), V_short_a(:)
  REAL(8) :: dr_vec(3), center(3), dr
  INTEGER :: isp
  INTEGER :: iargc

  Narg = iargc()
  IF( Narg /= 2 ) THEN 
    WRITE(*,*) 'ERROR: exactly two arguments must be given: N, N_a'
    STOP 
  ENDIF 

  CALL getarg( 1, arg_N )
  READ(arg_N, *) N_in

  CALL getarg( 2, arg_N_a )
  READ(arg_N_a, *) N_a

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
!  CALL info_LF3d()


  isp = 1  ! use the first species for pseudopotential

  
  CALL init_grid_atom_cube( center, 1.5d0, N_a )

  ALLOCATE( V_short_a(Npoints_a) )
  DO ip = 1,Npoints_a
    CALL calc_dr_periodic_1pnt( LL, center, grid_a(:,ip), dr_vec )
    dr = sqrt( dr_vec(1)**2 + dr_vec(2)**2 + dr_vec(3)**2 )
    V_short_a(ip) = hgh_eval_Vloc_R_short( Ps(isp), dr ) 
!    WRITE(*,*) dr, V_short_a(ip)
  ENDDO 
!  WRITE(*,*) sum(V_short_a)*dVol_a

  ALLOCATE( V_short(Npoints) )
  ALLOCATE( V_long(Npoints) )
  V_short(:) = 0.d0
  DO ip = 1, Npoints
    CALL calc_dr_periodic_1pnt( LL, center, lingrid(:,ip), dr_vec )
    dr = sqrt( dr_vec(1)**2 + dr_vec(2)**2 + dr_vec(3)**2 )
    V_short(ip) = hgh_eval_Vloc_R_short( Ps(isp), dr ) 
    V_long(ip) = hgh_eval_Vloc_R_long( Ps(isp), dr ) 
  ENDDO 
!  WRITE(*,*) 'sum(V_short) = ', sum(V_short)*dVol
!  WRITE(*,*) 'sum(V_long) = ', sum(V_long)*dVol

  WRITE(*,*) sum(V_short(:) + V_long(:))*dVol, sum(V_short_a)*dVol_a + sum(V_long)*dVol

!  WRITE(*,*) dVol_a, dVol


  ! Free memory
  DEALLOCATE( V_short )

  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()
  CALL dealloc_PsPot()
  CALL dealloc_atoms()

END PROGRAM 

