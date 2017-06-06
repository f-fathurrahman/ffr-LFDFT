#include "gen_gaussian_G.f90"

PROGRAM test_gaussian_G

  USE m_PsPot, ONLY : PsPot_Dir
  USE m_LF3d, ONLY : xyz2lin => LF3d_xyz2lin, &
                     lingrid => LF3d_lingrid, &
                     Npoints => LF3d_Npoints, &
                     LL => LF3d_LL
  USE m_atoms, ONLY : atpos => AtomicCoords
  USE m_constants, ONLY : ANG2BOHR, PI, EPS_SMALL
  !
  IMPLICIT NONE 
  !
  INTEGER :: Narg
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  CHARACTER(64) :: arg_N, arg_alpha
  INTEGER :: ip, N_in, ix, iy, iz
  REAL(8), ALLOCATABLE :: V_gauss(:), V_gauss_ss(:)
  REAL(8) :: dr_vec(3), center(3), dr
  REAL(8) :: alpha
  REAL(8) :: gauss3d_periodic

  Narg = iargc()
  IF( Narg /= 2 ) THEN 
    WRITE(*,*) 'ERROR: exactly two arguments must be given: N and alpha'
    STOP 
  ENDIF 

  CALL getarg( 1, arg_N )
  READ(arg_N, *) N_in

  CALL getarg( 2, arg_alpha )
  READ(arg_alpha, *) alpha

  CALL init_atoms_xyz('ATOM.xyz')
  ! so that coord given in xyz file is in bohr
  atpos(:,:) = atpos(:,:)/ANG2BOHR  

  center(:) = atpos(:,1)

  !
  NN = (/ N_in, N_in, N_in /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)
  CALL init_LF3d_p( NN, AA, BB )

  !
  ALLOCATE( V_gauss(Npoints) )
  WRITE(*,*) 'alpha = ', alpha
  CALL gen_gaussian_G( center, alpha, V_gauss )
  
  !
  ix = 1
  iz = 1
  WRITE(*,*) 'ix iz = ', ix, iz
  DO iy = 1,NN(2)
    ip = xyz2lin(ix,iy,iz)
    WRITE(N_in,'(3F22.12)') lingrid(2,ip), V_gauss(ip), &
                            gauss3d_periodic( LL, center, alpha, lingrid(:,ip) )
  ENDDO 

  ! Free memory
  DEALLOCATE( V_gauss )

  CALL dealloc_LF3d()
  CALL dealloc_atoms()

END PROGRAM


FUNCTION gauss3d_periodic( LL, center, alpha, r ) RESULT(f)
  IMPLICIT NONE 
  !
  REAL(8) :: center(3), LL(3)
  REAL(8) :: alpha
  REAL(8) :: r(3)
  REAL(8) :: dr2, dr_vec(3)
  REAL(8) :: f

  CALL calc_dr_periodic_1pnt( LL, center, r, dr_vec )

  dr2 = dr_vec(1)**2 + dr_vec(2)**2 + dr_vec(3)**2

  f = exp( -alpha*dr2 )

END FUNCTION 
