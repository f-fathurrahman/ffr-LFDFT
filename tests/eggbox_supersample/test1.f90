#include "m_supersample.f90"
#include "m_LF3d_ss.f90"
#include "m_grid_atom.f90"
#include "m_grid_atom_ss.f90"
#include "do_supersample.f90"
#include "sinc.f90"

PROGRAM test
  USE m_LF3d, ONLY : grid_x => LF3d_grid_x, &
                     lingrid => LF3d_lingrid, &
                     lin2xyz => LF3d_lin2xyz, &
                     Npoints => LF3d_Npoints
  USE m_grid_atom
  USE m_LF3d_ss, ONLY : lingrid_ss => LF3d_lingrid_ss, &
                        lin2xyz_ss => LF3d_lin2xyz_ss, &
                        Npoints_ss => LF3d_Npoints_ss
  USE m_grid_atom_ss
  IMPLICIT NONE 
  REAL(8) :: center(3)
  REAL(8) :: c1, c2
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  REAL(8) :: cutoff
  !
  INTEGER :: Nargs, iargc
  CHARACTER(32) :: arg_in
  INTEGER :: N_in, idx_center, idx_center_ss
  INTEGER :: ip, ip_a, ix, iy, iz

  Nargs = iargc()
  IF( Nargs /= 4 ) THEN 
    WRITE(*,*) 'ERROR: need 4 arguments'
    STOP 
  ENDIF 

  !
  ! Parse arguments
  !
  CALL getarg( 1, arg_in )
  READ(arg_in,*) N_in

  CALL getarg( 2, arg_in )
  READ(arg_in,*) c1
  
  CALL getarg( 3, arg_in )
  READ(arg_in,*) c2

  CALL getarg( 4, arg_in )
  READ(arg_in,*) cutoff

  WRITE(*,*) 'N_in   = ', N_in
  WRITE(*,*) 'c1     = ', c1
  WRITE(*,*) 'c2     = ', c2
  WRITE(*,*) 'cutoff = ', cutoff

  WRITE(99,'(1x,3F18.10)') c1, c2, cutoff
  
  center(:) = (/ c1, c2, 8.d0 /)

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


CONTAINS 

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

