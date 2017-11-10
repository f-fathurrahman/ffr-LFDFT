MODULE m_grid_atom_ss
  IMPLICIT NONE 
  INTEGER :: Ngrid_atom_ss
  INTEGER, ALLOCATABLE :: idx_grid_atom_ss(:)
END MODULE 

SUBROUTINE init_grid_atom_ss( center, cutoff )

  USE m_LF3d_ss, ONLY : lingrid_ss => LF3d_lingrid_ss, &
                                 Npoints_ss => LF3d_Npoints_ss
  USE m_LF3d, ONLY : LL => LF3d_LL
  USE m_grid_atom_ss, ONLY : Ngrid_atom_ss, idx_ss => idx_grid_atom_ss
  IMPLICIT NONE 
  REAL(8) :: center(3)
  REAL(8) :: cutoff
  INTEGER :: ip_ss, ii
  REAL(8) :: x_ss, y_ss, z_ss, dr_vec(3), dr
  !
  INTEGER :: calc_Ngrid_atom_ss
  
  Ngrid_atom_ss = calc_Ngrid_atom_ss( center, cutoff )
  WRITE(*,*) 'Ngrid_ss_atom = ', Ngrid_atom_ss

  ALLOCATE( idx_ss(Ngrid_atom_ss) )

  ii = 0
  DO ip_ss = 1,Npoints_ss
    x_ss = lingrid_ss(1,ip_ss)
    y_ss = lingrid_ss(2,ip_ss)
    z_ss = lingrid_ss(3,ip_ss)
    !
    CALL calc_dr_periodic_1pnt( LL, center, (/x_ss,y_ss,z_ss/), dr_vec )
    dr = sqrt( dr_vec(1)**2 + dr_vec(2)**2 + dr_vec(3)**2 )
    IF( dr < cutoff ) THEN 
      ii = ii + 1
      idx_ss(ii) = ip_ss
    ENDIF 
  ENDDO 
END SUBROUTINE 


FUNCTION calc_Ngrid_atom_ss( center, cutoff ) RESULT( Npts )
  USE m_LF3d_ss, ONLY : lingrid_ss => LF3d_lingrid_ss, &
                                 Npoints_ss => LF3d_Npoints_ss
  USE m_LF3d, ONLY : LL => LF3d_LL
  IMPLICIT NONE 
  REAL(8) :: center(3)
  REAL(8) :: cutoff
  INTEGER :: Npts, ip_ss
  REAL(8) :: x_ss, y_ss, z_ss, dr_vec(3), dr
  
  !! FIXME: This is only for periodic case
  !!

  Npts = 0
  DO ip_ss = 1,Npoints_ss
    x_ss = lingrid_ss(1,ip_ss)
    y_ss = lingrid_ss(2,ip_ss)
    z_ss = lingrid_ss(3,ip_ss)
    !
    CALL calc_dr_periodic_1pnt( LL, center, (/x_ss,y_ss,z_ss/), dr_vec )
    dr = sqrt( dr_vec(1)**2 + dr_vec(2)**2 + dr_vec(3)**2 )
    IF( dr < cutoff ) THEN 
      Npts = Npts + 1
    ENDIF 
  ENDDO 

END FUNCTION 

