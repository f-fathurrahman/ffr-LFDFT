SUBROUTINE init_grid_atom( center, cutoff )

  USE m_LF3d, ONLY : lingrid => LF3d_lingrid, &
                     Npoints => LF3d_Npoints
  USE m_LF3d, ONLY : LL => LF3d_LL
  USE m_grid_atom, ONLY : Ngrid_atom, idxa => idx_grid_atom
  IMPLICIT NONE 
  REAL(8) :: center(3)
  REAL(8) :: cutoff
  INTEGER :: ip, ii
  REAL(8) :: x, y, z, dr_vec(3), dr
  !
  INTEGER :: calc_Ngrid_atom
  
  Ngrid_atom = calc_Ngrid_atom(center, cutoff )
  WRITE(*,*) 'Ngrid_atom = ', Ngrid_atom

  ALLOCATE( idxa(Ngrid_atom) )

  ii = 0
  DO ip = 1,Npoints
    x = lingrid(1,ip)
    y = lingrid(2,ip)
    z = lingrid(3,ip)
    !
    CALL calc_dr_periodic_1pnt( LL, center, (/x,y,z/), dr_vec )
    dr = sqrt( dr_vec(1)**2 + dr_vec(2)**2 + dr_vec(3)**2 )
    IF( dr < cutoff ) THEN 
      ii = ii + 1
      idxa(ii) = ip
    ENDIF 
  ENDDO 
END SUBROUTINE 


FUNCTION calc_Ngrid_atom( center, cutoff ) RESULT( Npts )
  USE m_LF3d, ONLY : lingrid => LF3d_lingrid, &
                     Npoints => LF3d_Npoints
  USE m_LF3d, ONLY : LL => LF3d_LL
  IMPLICIT NONE 
  REAL(8) :: center(3)
  REAL(8) :: cutoff
  INTEGER :: Npts, ip
  REAL(8) :: x, y, z, dr_vec(3), dr
  
  Npts = 0
  DO ip = 1,Npoints
    x = lingrid(1,ip)
    y = lingrid(2,ip)
    z = lingrid(3,ip)
    !
    CALL calc_dr_periodic_1pnt( LL, center, (/x,y,z/), dr_vec )
    dr = sqrt( dr_vec(1)**2 + dr_vec(2)**2 + dr_vec(3)**2 )
    IF( dr < cutoff ) THEN 
      Npts = Npts + 1
    ENDIF 
  ENDDO 

END FUNCTION 

