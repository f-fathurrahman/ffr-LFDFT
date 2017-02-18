! ffr: January 2016

SUBROUTINE init_Vpsloc()
  USE m_globals, ONLY : Vpsloc, ATOMS, PSPOTS, Npoints, LF
  USE m_ps_hgh, ONLY : hgh_eval_Vloc
  IMPLICIT NONE
  REAL(8) :: r, r0(3)
  INTEGER :: ia, ip
 
  WRITE(*,'(/,1x,A)') 'Initializing local ionic pseudopotential'
  WRITE(*,*)          '----------------------------------------'

  ! Vpsloc should be already allocated before
  Vpsloc(:) = 0.d0
  
  DO ia = 1, ATOMS%Natoms
    r0 = ATOMS%positions(:,ia)
    DO ip = 1, Npoints
      r = norm2( LF%lingrid(:,ip)-r0(:) )
      Vpsloc(ip) = Vpsloc(ip) + hgh_eval_Vloc( PSPOTS( ATOMS%atmToSpecies(ia) ), r )
    ENDDO
  ENDDO
  ! This is sometimes useful to detect NaN due to 1/r term in the pseudopotential
  WRITE(*,*) 'sum(Vpsloc) = ', sum(Vpsloc)

END SUBROUTINE


SUBROUTINE init_Vpsloc_per()
  USE m_globals, ONLY : Vpsloc, ATOMS, PSPOTS, Npoints, LF
  USE m_ps_hgh, ONLY : hgh_eval_Vloc
  IMPLICIT NONE
  REAL(8) :: r, r0(3)
  INTEGER :: ia, ip
  REAL(8) :: Lx, Ly, Lz
  REAL(8) :: xx,xx1,xx2,xx3, yy,yy1,yy2,yy3, zz,zz1,zz2,zz3
 
  WRITE(*,'(/,1x,A)') 'Initializing local ionic pseudopotential (periodic)'
  WRITE(*,*)          '---------------------------------------------------'

  ! Vpsloc should be already allocated before
  Vpsloc(:) = 0.d0
  
  Lx = LF%Lx
  Ly = LF%Ly
  Lz = LF%Lz

  DO ia = 1, ATOMS%Natoms
    r0 = ATOMS%positions(:,ia)
    DO ip = 1, Npoints
      ! FIXME: only works for orthorombic cells
      xx1 = abs( LF%lingrid(1,ip) - r0(1) )
      xx2 = abs( r0(1) + Lx - LF%lingrid(1,ip) )
      xx3 = abs( LF%lingrid(1,ip) - r0(1) - Lx)
      xx  = minval( (/xx1,xx2,xx3/) )
      !
      yy1 = abs( LF%lingrid(2,ip) - r0(2) )
      yy2 = abs( r0(2) + Ly - LF%lingrid(2,ip) )
      yy3 = abs( LF%lingrid(2,ip) - r0(2) - Ly)
      yy  = minval( (/yy1,yy2,yy3/) )
      !
      zz1 = abs( LF%lingrid(3,ip) - r0(3) )
      zz2 = abs( r0(3) + Lz - LF%lingrid(3,ip) )
      zz3 = abs( LF%lingrid(3,ip) - r0(3) - Lz)
      zz  = minval( (/zz1,zz2,zz3/) )
      !
      r = sqrt( xx**2 + yy**2 + zz**2 )
      !
      Vpsloc(ip) = Vpsloc(ip) + hgh_eval_Vloc( PSPOTS( ATOMS%atmToSpecies(ia) ), r )
    ENDDO
  ENDDO
  ! This is sometimes useful to detect NaN due to 1/r term in the pseudopotential
  WRITE(*,*) 'sum(Vpsloc) = ', sum(Vpsloc)

END SUBROUTINE


