! efefer, 15 March 2016

! Various subroutines for generating initial (guess) rho


!------------------------------------------------------------------------------
SUBROUTINE init_rho_gaussian_per()
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : Rho, N, Nstate, LF, ATOMS
  IMPLICIT NONE
  INTEGER :: ip, ia, is
  REAL(8) :: r, sigma, r0(3), Zval
  REAL(8) :: deltaV, normChg
  REAL(8) :: integRho
  REAL(8) :: Lx,Ly,Lz
  REAL(8) :: xx1,xx2,xx3,xx, yy1,yy2,yy3,yy, zz1,zz2,zz3,zz

  Lx = LF%Lx
  Ly = LF%Ly
  Lz = LF%Lz

  WRITE(*,*) 'Lx,Ly,Lz = ', Lx,Ly,Lz

  deltaV = LF%LFx%h * LF%LFy%h * LF%LFz%h
  sigma = 0.75d0  ! XXX HARDCODED
  WRITE(*,'(/1x,A)') 'Generating initial Rho (gaussian)'
  WRITE(*,*)         '---------------------------------'
  normChg = 0.d0
  DO ia = 1, ATOMS%Natoms
    r0 = ATOMS%positions(:,ia)
    is = ATOMS%atmToSpecies(ia)
    Zval = ATOMS%Zv(is)
    normChg = normChg + Zval
    !WRITE(*,*) ia, ZZ  ! DEBUG
    DO ip = 1, N**3
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
      WRITE(*,'(A,4F10.5)') 'x, r = ', LF%lingrid(:,ip), r
      ! XXX For the moment only s-like rho
      Rho(ip) = Rho(ip) + ZZ*exp(-r**2/(2.d0*sigma**2))/(2.d0*pi*sigma**2)**1.5d0
    ENDDO
  ENDDO
  integRho = sum(Rho(:))*deltaV
  WRITE(*,*) 'Integrated Rho:', integRho 
  !WRITE(*,*) 'normChg:', normChg   ! DEBUG
  ! Renormalize
  Rho(:) = Rho(:)*normChg/integRho
  WRITE(*,*) 'Renormalized to:', sum(Rho)*deltaV

END SUBROUTINE


!------------------------------------------------------------------------------
SUBROUTINE init_rho_gaussian()
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : Rho, N, Nstate, LF, ATOMS
  IMPLICIT NONE
  INTEGER :: ip, ia, is
  REAL(8) :: r, sigma, r0(3), ZZ
  REAL(8) :: deltaV, normChg
  REAL(8) :: integRho

  deltaV = LF%LFx%h * LF%LFy%h * LF%LFz%h
  sigma = 0.75d0  ! XXX HARDCODED
  WRITE(*,'(/1x,A)') 'Generating initial Rho (gaussian)'
  WRITE(*,*)         '---------------------------------'
  normChg = 0.d0
  DO ia = 1, ATOMS%Natoms
    r0 = ATOMS%positions(:,ia)
    is = ATOMS%atmToSpecies(ia)
    ZZ = ATOMS%Zv(is)
    normChg = normChg + ZZ
    !WRITE(*,*) ia, ZZ  ! DEBUG
    DO ip = 1, N**3
      r = norm2( LF%lingrid(:,ip) - r0(:) )
      ! XXX For the moment only s-like rho
      Rho(ip) = Rho(ip) + ZZ*exp(-r**2/(2.d0*sigma**2))/(2.d0*pi*sigma**2)**1.5d0
    ENDDO
  ENDDO
  integRho = sum(Rho(:))*deltaV
  WRITE(*,*) 'Integrated Rho:', integRho 
  !WRITE(*,*) 'normChg:', normChg   ! DEBUG
  ! Renormalize
  Rho(:) = Rho(:)*normChg/integRho
  WRITE(*,*) 'Renormalized to:', sum(Rho)*deltaV

END SUBROUTINE

