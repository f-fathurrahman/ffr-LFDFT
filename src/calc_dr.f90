!! PURPOSE:
!!
!!  This subroutine calculates distance between grid points relative
!!  to a center.
!!
!! AUTHOR:
!!
!!   Fadjar Fathurrahman
!!

SUBROUTINE calc_dr( r0, Npoints, lingrid, dr )

  IMPLICIT NONE
  !
  REAL(8) :: r0(3)
  INTEGER :: Npoints
  REAL(8) :: lingrid(3,Npoints)
  REAL(8) :: dr(Npoints)
  !
  INTEGER :: ip
  REAL(8) :: dx2, dy2, dz2

  DO ip = 1, Npoints
    dx2 = abs( lingrid(1,ip) - r0(1) )**2
    dy2 = abs( lingrid(2,ip) - r0(2) )**2
    dz2 = abs( lingrid(3,ip) - r0(3) )**2
    dr(ip) = sqrt( dx2 + dy2 + dz2 )
  ENDDO

END SUBROUTINE

SUBROUTINE calc_dr_1pnt( r0, lingrid, dr )

  IMPLICIT NONE
  !
  REAL(8) :: r0(3)
  REAL(8) :: lingrid(3)
  REAL(8) :: dr
  !
  REAL(8) :: dx2, dy2, dz2

  dx2 = abs( lingrid(1) - r0(1) )**2
  dy2 = abs( lingrid(2) - r0(2) )**2
  dz2 = abs( lingrid(3) - r0(3) )**2
  dr  = sqrt( dx2 + dy2 + dz2 )

END SUBROUTINE

