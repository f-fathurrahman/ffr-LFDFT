!! PURPOSE:
!!
!!  This subroutine calculates distance between grid points relative
!!  to a center, taking into account periodicity.
!!
!! AUTHOR:
!!
!!   Fadjar Fathurrahman
!!
!! NOTE:
!!
!!   This will only works for orthorombic cell
!!

SUBROUTINE calc_dr_periodic_1pnt( LL, r0, lingrid, dr_vec )

  IMPLICIT NONE
  !
  REAL(8) :: LL(3)
  REAL(8) :: r0(3)
  REAL(8) :: lingrid(3)
  REAL(8) :: dr_vec(3)
  !
  REAL(8) :: Lx, Ly, Lz
  REAL(8) :: xx1, xx2, xx3
  REAL(8) :: yy1, yy2, yy3
  REAL(8) :: zz1, zz2, zz3

  Lx = LL(1)
  Ly = LL(2)
  Lz = LL(3)

  xx1 = lingrid(1) - r0(1)
  xx2 = r0(1) + Lx - lingrid(1)
  xx3 = lingrid(1) - r0(1) - Lx
  IF( abs(xx1) < abs(xx2) ) THEN
    dr_vec(1) = xx1
  ELSE
    dr_vec(1) = xx2
  ENDIF 
  IF( abs(dr_vec(1) ) > abs(xx3) ) dr_vec(1) = xx3

  
  yy1 = lingrid(2) - r0(2)
  yy2 = r0(2) + Ly - lingrid(2)
  yy3 = lingrid(2) - r0(2) - Ly
  IF( abs(yy1) < abs(yy2) ) THEN
    dr_vec(2) = yy1
  ELSE
    dr_vec(2) = yy2
  ENDIF 
  IF( abs(dr_vec(2) ) > abs(yy3) ) dr_vec(2) = yy3

  
  zz1 = lingrid(3) - r0(3)
  zz2 = r0(3) + Lz - lingrid(3)
  zz3 = lingrid(3) - r0(3) - Lz
  IF( abs(zz1) < abs(zz2) ) THEN
    dr_vec(3) = zz1
  ELSE
    dr_vec(3) = zz2
  ENDIF 
  IF( abs(dr_vec(3) ) > abs(zz3) ) dr_vec(3) = zz3

END SUBROUTINE

