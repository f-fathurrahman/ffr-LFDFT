
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! Adapted for ffr-LFDFT by Fadjar Fathurrahman (2017)

SUBROUTINE mixerifc( iscl, mtype, n, v, dv, nwork, work )
  USE m_options, ONLY : beta0, betamax, mixsdb, broydpm
  IMPLICIT NONE
  ! arguments
  INTEGER, INTENT(in) :: iscl
  INTEGER, INTENT(in) :: mtype,n
  REAL(8), INTENT(inout) :: v(n)
  REAL(8), INTENT(out) :: dv
  INTEGER, INTENT(inout) :: nwork
  REAL(8), INTENT(inout) :: work(*)

  ! local variables
  SELECT CASE(mtype)
  CASE(0)
  ! linear mixing
    IF(nwork <= 0) THEN
      nwork = n
      RETURN 
    ENDIF 
    call mixlinear(iscl,beta0,n,v,work,dv)
  CASE(1)
    ! adaptive linear mixing
    IF(nwork <= 0) THEN 
      nwork = 3*n
      RETURN
    ENDIF
    CALL mixadapt(iscl,beta0,betamax,n,v,work,work(n+1),work(2*n+1),dv)
  CASE(3)
  ! Broyden mixing
    IF(nwork <= 0) THEN 
      nwork = (4+2*mixsdb)*n+mixsdb**2
      RETURN 
    ENDIF
    CALL mixbroyden(iscl,n,mixsdb,broydpm(1),broydpm(2),v,work,work(2*n+1), &
                    work(4*n+1),work((4+mixsdb)*n+1),work((4+2*mixsdb)*n+1),dv)
  CASE DEFAULT 
    WRITE(*,*)
    WRITE(*,'("Error(mixerifc): mtype not defined : ",I8)') mtype
    WRITE(*,*)
    STOP 
  END SELECT
  RETURN 
END SUBROUTINE 

subroutine getmixdata(mtype,mixdescr)
implicit none
! arguments
integer, intent(in) :: mtype
character(*), intent(out) :: mixdescr
select case(mtype)
case(0)
  mixdescr='Linear mixing'
case(1)
  mixdescr='Adaptive linear mixing'
case(3)
  mixdescr='Broyden mixing, J. Phys. A: Math. Gen. 17, L317 (1984)'
case default
  write(*,*)
  write(*,'("Error(getmixdata): mixtype not defined : ",I8)') mtype
  write(*,*)
  stop
end select
return
end subroutine

