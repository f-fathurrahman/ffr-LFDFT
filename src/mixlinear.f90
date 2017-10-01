
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! slightly modified by Fadjar Fathurrahman for ffr-LFDFT

SUBROUTINE mixlinear( iscl, beta, n, nu, mu, d )
  IMPLICIT NONE 
  ! arguments
  INTEGER, INTENT(in) :: iscl
  REAL(8), INTENT(in) :: beta
  INTEGER, INTENT(in) :: n
  REAL(8), INTENT(inout) :: nu(n),mu(n)
  REAL(8), intent(out) :: d
  ! local variables
  INTEGER :: i
  real(8) :: t0,t1
  IF( n <= 0 ) RETURN 

  WRITE(*,*)
  WRITE(*,*) 'Mixing: ELK-linear: beta = ', beta
  WRITE(*,*)

  ! initialise mixer
  IF( iscl <= 0 ) THEN 
    mu(:) = nu(:)
    d = 1.d0
    RETURN 
  ENDIF 
  t0 = 1.d0 - beta
  d = 0.d0
  DO i=1,n
    t1 = nu(i) - mu(i)
    nu(i) = beta*nu(i) + t0*mu(i)
    d = d + t1**2
    mu(i) = nu(i)
  ENDDO 
  d = sqrt( d/dble(n) )
  RETURN 
END SUBROUTINE 

