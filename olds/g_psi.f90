!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define TEST_NEW_PRECONDITIONING
!
!-----------------------------------------------------------------------
subroutine g_psi (lda, n, m, psi, e)
  !-----------------------------------------------------------------------
  !
  !    This routine computes an estimate of the inverse Hamiltonian
  !    and applies it to m wavefunctions
  !
  USE m_globals, ONLY : h_diag
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = 8
  INTEGER :: lda, n, m
  ! input: the leading dimension of psi
  ! input: the real dimension of psi
  ! input: the number of bands
  REAL(DP) :: e (m)
  ! input: the eigenvectors
  REAL(DP) :: psi (lda, m)
  ! inp/out: the psi vector
  !
  !    Local variables
  !
  REAL(DP), parameter :: eps = 1.0d-4
  ! a small number
  REAL(DP) :: x, scala, denm
  INTEGER :: k, i
  ! counter on psi functions
  ! counter on G vectors
  !
  !
#ifdef TEST_NEW_PRECONDITIONING
  scala = 1.d0
  DO k = 1, m
    DO i = 1, n
      x = (h_diag(i) - e(k))*scala
      denm = 0.5_dp*(1.d0+x+sqrt(1.d0+(x-1)*(x-1.d0)))/scala
      psi (i, k) = psi (i, k) / denm
    ENDDO
  ENDDO
#else
  DO k = 1, m
    DO i = 1, n
      denm = h_diag (i) - e (k)
      IF (abs (denm) < eps) denm = sign (eps, denm)
      psi (i, k) = psi (i, k) / denm
    ENDDO
  ENDDO
#endif

  RETURN
END SUBROUTINE g_psi
