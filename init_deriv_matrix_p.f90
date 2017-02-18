! TODO: using full storage, not yet exploiting symmetry of D1jl and D2jl
!
! The matrices D1jl and D2jl are assumed to have been allocated somewhere else.
!------------------------------------------------------------------------------
SUBROUTINE init_deriv_matrix_p( N, L, D1jl, D2jl )
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: L
  REAL(8) :: D1jl(N,N)
  REAL(8) :: D2jl(N,N)
  ! Local
  INTEGER :: NPRIMED, nn
  INTEGER :: jj, ll
  REAL(8) :: tt1, tt2, tt3, tt4

  !
  ! Diagonal elements
  NPRIMED = (N-1)/2
  DO jj = 1, N
    D1jl(jj,jj) = 0d0
    D2jl(jj,jj) = -( 2.d0*PI/L )**2.d0 * NPRIMED * (NPRIMED+1)/3.d0
  ENDDO
  !
  ! Off diagonal elements
  DO jj = 1, N
    DO ll = jj+1, N
      !
      nn = jj - ll
      !
      tt1 = PI/L * (-1.d0)**nn
      tt2 = sin(PI*nn/N)
      !
      tt3 = (2.d0*PI/L)**2d0 * (-1.d0)**nn * cos(PI*nn/N)
      tt4 = 2.d0*sin(PI*nn/N)**2d0
      !
      D1jl(jj,ll) =  tt1/tt2
      D1jl(ll,jj) = -tt1/tt2
      !
      D2jl(jj,ll) = -tt3/tt4
      D2jl(ll,jj) = -tt3/tt4
    ENDDO
  ENDDO
  
END SUBROUTINE

