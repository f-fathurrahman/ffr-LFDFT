! Initialize periodic Lagrange function
! The LFs are defined on [A,B].
! `grid` is assumed to be already allocated properly elsewhere
!------------------------------------------------------------------------------
SUBROUTINE init_grid_1d_p( N, A, B, grid )
!------------------------------------------------------------------------------
  IMPLICIT NONE
  ! Arguments
  INTEGER :: N
  REAL(8) :: A, B
  REAL(8) :: grid(N)
  ! Local
  INTEGER :: ii

  ! Check if N is odd
  IF( mod(N,2)==0 ) THEN
    WRITE(*,*) 'N should be an odd number, this N=', N
    STOP
  ENDIF

  ! Check if B > A
  IF( B < A ) THEN
    WRITE(*,*) 'Error in init_LF1d_p: B should be larger that A', A, B
    STOP
  ENDIF

  ! Initialization of grid points
  DO ii=1,N
    grid(ii) = A + 0.5d0*(B-A)*(2*ii-1)/N
  ENDDO

END SUBROUTINE
