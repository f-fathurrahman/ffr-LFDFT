!
! Minimal subroutine to generate magnitude of G-vectors
!
SUBROUTINE init_gvec()
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : NN => LF3d_NN, &
                       LL => LF3d_LL, &
                       G2 => LF3d_G2, &
                       Npoints => LF3d_Npoints
  IMPLICIT NONE
  !
  REAL(8) :: gg(3)
  INTEGER :: i, j, k, ig, ii, jj, kk
  ! Function
  INTEGER :: mm_to_nn

  ALLOCATE( G2(Npoints) )
  
  ig = 0
  DO k = 0, NN(3)-1
  DO j = 0, NN(2)-1
  DO i = 0, NN(1)-1
    !
    ig = ig + 1
    !
    ii = mm_to_nn( i, NN(1) )
    jj = mm_to_nn( j, NN(2) )
    kk = mm_to_nn( k, NN(3) )
    !
    gg(1) = ii * 2.d0*PI/LL(1)
    gg(2) = jj * 2.d0*PI/LL(2)
    gg(3) = kk * 2.d0*PI/LL(3)
    !
    G2(ig) = gg(1)**2 + gg(2)**2 + gg(3)**2
  ENDDO
  ENDDO
  ENDDO

END SUBROUTINE
