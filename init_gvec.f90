!
! A subroutine to generate of G-vectors
! Not using minimal isotropic representation
!
SUBROUTINE init_gvec()
  USE m_constants, ONLY : PI
  USE m_LF3d, ONLY : NN => LF3d_NN, &
                     LL => LF3d_LL, &
                     G2 => LF3d_G2, &
                     Gv => LF3d_Gv, &
                     Npoints => LF3d_Npoints
  IMPLICIT NONE
  !
  INTEGER :: i, j, k, ig, ii, jj, kk
  ! Function
  INTEGER :: mm_to_nn

  ALLOCATE( G2(Npoints) )
  ALLOCATE( Gv(3,Npoints) )
  
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
    Gv(1,ig) = ii * 2.d0*PI/LL(1)
    Gv(2,ig) = jj * 2.d0*PI/LL(2)
    Gv(3,ig) = kk * 2.d0*PI/LL(3)
    !
    G2(ig) = Gv(1,ig)**2 + Gv(2,ig)**2 + Gv(3,ig)**2
  ENDDO
  ENDDO
  ENDDO

END SUBROUTINE


