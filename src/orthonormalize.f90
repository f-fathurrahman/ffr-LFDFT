SUBROUTINE orthonormalize( Ncols, v )
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  IMPLICIT NONE 
  INTEGER :: Ncols
  REAL(8) :: v(Npoints,Ncols)
  
  CALL ortho_gram_schmidt( v, Npoints, Npoints, Ncols )
  v(:,:) = v(:,:)/sqrt(dVol)

END SUBROUTINE

