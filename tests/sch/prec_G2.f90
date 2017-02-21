SUBROUTINE prec_G2( Ncols, v, prec_v )

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     NN => LF3d_NN, &
                     G2 => LF3d_G2
  IMPLICIT NONE 
  INTEGER :: Ncols
  REAL(8) :: v(Npoints,Ncols)
  REAL(8) :: prec_v(Npoints,Ncols)
  !
  COMPLEX(8), ALLOCATABLE :: ctmp(:)
  INTEGER :: ig, ic

  ALLOCATE( ctmp(Npoints) )

  DO ic = 1, Ncols
    ctmp(:) = v(:,ic)
    WRITE(*,*)
    WRITE(*,*) 'Before forward FT: ', sum(ctmp)
    CALL fft_fftw3( ctmp, NN(1), NN(2), NN(3), .FALSE. )
    WRITE(*,*) 'After forward FT: ', sum(ctmp)
    DO ig = 1, Npoints
      ctmp(ig) = ctmp(ig) / ( 1.d0 + G2(ig) )
    ENDDO
    WRITE(*,*) 'After prec: ', sum(ctmp)
    CALL fft_fftw3( ctmp, NN(1), NN(2), NN(3), .TRUE. )
    WRITE(*,*) 'After inverse FT: ', sum(ctmp)
    prec_v(:,ic) = real( ctmp(:), kind=8 )
    WRITE(*,*) ic, sum(v(:,ic)), sum(prec_v(:,ic))
  ENDDO

  DEALLOCATE( ctmp )
END SUBROUTINE
