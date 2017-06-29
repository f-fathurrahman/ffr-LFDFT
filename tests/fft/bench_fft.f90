PROGRAM bench_fft
  USE m_fftw3
  IMPLICIT NONE 
  INTEGER :: NN(3)
  INTEGER :: Npoints
  COMPLEX(8), ALLOCATABLE :: dat1(:)
  INTEGER :: Ntrials, i
  !
  INTEGER :: tstart, counts_per_second, tstop

  NN(:) = (/ 55, 55, 55 /)
  Npoints = NN(1)*NN(2)*NN(3)

  ALLOCATE( dat1(Npoints) )

  dat1(:) = cmplx( 1.d0, 2.3d0 )

  Ntrials = 1000

  CALL system_clock( tstart, counts_per_second )
  CALL init_fftw3_plans( dat1, NN(1), NN(2), NN(3) )
  DO i = 1,Ntrials
    CALL exec_fft_fftw3( dat1, NN(1), NN(2), NN(3), .FALSE. )  ! forward transform
    CALL exec_fft_fftw3( dat1, NN(1), NN(2), NN(3), .TRUE. )   ! backward transform
  ENDDO 
  CALL destroy_fftw3_plans()
  CALL system_clock( tstop )

  WRITE(*,*) 'Total elapsed time: ', dble(tstop - tstart)/counts_per_second, ' seconds.'

  DEALLOCATE( dat1 )

END PROGRAM 

