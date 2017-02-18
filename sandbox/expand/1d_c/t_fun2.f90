FUNCTION gaussian( sigma, mu, x )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  REAL(8) :: gaussian
  REAL(8) :: sigma, mu, x

  gaussian = 1.d0/(sigma*sqrt(2.d0*PI))*exp( -0.5d0*( (x-mu)/sigma )**2 )
END FUNCTION


!------------------------------------------------------------------------------
PROGRAM t_fun
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF1d
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: A, B, L
  !
  TYPE(LF1d_t) :: LF
  REAL(8), ALLOCATABLE :: coefs(:)
  !
  INTEGER :: ii, jj
  ! Functions
  REAL(8) :: gaussian
  !
  REAL(8) :: sigma, mu, mu2, h
  INTEGER :: NPTS_PLOT=300
  REAL(8) :: xx, yy
  INTEGER :: narg
  CHARACTER(32) :: filnam1, filnam2, buffer

  A = -5.d0
  B = 5.d0
  L = B - A

  narg = iargc()
  IF( narg /= 1) THEN
    WRITE(*,*) 'ERROR: need exactly one argument, i.e. number of basis functions'
    STOP
  ENDIF

  ! Get N from program's argument
  CALL getarg(1,buffer)
  READ(buffer,*) N
  WRITE(*,*) 'N = ', N

  ! Sanity check
  IF( N > NPTS_PLOT ) THEN
    WRITE(*,*) 'ERROR: N is larger than NPTS_PLOT.'
    WRITE(*,*) 'Please reduce N or modify NPTS_PLOT in the source.'
    STOP
  ENDIF

  ! Gaussian parameters
  sigma = 0.3d0
  mu    = -2.d0
  mu2   = 2.d0

  ! Initialize the basis functions
  CALL init_LF1d_c( LF, N, A, B )
  CALL info_LF1d( LF, .TRUE. )
  h = LF%h
  
  ! Calculate coefficients of expansion for the function
  ALLOCATE( coefs(N) )

  IF( N < 10 ) THEN
    WRITE(filnam1,'(A,I1)') 'COARSE_N_', N
    WRITE(filnam2,'(A,I1)') 'DENSE_N_', N
  ELSEIF( N < 100 ) THEN
    WRITE(filnam1,'(A,I2)') 'COARSE_N_', N
    WRITE(filnam2,'(A,I2)') 'DENSE_N_', N
  ELSEIF( N < 1000 ) THEN
    WRITE(filnam1,'(A,I3)') 'COARSE_N_', N
    WRITE(filnam2,'(A,I3)') 'DENSE_N_', N
  ELSE
    WRITE(*,'(A,I10)') 'ERROR: N is not supported: ', N
    STOP
  ENDIF

  OPEN(unit=11,file=filnam1,action='write')
  DO ii=1,N
    coefs(ii) =  sqrt(h)*( gaussian( sigma, mu, LF%grid(ii) ) + &
       gaussian( sigma, mu2, LF%grid(ii) ) )
    WRITE(11,*) LF%grid(ii), coefs(ii)/sqrt(h)
  ENDDO
  CLOSE(11)
  
  WRITE(*,*) 'integration test: ', sum( coefs(:)/sqrt(h) )*h

  ! Write to file for plotting
  OPEN(unit=12,file=filnam2,action='write')
  DO ii=1,NPTS_PLOT
    !
    xx = A + (ii-1)*L/(NPTS_PLOT-1)
    yy = 0.d0
    DO jj=1,N
      yy = yy + eval_LF1d_c( LF, jj, xx )*coefs(jj)
    ENDDO
    WRITE(12,*) xx, yy, gaussian( sigma, mu, xx ) + gaussian( sigma, mu2, xx )
  ENDDO
  CLOSE(12)

  ! Free memory
  CALL dealloc_LF1d( LF )
  DEALLOCATE( coefs )


END PROGRAM

