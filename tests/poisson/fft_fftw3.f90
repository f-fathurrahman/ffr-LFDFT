! Fadjar Fathurrahman (20910015), May 2011

! In place FFT using FFTW3
!-----------------------------------------------------
SUBROUTINE fft_fftw3( zdata, NX, NY, NZ, t_inv)
!-----------------------------------------------------
  IMPLICIT NONE
  INCLUDE 'fftw3.f'
  ! ARGUMENTS
  INTEGER :: NX, NY, NZ
  COMPLEX(8) :: zdata(NX,NY,NZ)
  LOGICAL :: t_inv

  ! TODO: Probably it is better to save PLAN_BACKWARD and PLAN_FORWARD
  !       as global variables
  INTEGER(4) :: PLAN_BACKWARD
  INTEGER(4) :: PLAN_FORWARD
 
  IF( t_inv ) THEN ! backward transform
    CALL dfftw_plan_dft_3d( PLAN_BACKWARD, NZ,NY,NX, zdata, &
                            zdata, FFTW_BACKWARD, FFTW_ESTIMATE )
    CALL dfftw_execute( PLAN_BACKWARD )
    CALL dfftw_destroy( PLAN_BACKWARD )
  !
  ELSE ! forward transform
    CALL dfftw_plan_dft_3d( PLAN_FORWARD, NZ, NY, NX, zdata, &
                            DATAIN, FFTW_FORWARD,FFTW_ESTIMATE )
    CALL dfftw_execute( PLAN_FORWARD )
    ! SCALE THE RESULT
    zdata = zdata / ( NX*NY*NZ ) !TODO: Use zdscal
    CALL dfftw_destroy( PLAN_FORWARD )
  ENDIF

END SUBROUTINE

