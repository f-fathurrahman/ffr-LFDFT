! Fadjar Fathurrahman (20910015), May 2011
! Updated: February 2017
! Updated: June 2017 (using module)

MODULE m_fftw3

  IMPLICIT NONE 

  INTEGER(8) :: PLAN_FORWARD
  INTEGER(8) :: PLAN_BACKWARD

CONTAINS 

SUBROUTINE init_fftw3_plans( zdata, Nx, Ny, Nz )
  IMPLICIT NONE 
  include 'fftw3.f'
  INTEGER :: Nx, Ny, Nz
  COMPLEX(8) :: zdata(Nx,Ny,Nz)

  CALL dfftw_plan_dft_3d( PLAN_FORWARD, Nx, Ny, Nz, zdata, &
                          zdata, FFTW_FORWARD,FFTW_ESTIMATE )

  CALL dfftw_plan_dft_3d( PLAN_BACKWARD, Nx, Ny, Nz, zdata, &
                          zdata, FFTW_BACKWARD,FFTW_ESTIMATE )
END SUBROUTINE 


SUBROUTINE destroy_fftw3_plans()
  IMPLICIT NONE 
  CALL dfftw_destroy_plan( PLAN_FORWARD )
  CALL dfftw_destroy_plan( PLAN_BACKWARD )
END SUBROUTINE 

!
! Calling FFTW from legacy Fortran.
! In place FFT using FFTW3: `zdata` is overwritten
!-----------------------------------------------------
SUBROUTINE exec_fft_fftw3( zdata, Nx, Ny, Nz, t_inv)
!-----------------------------------------------------
  IMPLICIT NONE
  INCLUDE 'fftw3.f'
  ! Arguments
  INTEGER :: Nx, Ny, Nz
  COMPLEX(8) :: zdata(Nx, Ny, Nz)
  LOGICAL :: t_inv
  REAL(8) :: scal

  IF( t_inv ) THEN ! backward transform
    CALL dfftw_execute_dft( PLAN_BACKWARD, zdata, zdata )
  !
  ELSE ! forward transform
    CALL dfftw_execute_dft( PLAN_FORWARD, zdata, zdata )
    ! Scale the result
    scal = 1.d0/REAL( Nx*Ny*Nz )
    CALL zdscal( Nx*Ny*Nz, scal, zdata, 1 )
  ENDIF

END SUBROUTINE


END MODULE 

