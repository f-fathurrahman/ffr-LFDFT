! Similar to calc_evals, but without updating potentials
SUBROUTINE calc_Ebands( Nstates, v, evals, Ebands )
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, dVol => LF3d_dVol
  IMPLICIT NONE 
  INTEGER :: Nstates
  REAL(8) :: v(Npoints,Nstates)
  REAL(8) :: evals(Nstates)  ! output
  REAL(8) :: Ebands
  !
  REAL(8), ALLOCATABLE :: Hv(:,:)
  REAL(8), ALLOCATABLE :: Hred(:,:)
  INTEGER :: lwork, info
  REAL(8), ALLOCATABLE :: work(:)

  ALLOCATE( Hv(Npoints,Nstates) )
  ALLOCATE( Hred(Nstates,Nstates) )
  lwork = max(1, 3*Nstates-1)
  ALLOCATE( work(lwork) )

  CALL op_H( Nstates, v, Hv )

  CALL dgemm( 'T', 'N', Nstates, Nstates, Npoints, 1.d0, v,Npoints, Hv,Npoints, 0.d0, Hred, Nstates )

  CALL dscal( Nstates*Nstates, dVol, Hred, 1 )

  ! only calculate 
  CALL dsyev( 'N', 'U', Nstates, Hred, Nstates, evals, work, lwork, info )
  IF( info /= 0 ) THEN 
    WRITE(*,*) 'ERROR in dsyev: info = ', info
    STOP 
  ENDIF 

  Ebands = sum(evals)

  DEALLOCATE( work )
  DEALLOCATE( Hv )
  DEALLOCATE( Hred )
END SUBROUTINE 
