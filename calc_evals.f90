SUBROUTINE calc_evals( Nstates, Focc, v, evals )
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, dVol => LF3d_dVol
  IMPLICIT NONE 
  INTEGER :: Nstates
  REAL(8) :: Focc(Nstates)
  REAL(8) :: v(Npoints,Nstates)
  REAL(8) :: evals(Nstates)
  !
  REAL(8), ALLOCATABLE :: Hv(:,:)
  REAL(8), ALLOCATABLE :: Hred(:,:)
  INTEGER :: lwork, info
  REAL(8), ALLOCATABLE :: work(:)

  ALLOCATE( Hv(Npoints,Nstates) )
  ALLOCATE( Hred(Nstates,Nstates) )
  lwork = max(1, 3*Nstates-1)
  ALLOCATE( work(lwork) )

  CALL calc_rhoe( Focc, v )
  CALL update_potentials()
  CALL calc_betaNL_psi( Nstates, v )
  CALL op_H( Nstates, v, Hv )

  CALL dgemm( 'T', 'N', Nstates, Nstates, Npoints, 1.d0, v,Npoints, Hv,Npoints, 0.d0, Hred, Nstates )

  CALL dscal( Nstates*Nstates, dVol, Hred, 1 )

  ! only calculate 
  CALL dsyev( 'N', 'U', Nstates, Hred, Nstates, evals, work, lwork, info )
  IF( info /= 0 ) THEN 
    WRITE(*,*) 'ERROR in dsyev: info = ', info
    STOP 
  ENDIF 

  DEALLOCATE( work )
  DEALLOCATE( Hv )
  DEALLOCATE( Hred )

END SUBROUTINE 

