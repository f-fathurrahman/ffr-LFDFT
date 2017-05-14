!! PURPOSE:
!!
!!   This subroutine solves Schrodinger equation using iterative
!!   diagonalization. Currently this subroutine calls `diag_davidson_qe`.
!!
!! AUTHOR:
!!
!!   Fadjar Fathurrahman
!!
!! MODIFIES:
!!
!!   Global variables `KS_evecs` and `KS_evals`.
!!
!! IMPORTANT
!!
!!   `KS_evecs` should be initialized outside this subroutine.
!!   In the subsequent calls, `KS_evecs` is used as initial guesses.

SUBROUTINE Sch_solve_diag()

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_states, ONLY : Nstates, Focc, &
                       evecs => KS_evecs, &
                       evals => KS_evals
  USE m_options, ONLY : ethr => DIAG_DAVIDSON_QE_ETHR, IALG_DIAG
  IMPLICIT NONE 
  INTEGER, ALLOCATABLE :: btype(:)
  INTEGER :: dav_iter
  INTEGER :: notcnv
  INTEGER :: ist

  ALLOCATE( btype(Nstates) )
  btype(:) = 1 ! assume all bands are occupied
  DO ist = 1,Nstates
    IF( Focc(ist) <= 1d-13 ) btype(ist) = 0
  ENDDO 

  !WRITE(*,*)
  !WRITE(*,*) 'Solving Schrodinger equation with Davidson iterative diagonalization'
  !WRITE(*,*)
 
  IF( IALG_DIAG == 1 ) THEN 
    CALL diag_davidson_qe( Npoints, Nstates, 3*Nstates, evecs, ethr, &
                           evals, btype, notcnv, dav_iter )
  ELSEIF( IALG_DIAG == 2 ) THEN 
    CALL diag_davidson( evals, evecs, ethr )
  ELSEIF( IALG_DIAG == 3 ) THEN 
    CALL diag_lobpcg( evals, evecs, ethr )
  ENDIF 
    
  !WRITE(*,'(1x,A,ES18.10,A,I4)') 'Davidson_QE: ethr = ', ethr, ' dav_iter = ', dav_iter

  WRITE(*,*)
  WRITE(*,*) 'Eigenvalues:'
  WRITE(*,*)
  DO ist = 1, Nstates
    WRITE(*,'(1x,I4,F18.10)') ist, evals(ist)
  ENDDO
  WRITE(*,*)

  ! normalize evecs properly
  !evecs(:,:) = evecs(:,:)/sqrt(dVol)
  
  CALL ortho_check( Npoints, Nstates, dVol, evecs )
  STOP 

  DEALLOCATE( btype )
END SUBROUTINE 

