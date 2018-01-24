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
  USE m_options, ONLY : ethr => ETHR_EVALS, I_ALG_DIAG
  USE m_input_vars, ONLY : ortho_check_after_diag
  IMPLICIT NONE 
  INTEGER, ALLOCATABLE :: btype(:)
  INTEGER :: dav_iter
  INTEGER :: notcnv
  INTEGER :: ist

  ! This is required by \texttt{diag_davidson_qe}.
  ALLOCATE( btype(Nstates) )
  btype(:) = 1 ! by default all bands are treated with the same accuracy
  ! For states with small occupancy, we treat is as unoccupied bands
  DO ist = 1,Nstates
    IF( Focc(ist) <= 1d-13 ) btype(ist) = 0
  ENDDO 

  ! Note that, for some diagonalization routines, evecs are assumed to
  ! be normalized according to:
  !     sum( evecs(:,ist)*evecs(:,ist) ) = 1
  !
  ! This will be treated differently below.

  IF( I_ALG_DIAG == 1 ) THEN 

    evecs = evecs(:,:)*sqrt(dVol)  ! normalize
    CALL diag_davidson_qe( Npoints, Nstates, 3*Nstates, evecs, ethr, &
                           evals, btype, notcnv, dav_iter )
    evecs(:,:) = evecs(:,:)/sqrt(dVol)

  ELSEIF( I_ALG_DIAG == 2 ) THEN 
 
    ! No need to (re)normalize the eigenvectors
    ! Factor dVol is already incorpared in \texttt{diag_davidson}
    CALL diag_davidson( evals, evecs, ethr, .FALSE. )

  ELSEIF( I_ALG_DIAG == 3 ) THEN 

    evecs = evecs(:,:)*sqrt(dVol)  ! normalize
    CALL diag_lobpcg( evals, evecs, ethr, .FALSE. )
    evecs(:,:) = evecs(:,:)/sqrt(dVol)

  ENDIF 
  
  WRITE(*,*)
  WRITE(*,*) 'End of iterative diagonalization'
  WRITE(*,*)
  WRITE(*,*) 'Final eigenvalues (in Hartree)'
  WRITE(*,*)
  DO ist = 1, Nstates
    WRITE(*,'(1x,I4,F18.10)') ist, evals(ist)
  ENDDO

  IF( ortho_check_after_diag ) THEN 
    CALL ortho_check( Npoints, Nstates, dVol, evecs )
  ENDIF 

  DEALLOCATE( btype )

  flush(6)

END SUBROUTINE 

