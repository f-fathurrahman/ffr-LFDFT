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
  USE m_states, ONLY : Nstates, &
                       evecs => KS_evecs, &
                       evals => KS_evals
  USE m_options, ONLY : ethr => DIAG_DAVIDSON_QE_ETHR
  IMPLICIT NONE 
  INTEGER, ALLOCATABLE :: btype(:)
  INTEGER :: dav_iter
  INTEGER :: notcnv
  INTEGER :: ist

  ! starting from random eigenvectors
!  DO ist = 1, Nstates
!    DO ip = 1, Npoints
!      CALL random_number( evecs(ip,ist) )
!    ENDDO
!  ENDDO
!  CALL orthonormalize( Nstates, evecs )

  ALLOCATE( btype(Nstates) )
  btype(:) = 1 ! all bands are occupied

  WRITE(*,*)
  WRITE(*,*) 'Solving Schrodinger equation with Davidson iterative diagonalization'
  WRITE(*,*)
    
  CALL diag_davidson_qe( Npoints, Nstates, 4*Nstates, evecs, ethr, &
                         evals, btype, notcnv, dav_iter )
    
  WRITE(*,'(1x,A,ES18.10,A,I4)') 'Davidson_QE: ethr = ', ethr, ' dav_iter = ', dav_iter

  WRITE(*,*)
  WRITE(*,*) 'Eigenvalues:'
  WRITE(*,*)
  DO ist = 1, Nstates
    WRITE(*,'(1x,I4,F18.10)') ist, evals(ist)
  ENDDO

  ! normalize evecs properly
  evecs(:,:) = evecs(:,:)/sqrt(dVol)

  DEALLOCATE( btype )
END SUBROUTINE 

