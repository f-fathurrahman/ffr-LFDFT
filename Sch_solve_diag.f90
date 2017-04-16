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

  IMPLICIT NONE 
  INTEGER, ALLOCATABLE :: btype(:)
  INTEGER, PARAMETER :: Nstages = 3
  REAL(8) :: ethr(Nstages)
  INTEGER :: dav_iter(Nstages)
  INTEGER :: notcnv
  INTEGER :: istage, ist !, ip
  INTEGER :: NiterTot

  ! starting from random eigenvectors
!  DO ist = 1, Nstates
!    DO ip = 1, Npoints
!      CALL random_number( evecs(ip,ist) )
!    ENDDO
!  ENDDO
!  CALL orthonormalize( Nstates, evecs )

  ALLOCATE( btype(Nstates) )
  btype(:) = 1 ! all bands are occupied

  ethr(1) = 1.d-1
  ethr(2) = 1.d-3
  ethr(3) = 1.d-6

  NiterTot = 0

  WRITE(*,*)
  WRITE(*,*) 'Solving Schrodinger equation with Davidson iterative diagonalization'
  WRITE(*,*)
  DO istage = 1, Nstages
    
    CALL diag_davidson_qe( Npoints, Nstates, 4*Nstates, evecs, ethr(istage), &
                           evals, btype, notcnv, dav_iter(istage) )
    
    WRITE(*,fmt=999) 'Davidson stage: ', istage, ' ethr = ', ethr(istage), ' dav_iter = ', dav_iter(istage)

    NiterTot = NiterTot + dav_iter(istage)

  ENDDO
  WRITE(*,'(1x,A,I10)') 'NiterTot = ', NiterTot
 
  999 FORMAT(1x,A,I2,A,E9.4,A,I5)

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

