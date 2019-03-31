!!>
!!> \section{Subroutine \texttt{KS\_solve\_SCF}}
!!> 
!!> This is main computational routine for SCF solution of Kohn-Sham equations.
!!>
SUBROUTINE KS_solve_SCF()
  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_states, ONLY : Nstates, Focc, &
                       evecs => KS_evecs
  USE m_options, ONLY : ETHR_EVALS_LAST, &
                        ethr => ETHR_EVALS, &
                        MIXTYPE, &
                        SCF_ETOT_CONV_THR
  USE m_hamiltonian, ONLY : Rhoe
  USE m_energies, ONLY : Etot => E_total
!!> These mixing-related variables should be put into specialized modules
  USE m_options, ONLY : beta0, betamax, mixsdb, broydpm, SCF_betamix, SCF_NiterMax
  IMPLICIT NONE
  !
  INTEGER :: iterSCF
  REAL(8) :: Etot_old, dEtot
  REAL(8) :: dr2
  REAL(8), ALLOCATABLE :: Rhoe_old(:)

  REAL(8) :: integRho
  REAL(8), ALLOCATABLE :: beta_work(:), f_work(:)
  REAL(8), ALLOCATABLE :: workmix(:)
  INTEGER :: nwork
  INTEGER :: Nconverged

  integRho = sum(Rhoe)*dVol
  WRITE(*,*)
  WRITE(*,'(1x,A,F18.10)') 'Initial guess: integRho = ', integRho

  ALLOCATE( Rhoe_old(Npoints) )

  beta0 = SCF_betamix
  betamax = 1.d0
  ! Broyden parameters recommended by M. Meinert
  mixsdb = 4
  broydpm(1) = 0.4d0
  broydpm(2) = 0.15d0

  ALLOCATE( beta_work(Npoints) )
  ALLOCATE( f_work(Npoints) )

  Etot_old = 0.d0
  Rhoe_old(:) = Rhoe(:)

  ! Allocate memory for ELK mixing subroutines
  !
  ! Linear mixing
  IF( MIXTYPE == 0 ) THEN 
    nwork = Npoints
    ALLOCATE( workmix(nwork) )
    WRITE(*,*)
    WRITE(*,*) 'Using linear mixing with beta = ', beta0
    WRITE(*,*)
  !
  ! Adaptive linear mixing
  ELSEIF( MIXTYPE == 1 ) THEN 
    nwork = 3*Npoints
    ALLOCATE( workmix(nwork) )
    WRITE(*,*)
    WRITE(*,*) 'Using adaptive linear mixing'
    WRITE(*,*)
  !
  ! Broyden mixing
  ELSEIF( MIXTYPE == 3 ) THEN 
    nwork = (4+2*mixsdb)*Npoints + mixsdb**2
    ALLOCATE( workmix(nwork) )
    WRITE(*,*)
    WRITE(*,*) 'Using Broyden mixing'
    WRITE(*,*) 'Broyden mixing nwork = ', nwork
    WRITE(*,*)
  ELSE 
    WRITE(*,*)
    WRITE(*,*) 'WARNING: Unknown MIXTYPE: ', MIXTYPE
    WRITE(*,*) 'Switching to default: adaptive linear mixing'
    MIXTYPE = 1
    nwork = 3*Npoints
    ALLOCATE( workmix(nwork) )
  ENDIF 

  flush(6)

!!> Initial value for dr2 (convergence criteria for rhoe ???)
  dr2 = 1.d0

  Nconverged = 0

!!>
!!> SCF iterations begins here
!!>
  DO iterSCF = 1, SCF_NiterMax

    WRITE(*,*)
    WRITE(*,'(1x,A,I5,A)') '*** Begin SCF iter: ', iterSCF, '***'

!!>
!!> Determine convergence criteria for iterative diagonalization
!!>
    IF( iterSCF==1 ) THEN
      ethr = 1.d-1
    ELSE 
      IF( iterSCF == 2 ) ethr = 1.d-2
      ethr = ethr/5.d0
      ethr = max( ethr, ETHR_EVALS_LAST )
    ENDIF 
    ethr = 1.d-8  ! temporary fix for getting accurate results

!!> Call the driver for iterative diagoanalization routine
    CALL Sch_solve_diag()

!!> Calculate electron density
    CALL calc_rhoe( Focc, evecs )

!!> Mix electron density
    CALL mixerifc( iterSCF, MIXTYPE, Npoints, Rhoe, dr2, nwork, workmix )

!!> Make sure that there are no negative densities
    CALL normalize_rhoe( Npoints, Rhoe )  ! make sure no negative or very small rhoe

!!> Update potentials and evaluate total energy.
    CALL update_potentials()
    CALL calc_betaNL_psi( Nstates, evecs )
    CALL calc_energies( evecs ) ! update the potentials or not ?
    
    dEtot = abs(Etot - Etot_old)

    WRITE(*,*)
    WRITE(*,'(1x,A,I5,F18.10,2ES18.10)') 'SCF: ', iterSCF, Etot, dEtot, dr2

!!> Check total energy convergence. Exit SCF loop if convergence is achieved.
    IF( dEtot < SCF_ETOT_CONV_THR ) THEN
      Nconverged = Nconverged + 1
    ELSE 
      Nconverged = 0
    ENDIF 

    IF( Nconverged >= 2 ) THEN 
      WRITE(*,*)
      WRITE(*,'(1x,A,I5,A)') 'SCF converged after ', iterSCF, ' iterations.'
      WRITE(*,'(1x,A,ES18.10)') 'SCF_ETOT_CONV_THR = ', SCF_ETOT_CONV_THR
      WRITE(*,'(1x,A,ES18.10)') 'dEtot             = ', dEtot
      EXIT 
    ENDIF 

!!> Convergence is not achieved. Prepare for the next SCF iteration.
!!> Save previous Rhoe and total energy.

    Etot_old = Etot
    Rhoe_old(:) = Rhoe(:)
    flush(6)
  ENDDO

  IF( Nconverged < 2 ) THEN 
    WRITE(*,*)
    WRITE(*,*) 'WARNING: SCF is not converged'
    WRITE(*,*)
  ENDIF 

!!> Free memory

  DEALLOCATE( workmix )
  DEALLOCATE( Rhoe_old )
  DEALLOCATE( beta_work )
  DEALLOCATE( f_work )

END SUBROUTINE 

