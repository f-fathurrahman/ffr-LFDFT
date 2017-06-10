! nocc --> Nelectrons
! nev --> Nstates

SUBROUTINE calc_occupations( evals, Tbeta, efermi )

  USE m_states, ONLY : Nstates, Focc, Nelectrons
  IMPLICIT NONE 
  REAL(8) :: evals(Nstates)
  REAL(8) :: Tbeta
  REAL(8) :: efermi
  !
  REAL(8), PARAMETER :: TOL = 1d-15
  INTEGER, PARAMETER :: MAXITER = 100
  INTEGER :: ilb, ulb, iter, iub, ist
  REAL(8) :: lb, ub, flb, fub
  REAL(8), ALLOCATABLE :: Focc_ub(:), Focc_lb(:)
  REAL(8) :: occsum

  ALLOCATE( Focc_ub(Nstates) )
  ALLOCATE( Focc_lb(Nstates) )

  IF( Nstates > Nelectrons ) THEN 
    !
    ! use bisection to find efermi such that 
    !       sum_i fermidirac(ev(i)) = Nelectrons
    !
    ilb = int(Nelectrons) - 1
    iub = int(Nelectrons) + 1
    lb = evals(ilb)
    ub = evals(iub)
    !
    ! make sure flb < Nelectrons and fub > Nelectrons
    !
    CALL fermi_dirac(Nstates, evals, lb, Tbeta, Focc_lb ) 
    flb = sum(Focc_lb)
    !
    CALL fermi_dirac(Nstates, evals, ub, Tbeta, Focc_ub )
    fub = sum(Focc_ub)

    DO WHILE( (Nelectrons-flb)*(fub-Nelectrons) < 0.d0 )

      WRITE(*,*) 'calc_occupations initial bounds are off'
      WRITE(*,*) 'flb, fub, nocc = ', flb, fub, Nelectrons

      IF( flb > Nelectrons ) THEN 
        IF( ilb > 1 ) THEN 
          ilb = ilb - 1
          lb = evals(ilb)
          CALL fermi_dirac( Nstates, evals, lb, Tbeta, Focc_lb )
          flb = sum( Focc_lb )
        ELSE 
          WRITE(*,*) 'calc_occupations cannot find a lower bound for efermi'
          WRITE(*,*) 'something is wrong'
          EXIT 
        ENDIF 
      ENDIF 
      
      IF( fub < Nelectrons ) THEN 
        IF( iub < int(Nstates) ) THEN 
          iub = iub + 1
          ub  = evals(iub)
          CALL fermi_dirac( Nstates, evals, ub, Tbeta, Focc_ub )
          fub = sum(Focc_ub)
        ELSE 
          WRITE(*,*) 'getocc: cannot find an upper bound for efermi,'
          WRITE(*,*) 'something is wrong, try increasing the number of wavefunctions in X0'
          EXIT 
        ENDIF 
      ENDIF 

    ENDDO 
   
    WRITE(*,*) 'flb, fub = ', flb, fub
    efermi = ( lb + ub )/2.d0
    CALL fermi_dirac( Nstates, evals, efermi, Tbeta, Focc )
    occsum = sum(Focc)
    iter = 1

    DO WHILE( abs(occsum - Nelectrons) > tol .AND. iter < maxiter )
      WRITE(*,*) 'iter, efermi, sum = ', iter, efermi, occsum
      WRITE(*,*) 'lb, ub = ', lb, ub
      !
      IF( occsum < Nelectrons ) THEN 
        lb = efermi
      ELSE 
        ub = efermi
      ENDIF 
      !
      efermi = ( lb + ub )/2.d0
      CALL fermi_dirac( Nstates, evals, efermi, Tbeta, Focc )
      occsum = sum(Focc)
      iter = iter + 1
    ENDDO 

  ELSEIF( Nstates == Nelectrons ) THEN 
    DO ist = 1,Nstates
      Focc(ist) = 1.d0
    ENDDO 
    efermi = evals(int(Nelectrons))

  ELSE
    WRITE(*,*) 'ERROR in calc_occupations:'
    WRITE(*,*) 'The number of eigenvalues in ev should be larger than nocc'
    STOP 

  ENDIF 

  DEALLOCATE( Focc_lb )
  DEALLOCATE( Focc_ub )

END SUBROUTINE 

