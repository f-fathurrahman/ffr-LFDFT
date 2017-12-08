!!>
!!> \section{Subroutine \texttt{update\_potentials}}
!!>
!!> This subroutine updates (calculates) Hartree and XC potentials for
!!> a given electronic density.
!!>
SUBROUTINE update_potentials()
  
  USE m_options, ONLY : I_POISSON_SOLVE
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     LF3d_TYPE, LF3d_PERIODIC
  USE m_hamiltonian, ONLY : Rhoe, V_Hartree, V_xc

  IMPLICIT NONE 

  IF ( LF3d_TYPE == LF3d_PERIODIC ) THEN
    CALL Poisson_solve_fft( Rhoe, V_Hartree )
  ELSE 
    !CALL Poisson_solve_pcg( Rhoe, V_Hartree )
    !CALL Poisson_solve_fft_MT( Rhoe, V_Hartree )
    IF( I_POISSON_SOLVE==1 ) THEN 
      CALL Poisson_solve_ISF( Rhoe, V_Hartree )  ! for Lagrange-sinc functions
    ELSEIF( I_POISSON_SOLVE == 2 ) THEN 
      CALL Poisson_solve_DAGE( Rhoe, V_Hartree ) ! for Lagrange-sinc functions
    ELSE 
      WRITE(*,*)
      WRITE(*,*) 'ERROR: Unknown I_POISSON_SOLVE = ', I_POISSON_SOLVE
      STOP 
    ENDIF 
  ENDIF

  CALL calc_Exc_Vxc()

!  WRITE(*,*) 'sum(V_Hartree) = ', sum(V_Hartree)

END SUBROUTINE 

