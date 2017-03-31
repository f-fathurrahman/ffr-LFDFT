SUBROUTINE update_potentials()

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     LF3d_TYPE, LF3d_PERIODIC
  USE m_hamiltonian, ONLY : Rhoe, V_Hartree, V_xc

  IMPLICIT NONE 
  REAL(8), ALLOCATABLE :: epsxc(:), depsxc(:)
  
  ALLOCATE( epsxc(Npoints) )
  ALLOCATE( depsxc(Npoints) )

  IF ( LF3d_TYPE == LF3d_PERIODIC ) THEN
    CALL Poisson_solve_fft( Rhoe, V_Hartree )
  ELSE 
    CALL Poisson_solve_cg( Rhoe, V_Hartree )
  ENDIF

  CALL excVWN( Npoints, Rhoe, epsxc )
  CALL excpVWN( Npoints, Rhoe, depsxc )

  V_xc(:) = epsxc(:) + Rhoe(:)*depsxc(:)

  DEALLOCATE( epsxc )
  DEALLOCATE( depsxc )

END SUBROUTINE 
