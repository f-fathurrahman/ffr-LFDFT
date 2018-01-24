! Gaussian potential in R-space:
! V(r) = sum_{i} A_{i} \exp( -\alpha_{i} r^{2} )
!
SUBROUTINE init_V_ps_loc_gaussian( Nparams, A, alpha )
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     lingrid => LF3d_lingrid
  USE m_hamiltonian, ONLY : V_ps_loc
  USE m_atoms, ONLY : Nspecies, atpos => AtomicCoords, Natoms, atm2species
  IMPLICIT NONE 
  !
  INTEGER :: Nparams
  REAL(8) :: A(Nparams)
  REAL(8) :: alpha(Nparams)
  !
  INTEGER :: ip, isp, ia
  REAL(8) :: r

  IF( Nspecies /= Nparams ) THEN 
    WRITE(*,*) 'ERROR in init_V_ps_loc_gaussian'
    WRITE(*,'(1x,A,2I4)') 'Nspecies /= Nparams : ', Nspecies, Nparams
    STOP 
  ENDIF 

  WRITE(*,*)
  WRITE(*,*) 'Initializing V_ps_loc with Gaussian potential (real space)'

  DO ia = 1,Natoms
    isp = atm2species(ia)
    DO ip = 1, Npoints
      CALL calc_dr_1pnt( atpos(:,ia), lingrid(:,ip), r )
      V_ps_loc(ip) = V_ps_loc(ip) - A(isp)*exp( -alpha(isp)*r**2 )
    ENDDO 
  ENDDO 

  WRITE(*,*) 'sum(V_ps_loc) = ', sum(V_ps_loc)

  flush(6)

END SUBROUTINE 

