SUBROUTINE op_V_ps_NL( Nstates, psi, Vpsi )
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_PsPot
  USE m_atoms, ONLY : atm2species, Natoms
  IMPLICIT NONE 
  INTEGER :: Nstates
  REAL(8) :: psi(Npoints,Nstates)
  REAL(8) :: Vpsi(Npoints,Nstates)
  INTEGER :: ia, isp, ibeta, ist, iprjl

  Vpsi(:,:) = 0.d0

  DO ist = 1,Nstates
    DO ia = 1,Natoms
      isp = atm2species(ia)
      DO ibeta = 1,NbetaNL
        ! iprjl = beta2prjl(ibeta,ia)
        iprjl = 1
        Vpsi(:,ist) = w_NL(isp,iprjl)*betaNL(:,ibeta)*betaNL_psi(ia,ist,ibeta)*dVol
      ENDDO 
    ENDDO 
    Vpsi(:,ist) = 2.d0*Vpsi(:,ist)
  ENDDO 

END SUBROUTINE 

