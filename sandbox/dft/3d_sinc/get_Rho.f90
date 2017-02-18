SUBROUTINE get_Rho( Nbasis, Nstates, evecs, Rho )
  USE m_globals, ONLY : LF, deltaV
  IMPLICIT NONE
  INTEGER :: Nbasis, Nstates
  REAL(8) :: evecs(Nbasis,Nstates)
  REAL(8) :: Rho(Nbasis)
  !
  INTEGER :: is
  REAL(8) :: ddot

  Rho(:) = 0.d0

  WRITE(*,'(/,1x,A)') 'Calculating rho (electron density):'
  WRITE(*,*)          '-----------------------------------'
  DO is = 1, Nstates
    WRITE(*,'(1x,A,I5,F18.10)') 'is, norm:', is, sqrt( ddot( Nbasis, evecs(:,is),1, evecs(:,is),1 ) )
    Rho(:) = Rho(:) + 2.d0 * evecs(:,is)*evecs(:,is)  ! assume doubly occupied states
  ENDDO
  WRITE(*,'(/,1x,A,F18.10)') 'Integrated Rho: ', sum(Rho) !
  ! Renormalize Rho (?)
  Rho(:) = Rho(:) / deltaV
END SUBROUTINE

