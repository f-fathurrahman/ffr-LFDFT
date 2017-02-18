SUBROUTINE calc_rho()
  USE m_globals, ONLY : Nstate, evecs, Focc, rho, LF, N
  IMPLICIT NONE
  !
  INTEGER :: is
  REAL(8) :: ddot

  Rho(:) = 0.d0

  WRITE(*,'(/,1x,A)') 'Calculating rho (electron density):'
  WRITE(*,*)          '-----------------------------------'
  DO is = 1, Nstate
    WRITE(*,'(1x,A,I5,F18.10)') 'is, norm:', is, sqrt( ddot( N**3, evecs(:,is),1, evecs(:,is),1 ) )
    Rho(:) = Rho(:) + Focc(is) * evecs(:,is)*evecs(:,is)
  ENDDO
  WRITE(*,'(/,1x,A,F18.10)') 'Integrated Rho: ', sum(Rho) !
END SUBROUTINE

