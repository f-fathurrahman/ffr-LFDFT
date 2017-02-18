SUBROUTINE calc_rho()
  USE m_globals, ONLY : Nstate, evecs, Focc, rho, LF, N
  IMPLICIT NONE
  !
  INTEGER :: is
  REAL(8) :: ddot

  Rho(:) = 0.d0

  DO is = 1, Nstate
    WRITE(*,*) 'is, norm:', is, sqrt( ddot( N**3, evecs(:,is),1, evecs(:,is),1 ) )
    Rho(:) = Rho(:) + Focc(is) * evecs(:,is)*evecs(:,is)
  ENDDO
  WRITE(*,*) 'Integrated Rho: ', sum(Rho) !
END SUBROUTINE



SUBROUTINE init_rho_homogeneous()
  USE m_globals, ONLY : Rho, N, Nstate, Nelec
  IMPLICIT NONE

  WRITE(*,*) 'Generating initial Rho'
  Rho(:) = Nelec/N**3
  WRITE(*,*) 'Nelec/N**3 = ', Nelec/N**3
  WRITE(*,*) 'Integrated Rho:', sum(Rho)  ! no deltaV ?
END SUBROUTINE


SUBROUTINE init_rho_gaussian()
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : Rho, N, Nstate, LF
  IMPLICIT NONE
  INTEGER :: ip
  REAL(8) :: r, sigma
  REAL(8) :: deltaV, normChg

  deltaV = LF%LFx%h * LF%LFy%h * LF%LFz%h
  sigma = 0.75d0
  WRITE(*,*) 'Generating initial Rho (gaussian)'
  DO ip = 1, N**3
    r = norm2( LF%lingrid(:,ip) )
    !WRITE(*,'(1x,I4,4F7.3)') ip, LF%lingrid(:,ip), r
    Rho(ip) = exp(-r**2/(2*sigma**2))/(2*pi*sigma**2)**1.5d0
    !Rho(ip) = exp(-r**2/(sigma**2))/(sigma**3*sqrt(2d0*pi)**1.5)
  ENDDO
  !normChg = sum(Rho)
  !Rho = Rho/normChg
  WRITE(*,*) 'Integrated Rho:', sum(Rho)*deltaV
END SUBROUTINE


