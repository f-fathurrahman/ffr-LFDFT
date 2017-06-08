SUBROUTINE gen_random_evecs()
  USE m_states, ONLY : Nstates, evecs => KS_evecs
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  IMPLICIT NONE 
  INTEGER :: ist, ip
  ! Initialize to random wavefunction
  DO ist = 1, Nstates
    DO ip = 1, Npoints
      CALL random_number( evecs(ip,ist) )
    ENDDO
  ENDDO
  CALL orthonormalize( Nstates, evecs )
END SUBROUTINE 

