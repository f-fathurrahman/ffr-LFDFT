! FIXME: Not yet finalized. In this subroutine, we are supposed to test
! the implementation of get_grad() against finite difference result.
!-------------------------
SUBROUTINE test_grad(Ncol)
!-------------------------
  USE m_globals, ONLY : N
  IMPLICIT NONE
  !
  INTEGER :: Ncol
  REAL(8) :: v(N**3,Ncol)
  !
  INTEGER :: ic
  REAL(8) :: Etot, Ekin, Epot
  !
  REAL(8) :: ddot

  WRITE(*,*)
  WRITE(*,*) 'Initializing random vectors:'
  DO ic=1,Ncol
    CALL r8_rand_vec(N**3, v(:,ic))
    WRITE(*,*) ic, ddot( N**3, v(:,ic),1, v(:,ic),1 )
  ENDDO

  CALL ortho_gram_schmidt( v, N**3, N**3, Ncol)
  WRITE(*,*)
  WRITE(*,*) 'Norm:'
  DO ic=1,Ncol
    WRITE(*,*) ic, ddot( N**3, v(:,ic),1, v(:,ic),1 )
  ENDDO
  WRITE(*,*)
  WRITE(*,*) 'Norm wrt to col #1:'
  DO ic=2,Ncol
    WRITE(*,*) ic, ddot( N**3, v(:,ic),1, v(:,1),1 )
  ENDDO

  CALL get_Etot( Ncol, v, Ekin, Epot, Etot )
  WRITE(*,*) 'Ekin = ', Ekin
  WRITE(*,*) 'Epot = ', Epot
  WRITE(*,*) 'Etot = ', Etot

END SUBROUTINE
