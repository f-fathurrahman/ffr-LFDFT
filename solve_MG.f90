! efefer, 25 January 2016

! Using energy functional defined by Mauri and Galli
SUBROUTINE solve_MG()
  USE m_globals, ONLY : Nstate, N, evecs
  IMPLICIT NONE
  !
  REAL(8), ALLOCATABLE :: Hmn(:,:), Qmn(:,:), Smn(:,:)
  !REAL(8) :: eta
  REAL(8), ALLOCATABLE :: Hpsi(:)
  !
  INTEGER :: ic

  ALLOCATE( Hmn(Nstate,Nstate), Qmn(Nstate,Nstate), Smn(Nstate,Nstate) )
  ALLOCATE( Hpsi(N**3) )

  DO ic = 1, Nstate
    CALL r8_rand_vec( N**3, evecs(:,ic) )
  ENDDO

  CALL setupMatrix( )

  DEALLOCATE( Hmn, Qmn, Smn )
  DEALLOCATE( Hpsi )

  STOP


CONTAINS

  SUBROUTINE setupMatrix( )
    IMPLICIT NONE
    !
    INTEGER :: ic, jc
    REAL(8) :: ddot, deltaMN

    DO jc = 1, Nstate
      CALL apply_Ham( evecs(:,jc), Hpsi(:) )
      DO ic = jc, Nstate
        !
        Hmn(ic,jc) = ddot( N**3, evecs(:,ic),1, Hpsi(:),1 )
        Hmn(jc,ic) = Hmn(ic,jc)
        !
        Smn(ic,jc) = ddot( N**3, evecs(:,ic),1, evecs(:,jc),1 )
        Smn(jc,ic) = Smn(ic,jc)
        !
        Qmn(ic,jc) = 2.d0*deltaMN(ic,jc) - Smn(ic,jc)
        Qmn(jc,ic) = Qmn(ic,jc)
        !
        WRITE(*,*)
        WRITE(*,*) ic,jc
        WRITE(*,'(1x,A,3F18.10)') 'Hmn, Smn, Qmn = ', Hmn(ic,jc), Smn(ic,jc), Qmn(ic,jc)
      ENDDO
    ENDDO

  END SUBROUTINE


  SUBROUTINE calc_Etot_MG
  END SUBROUTINE

END SUBROUTINE
