
! described in CPC-134-33 (2001)
SUBROUTINE minimE_pcg_Gan( Niter, restart )
  USE m_globals, ONLY : N, evecs, Ncol => Nstate
  IMPLICIT NONE
  !
  INTEGER :: Niter
  LOGICAL :: restart
  !
  REAL(8), ALLOCATABLE :: Yk(:,:), Fk(:,:), Bk(:,:), Gk(:,:)
  REAL(8), ALLOCATABLE :: Ak(:,:), Dk(:,:), Fk_old(:,:), Ak_old(:,:)
  REAL(8), ALLOCATABLE :: Pk(:,:)
  REAL(8), ALLOCATABLE :: gammak1(:), gammak2(:), gammak(:), gammak1_old(:), lambda_opt(:)
  INTEGER :: ic, iter
  REAL(8) :: Etot, Ekin, Epot, Etot_old
  !
!  REAL(8) :: ddot

  ALLOCATE( Yk(N**3,Ncol), Fk(N**3,Ncol), Bk(N**3,Ncol), Fk_old(N**3,Ncol) )
  ALLOCATE( Ak(N**3,Ncol), Dk(N**3,Ncol), Gk(N**3,Ncol), Ak_old(N**3,Ncol) )
  ALLOCATE( Pk(Ncol,Ncol) )
  ALLOCATE( gammak1(Ncol), gammak2(Ncol), gammak(Ncol), gammak1_old(Ncol), lambda_opt(Ncol) )

  ! initialize v
  IF( .NOT. restart ) THEN
    DO ic = 1, Ncol
      CALL r8_rand_vec( N**3, evecs(:,ic) )
    ENDDO
    CALL ortho_gram_schmidt( evecs, N**3, N**3, Ncol )
  ELSE
    READ(112) evecs
    ! No need to orthonormalize, v should already be ortonormal
  ENDIF

  CALL get_Etot( Ncol, evecs, Ekin, Epot, Etot)
  Etot_old = Etot

  DO iter = 1, Niter
    DO ic = 1, Ncol
      CALL apply_Ham( evecs(:,ic), Yk(:,ic) )
    ENDDO
    Pk = matmul( transpose(evecs), Yk )
    Fk = Yk - matmul( evecs, Pk )
    !
    !DO ic = 1, Ncol
    !  CALL apply_PrecCG( N**3, Fk(:,ic), Bk(:,ic), 200 )
    !ENDDO
    Bk = Fk  ! No preconditioning
    !
    Pk = matmul( transpose(evecs), Bk )
    !
    Gk = Bk - matmul( evecs, Pk )
    !
    Pk = matmul( transpose(Gk), Fk )
    gammak1 = 0.d0
    DO ic = 1, Ncol
      gammak1 = gammak1 + Pk(ic,ic)
    ENDDO
    !
    IF( iter > 1 ) THEN
      Pk = matmul( transpose(Gk), Fk_old )
      gammak2 = 0.d0
      DO ic = 1, Ncol
       gammak2 = gammak2 + Pk(ic,ic)
      ENDDO
    ELSE
      gammak2 = 0.d0
    ENDIF
    !
    IF( iter > 1 ) THEN
      gammak = ( gammak1 - gammak2 ) / gammak1_old
      DO ic = 1, Ncol
        Ak(:,ic) = -Gk(:,ic) + gammak(ic)*Ak_old(:,ic)
      ENDDO
    ELSE
      gammak = gammak1
      Ak = -Gk
    ENDIF
    Pk = matmul( transpose(evecs), Ak )
    Dk = Ak - matmul(evecs,Pk)
    !CALL linmin( Ncol, evecs, Dk, lambda_opt )
    lambda_opt = 1.d-4
    DO ic = 1, Ncol
      evecs(:,ic) = evecs(:,ic) + lambda_opt(ic) * Dk(:,ic)
    ENDDO
    !
    CALL ortho_gram_schmidt( evecs, N**3, N**3, Ncol )
    CALL get_Etot( Ncol, evecs, Ekin, Epot, Etot )
    WRITE(*,*)
    WRITE(*,'(1x,A,I8,2F18.10)') 'PCG G.A.N: iter, Etot, deltaE = ', iter, Etot, Etot_old-Etot
    !
    Etot_old = Etot
    gammak1_old = gammak1
    Fk_old = Fk
    Ak_old = Fk
  ENDDO

  
  DEALLOCATE( Yk, Fk, Bk, Ak, Dk, Fk_old, Ak_old )
  DEALLOCATE( Pk )
  DEALLOCATE( gammak1, gammak2, gammak, gammak1_old )


END SUBROUTINE



SUBROUTINE linmin( Ncol, Xk, Dk, lambda_opt )
  USE m_globals, ONLY : N
  IMPLICIT NONE
  INTEGER :: Ncol
  REAL(8) :: Xk(N**3,Ncol), Dk(N**3,Ncol)
  REAL(8) :: lambda_opt(Ncol)
  !
  REAL(8) :: a, b, c, discr, x1, x2
  REAL(8), ALLOCATABLE :: HX(:), HD(:)
  INTEGER :: ic
  REAL(8) :: XX, DHX, XHX, DX, DHD, DD
  REAL(8) :: XD, XHD
  !
  REAL(8) :: ddot

  ALLOCATE( Hx(N**3), HD(N**3) )

  DO ic = 1, Ncol
    CALL apply_Ham( Xk(:,ic), Hx(:) )
    CALL apply_Ham( Dk(:,ic), HD(:) )
    !
    XX = ddot( N**3, Xk(:,ic),1, Xk(:,ic),1 )
    DX = ddot( N**3, Dk(:,ic),1, Xk(:,ic),1 )
    XD = DX
    DD = ddot( N**3, Dk(:,ic),1, Dk(:,ic),1 )
    !
    DHX = ddot( N**3, Dk(:,ic),1, HX(:),1 )
    XHD = DHX
    XHX = ddot( N**3, Xk(:,ic),1, HX(:),1 )
    DHD = ddot( N**3, Dk(:,ic),1, HD(:),1 )

    a = XX * DHX - XHX * DX
    b = XX * DHD - XHX * DD
    c = XD * DHD - XHD * DD
    WRITE(*,*) 'linmin: a,b,c = ', a, b, c
    
    discr = b**2 - 4.d0*a*c
    IF( discr >= 0.d0 ) THEN
      x1 = ( -b**2 + sqrt( discr ) ) / (2.d0*c)
      x2 = ( -b**2 - sqrt( discr ) ) / (2.d0*c)
      lambda_opt(ic) = x1
      WRITE(*,'(1x,A,I5,2F18.10)') 'ic, lambda_opt = ', ic, x1,x2
    ELSE
      WRITE(*,*) 'Negative discriminant in linmin = ', discr
      lambda_opt(ic) = 1.d-5
    ENDIF

  ENDDO

  DEALLOCATE( HX, HD )
END SUBROUTINE
