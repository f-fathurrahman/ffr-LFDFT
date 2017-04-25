SUBROUTINE diag_davidson( evals, v, TOLERANCE )
  
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_states, ONLY : Nstates
  IMPLICIT NONE
  ! Arguments
  REAL(8) :: v(Npoints,Nstates)
  REAL(8) :: evals(Nstates)
  ! Local variable
  REAL(8) :: RNORM
  INTEGER :: ist, istEP, MAX_DIR, NCONV
  LOGICAL :: IS_CONVERGED
  REAL(8) :: MACHINE_ZERO, TOLERANCE
  REAL(8), ALLOCATABLE :: RES_TOL(:), RES_NORM(:), evals_red(:)
  REAL(8), ALLOCATABLE :: cmat(:,:), H_MAT(:,:), O_MAT(:,:), evecs(:,:)
  REAL(8), ALLOCATABLE :: HV(:,:), R(:,:), HR(:,:), xtemp(:,:)
  REAL(8), ALLOCATABLE :: evals_old(:)
  ! BLAS function
  REAL(8) :: ddot

  ALLOCATE( RES_TOL(Nstates) )
  ALLOCATE( RES_NORM(Nstates) )
  ALLOCATE( cmat(Nstates,Nstates) )

  ALLOCATE( H_MAT(2*Nstates,2*Nstates) )
  ALLOCATE( O_MAT(2*Nstates,2*Nstates) )
  ALLOCATE( evecs(2*Nstates,2*Nstates) )
  ALLOCATE( evals_red(2*Nstates) )
  ALLOCATE( evals_old(Nstates) )

  ALLOCATE( HV(Npoints,Nstates) )
  ALLOCATE( R(Npoints,Nstates) )
  ALLOCATE( HR(Npoints,Nstates) )
  ALLOCATE( xtemp(Npoints,Nstates) )

  !CALL ortho_gram_schmidt( V, Npoints, Npoints, Nstates )
  V(:,:) = sqrt(dVol)*V(:,:)

  ! Apply Hamiltonian
  CALL op_H( Nstates, V, HV ) 

  ! Calculate Rayleigh quotient
  DO ist=1,Nstates
    evals(ist) = ddot(Npoints, V(:,ist),1, HV(:,ist),1)
  ENDDO

  ! Calculate matrix of residual vector
  DO ist=1,Nstates
    R(:,ist) = evals(ist)*V(:,ist) - HV(:,ist)
    RES_TOL(ist) = SQRT( ddot(Npoints, R(:,ist),1, R(:,ist),1) )
  ENDDO

  istep = 1
  IS_CONVERGED = .FALSE.
  MAX_DIR = 100
  MACHINE_ZERO = 2.220446049250313D-16
  
  RNORM = 1.D0

  evals_old(:) = 0.d0

  DO WHILE ( (istep <= MAX_DIR) .AND. (.NOT.IS_CONVERGED) )
    
    WRITE(*,'(I8,ES18.10)') istep, RNORM
    RES_NORM = 1.D0

    !WHERE(MACHINE_ZERO < RES_TOL) RES_NORM = 1.d0/RES_TOL
    !WRITE(*,*) 'RES_NORM:', RES_NORM
    
    DO ist = 1,Nstates
      IF(MACHINE_ZERO < RES_TOL(ist)) RES_NORM(ist) = 1.D0/RES_TOL(ist)
      !WRITE(*,*) ist, RES_NORM(ist)
    END DO

    ! Scale the residual vectors
    DO ist=1,Nstates
      R(:,ist) = RES_NORM(ist)*R(:,ist)
    ENDDO

    ! Apply preconditioner
    DO ist=1,Nstates
      call prec_ilu0_inplace( R(:,ist) )
    ENDDO

! Construct the reduced hamiltonian. The reduced hamiltonian has dimensions
!  2nb x 2nb and is constructed by filling in four nb x nb blocks one at a time:
! __ 
!|  |
!| <v|H|v>   <v|H|r>  |
!    h_mat = |  |
!| *******   <r|H|r>  | 
!|__|

    CALL op_H( Nstates, R, HR )

    IF(istep == 1) THEN
      CALL DGEMM('T','N',Nstates,Nstates,Npoints, 1.d0, V, Npoints, HV,Npoints, 0.d0,cmat,Nstates)
      H_MAT(1:Nstates,1:Nstates) = cmat
    ELSE
      H_MAT(1:Nstates,1:Nstates) = 0.d0
      DO ist = 1,Nstates
        H_MAT(ist,ist) = evals(ist)
      ENDDO
    ENDIF

    ! <v|H|r> --> cmat
    CALL DGEMM('T','N',Nstates,Nstates,Npoints, 1.d0,V,Npoints, HR,Npoints, 0.d0,cmat,Nstates)
    H_MAT(1:Nstates,Nstates+1:2*Nstates) = cmat
    CALL DGEMM('T','N',Nstates,Nstates,Npoints, 1.d0,HR,Npoints, V,Npoints, 0.d0,cmat,Nstates)
    H_MAT(Nstates+1:2*Nstates,1:Nstates) = cmat
    ! <r|H|r> --> cmat
    CALL DGEMM('T','N',Nstates,Nstates,Npoints, 1.d0,R,Npoints, HR,Npoints, 0.d0,cmat,Nstates)
    H_MAT(Nstates+1:2*Nstates,Nstates+1:2*Nstates) = cmat


! Construct the reduced overlap matrix which has dimenstions 2nb x 2nb
!   and is constructed by filling in four nb x nb blocks one at a time:
! _   _ 
!|     |
!|  <v|v>   <v|r>  |
!    o_mat = |     |
!|  *****   <r|r>  | 
!|_   _|

    O_MAT(1:Nstates,1:Nstates) = 0.d0
    DO ist = 1,Nstates
      O_MAT(ist,ist) = 1.d0
    END DO
    ! <v|r> --> cmat
    CALL DGEMM('T','N',Nstates,Nstates,Npoints, 1.d0,V,Npoints, R,Npoints, 0.d0,cmat,Nstates)
    O_MAT(1:Nstates,Nstates+1:2*Nstates) = cmat
    CALL DGEMM('T','N',Nstates,Nstates,Npoints, 1.d0,R,Npoints, V,Npoints, 0.d0,cmat,Nstates)
    O_MAT(Nstates+1:2*Nstates,1:Nstates) = cmat
    ! <r|r> --> cmat
    CALL DGEMM('T','N',Nstates,Nstates,Npoints, 1.d0,R,Npoints, R,Npoints, 0.d0,cmat,Nstates)
    O_MAT(Nstates+1:2*Nstates,Nstates+1:2*Nstates) = cmat

    CALL rdiaghg( 2*Nstates, 2*Nstates, H_MAT, O_MAT, 2*Nstates, evals_red, evecs )

    evals = evals_red(1:Nstates)
    cmat = evecs(1:Nstates,1:Nstates)
    ! V*cmat --> V
    CALL DGEMM('N','N',Npoints,Nstates,Nstates, 1.d0,V,Npoints, cmat,Nstates, 0.d0,xtemp,Npoints)
    V = xtemp
    ! HV = HV*cmat
    CALL DGEMM('N','N',Npoints,Nstates,Nstates, 1.d0,HV,Npoints, cmat,Nstates, 0.d0,xtemp,Npoints)
    HV = xtemp
    !
    cmat = evecs(Nstates+1:2*Nstates,1:Nstates)
    ! V = V + R*cmat
    CALL DGEMM('N','N',Npoints,Nstates,Nstates, 1.d0,R,Npoints, cmat,Nstates, 1.d0,V,Npoints)
    ! HV = HV + HR*cmat
    CALL DGEMM('N','N',Npoints,Nstates,Nstates, 1.d0,HR,Npoints, cmat,Nstates, 1.d0,HV,Npoints)

    ! Calculate matrix of residual vector
    DO ist=1,Nstates
      R(:,ist) = evals(ist)*V(:,ist) - HV(:,ist)
      RES_TOL(ist) = SQRT( ddot(Npoints, R(:,ist),1, R(:,ist),1) )
    ENDDO

    WRITE(*,*)
    WRITE(*,*) 'Eigenvalues convergence:'
    DO ist = 1,Nstates
      WRITE(*,'(1X,I5,F18.10,ES18.10)') ist, evals(ist), abs( evals(ist)-evals_old(ist) )
    ENDDO 

    RNORM = SUM( abs(evals - evals_old) )/REAL(Nstates, kind=8)

    IS_CONVERGED = rnorm <= TOLERANCE

    evals_old(:) = evals(:)
    istep = istep + 1
  END DO

  WRITE(*,*) 'End od block-Davidson iteration: rnorm = ', RNORM
  
  DEALLOCATE(RES_TOL)
  DEALLOCATE(RES_NORM)
  DEALLOCATE(cmat)
  DEALLOCATE(H_MAT)
  DEALLOCATE(O_MAT)
  DEALLOCATE(evecs)
  DEALLOCATE(evals_red)
  DEALLOCATE(HV)
  DEALLOCATE(R)
  DEALLOCATE(HR)
  DEALLOCATE(xtemp)
END SUBROUTINE


