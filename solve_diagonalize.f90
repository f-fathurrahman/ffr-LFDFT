SUBROUTINE init_h_diag()
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : N, LF, h_diag, Vpsloc
  IMPLICIT NONE
  INTEGER :: ip, i, j, k
  REAL(8) :: L

  L = LF%LFx%L  ! special case for Lx=Ly=Lz=L

  ALLOCATE( h_diag(N**3) )
  DO ip=1,N**3
    i = LF%lin2xyz(1,ip)
    j = LF%lin2xyz(2,ip)
    k = LF%lin2xyz(3,ip)
    h_diag(ip) = Vpsloc(ip) + ( LF%LFx%D2jl(i,i) + LF%LFy%D2jl(j,j) + LF%LFz%D2jl(k,k) )*-0.5d0
    !Kx = i*PI/L
    !Ky = j*PI/L
    !Kz = k*PI/L
    !h_diag(ip) = 0.5d0*( Kx**2 + Ky**2 + Kz**2 ) + Vpsloc(ip)
    !h_diag(ip) = Vpsloc(ip)
  ENDDO
END SUBROUTINE




SUBROUTINE solve_diagonalize()
  USE m_globals, ONLY : Nstate, evals, evecs, Npoints
  IMPLICIT NONE
  INTEGER, ALLOCATABLE :: btype(:)
  INTEGER :: avg_iter, gstart, notcnv
  INTEGER :: lrot
  REAL(8) :: ethr
  INTEGER :: ist
  REAL(8) :: nrm
  REAL(8) :: ddot

  !CALL init_h_diag()

  WRITE(*,'(/,1x,A)') 'Solving Schrodinger equation via diagonalization:'
  WRITE(*,*)          '-------------------------------------------------'

  gstart = 2
  ethr = 1.d-6  ! FIXME this should not be setup dynamically during SCF iterations
  ALLOCATE( btype(Nstate) )
  btype(:) = 1 ! all bands are occupied
  lrot = .FALSE.
 
  WRITE(*,*) 'Generating random vectors ...'
  DO ist = 1, Nstate
    CALL r8_rand_vec( Npoints, evecs(:,ist) )
  ENDDO

  WRITE(*,*) 'Orthogonalizing initial vectors ...'
  CALL ortho_gram_schmidt( evecs, Npoints, Npoints, Nstate )

  !WRITE(*,*) 'Calling regterg ...'
  !CALL dbg_regterg( Npoints, Npoints, Nstate, 8*Nstate, evecs, ethr, &
  !                  gstart, evals, btype, notcnv, lrot, avg_iter )

  WRITE(*,*) 'Calling rcgdiagg ...'
  CALL rcgdiagg( Npoints, Npoints, Nstate, evecs, evals, btype, &
                 ethr, 1000, .true., notcnv, avg_iter )
                 

  WRITE(*,'(1x,A,I8,A)') 'Diagonalization converged in ', avg_iter, ' iterations.'
  WRITE(*,*) 'Eigenvalues:'
  DO ist = 1, Nstate
    WRITE(*,'(1x,I5,F18.10)') ist, evals(ist)
  ENDDO

  ! renormalize
  DO ist = 1, Nstate
    nrm = ddot( Npoints, evecs(:,ist),1, evecs(:,ist),1 )
    evecs(:,ist) = evecs(:,ist)/sqrt(nrm)
  ENDDO

  DEALLOCATE( btype )

END SUBROUTINE

