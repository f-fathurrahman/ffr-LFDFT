SUBROUTINE init_h_diag()
  USE m_constants, ONLY : PI
  USE m_globals, ONLY : N, LF, h_diag, Vpot
  IMPLICIT NONE
  INTEGER :: ip, i, j, k
  REAL(8) :: Kx,Ky,Kz,L

  L = LF%LFx%L  ! special case for Lx=Ly=Lz=L

  ALLOCATE( h_diag(N**3) )
  DO ip=1,N**3
    i = LF%lin2xyz(1,ip)
    j = LF%lin2xyz(2,ip)
    k = LF%lin2xyz(3,ip)
    h_diag(ip) = Vpot(ip) + ( LF%LFx%D2jl(i,i) + LF%LFy%D2jl(j,j) + LF%LFz%D2jl(k,k) )*-0.5d0
    !Kx = i*PI/L
    !Ky = j*PI/L
    !Kz = k*PI/L
    !h_diag(ip) = 0.5d0*( Kx**2 + Ky**2 + Kz**2 ) + Vpot(ip)
    !h_diag(ip) = Vpot(ip)
  ENDDO
END SUBROUTINE




SUBROUTINE solve_diagonalize()
  USE m_globals
  IMPLICIT NONE
  INTEGER, ALLOCATABLE :: btype(:)
  INTEGER :: dav_iter, gstart, notcnv
  INTEGER :: lrot
  REAL(8) :: ethr
  INTEGER :: ist
  REAL(8) :: nrm
  REAL(8) :: ddot

  !CALL init_h_diag()

  gstart = 2
  ethr = 1.d-4
  ALLOCATE( btype(Nstate) )
  btype(:) = 1 ! all bands are occupied
  lrot = .FALSE.
 
  DO ist = 1, Nstate
    CALL r8_rand_vec( N**3, evecs(:,ist) )
  ENDDO
  CALL ortho_gram_schmidt( evecs, N**3, N**3, Nstate )

  CALL regterg( N**3, N**3, Nstate, 5*Nstate, evecs, ethr, &
                    gstart, evals, btype, notcnv, lrot, dav_iter )
  WRITE(*,*) 'dav_iter = ', dav_iter
  DO ist = 1, Nstate
    WRITE(*,'(1x,I5,F18.10)') ist, evals(ist)
  ENDDO

  ! renormalize
  DO ist = 1, Nstate
    nrm = ddot( N**3, evecs(:,ist),1, evecs(:,ist),1 )
    evecs(:,ist) = evecs(:,ist)/sqrt(nrm)
  ENDDO

  DEALLOCATE( btype )

END SUBROUTINE

