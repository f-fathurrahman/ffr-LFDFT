SUBROUTINE avgPrec( vec, Kvec )
  USE m_globals, ONLY : N, LF
  IMPLICIT NONE
  !
  REAL(8) :: vec(N**3), Kvec(N**3)
  INTEGER :: ip, ipp, i, j, k, di, dj, dk
  REAL(8) :: avg
  INTEGER :: nn ! number of neighbour points

  DO ip = 1, N**3
    i = LF%lin2xyz(1,ip)
    j = LF%lin2xyz(2,ip)
    k = LF%lin2xyz(3,ip)
    avg = 0.d0
    ! For the case of Nx=Ny=Nz=N
    !
    nn = 0
    DO di=-1,1
    DO dj=-1,1
    DO dk=-1,1
      IF(i+di<=N .AND. i+di>=1 .AND. j+dj<=N .AND. j+dj>=1 .AND. k+dk<=N .AND. k+dk>=1 ) THEN
        ipp = LF%xyz2lin(i+di,j+dj,k+dk)
        avg = avg + vec(ipp)
        nn = nn + 1
      ENDIF
    ENDDO
    ENDDO
    ENDDO
    !WRITE(*,*) 'ip, nn', ip, nn
    Kvec(ip) = avg/nn
  ENDDO
  !STOP
END SUBROUTINE


