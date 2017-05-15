SUBROUTINE ortho_check( Npoints, Ncols, dVol, v )
  IMPLICIT NONE 
  INTEGER :: Npoints
  INTEGER :: Ncols
  REAL(8) :: dVol
  REAL(8) :: v(Npoints,Ncols)
  !
  INTEGER :: ic
  REAL(8) :: nrm
  !
  REAL(8) :: ddot

  WRITE(*,*)
  WRITE(*,*) 'Checking orthonormalization'
  WRITE(*,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^'
  WRITE(*,*)
  WRITE(*,*) 'Norms:'
  DO ic = 1, Ncols
    nrm = ddot( Npoints, v(:,ic), 1, v(:,ic),1 ) * dVol
    WRITE(*,'(1x,I8,F18.10)') ic, nrm
  ENDDO
  
  WRITE(*,*)
  WRITE(*,*) 'Check wrt to vector 1'
  DO ic = 2,Ncols
    nrm = ddot( Npoints, v(:,ic), 1, v(:,1),1 ) * dVol
    WRITE(*,'(1x,I8,F18.10)') ic, nrm
  ENDDO
  WRITE(*,*)

END SUBROUTINE 
