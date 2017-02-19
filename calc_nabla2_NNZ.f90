! Set the value of `nabla2_NNZ`
SUBROUTINE calc_nabla2_NNZ()
  USE m_LF3d, ONLY : lin2xyz => LF3d_lin2xyz, &
                     Npoints => LF3d_Npoints
  USE m_nabla2_sparse, ONLY : NNZ => nabla2_NNZ
  IMPLICIT NONE
  !
  LOGICAL :: Tnz(3)
  INTEGER :: ip1, ip2, i1,i2,j1,j2,k1,k2
  
  WRITE(*,'(/,1x,A)') 'Calculating NNZ'

  NNZ = 0
  DO ip2 = 1, Npoints
    DO ip1 = 1, Npoints ! use `ip1 = ip2, Npoints` if want to exploit symmetry
      Tnz(:) = .FALSE.
      i1 = lin2xyz(1,ip1)
      i2 = lin2xyz(1,ip2)
      !
      j1 = lin2xyz(2,ip1)
      j2 = lin2xyz(2,ip2)
      !
      k1 = lin2xyz(3,ip1)
      k2 = lin2xyz(3,ip2)
      !
      IF( j1 == j2 .AND. k1 == k2 ) Tnz(1) = .TRUE.
      IF( i1 == i2 .AND. k1 == k2 ) Tnz(2) = .TRUE.
      IF( i1 == i2 .AND. j1 == j2 ) Tnz(3) = .TRUE.
      !WRITE(*,*) ip1,ip2, any(Tnz)
      IF( ANY(Tnz) ) THEN 
        NNz = NNz + 1
        WRITE(*,'(1x,2I5,3L3)') ip1, ip2, Tnz(:)
      ENDIF
    ENDDO
    WRITE(*,'(1x,A,I5,A,I5)') 'rowIdx(', ip2+1, ')', NNZ + 1
  ENDDO

  WRITE(*,*) 'NNZ = ', NNZ
  WRITE(*,*) 'Actual size of matrix: ', Npoints**2
  WRITE(*,*) 'Percentage = ', ( dble(NNZ) ) / dble(Npoints)**2 * 100d0

END SUBROUTINE
